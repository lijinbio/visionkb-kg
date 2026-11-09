"""KGX -> neo4j-admin import CSVs (contracts §5): nodes_<Label>.csv, rels_<type>.csv, import.args, post_import.cypher."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

from .kgx import read_jsonl

ARRAY_DELIM = ";"
FULLTEXT_EXCLUDE = {"StudyResult", "NamedThing"}


def pascal(name: str) -> str:
    return "".join(w[:1].upper() + w[1:] for w in re.split(r"[ _]+", name) if w)


def local(curie: str) -> str:
    return curie.split(":", 1)[1]


class LabelMaker:
    def __init__(self, tk):
        self.tk = tk
        self.cache = {}

    def labels(self, categories: list[str]) -> list[str]:
        key = tuple(categories)
        if key in self.cache:
            return self.cache[key]
        out = []
        for c in categories:
            if c.startswith("visionkb:"):
                out.append(local(c))
        for c in categories:
            if not c.startswith("biolink:"):
                continue
            for a in self.tk.get_ancestors(c, reflexive=True, mixin=False, formatted=False):
                if a == "entity":
                    continue
                lab = pascal(a)
                if lab not in out:
                    out.append(lab)
        if "NamedThing" in out:
            out.remove("NamedThing")
        out.append("NamedThing")
        self.cache[key] = out
        return out


def _scalar(v) -> str:
    s = str(v)
    return s.replace("\r", " ").replace("\n", " ")


def _type(values) -> str:
    kinds = set()
    for v in values:
        if v is None:
            continue
        if isinstance(v, bool):
            kinds.add("boolean")
        elif isinstance(v, int):
            kinds.add("long")
        elif isinstance(v, float):
            kinds.add("double")
        elif isinstance(v, list):
            kinds.add("string[]")
        else:
            kinds.add("string")
    if kinds <= {"long"}:
        return "long"
    if kinds <= {"long", "double"}:
        return "double"
    if kinds == {"boolean"}:
        return "boolean"
    if "string[]" in kinds:
        return "string[]"
    return "string"


def _cell(v, typ: str) -> str:
    if v is None:
        return ""
    if typ == "string[]":
        vals = v if isinstance(v, list) else [v]
        return ARRAY_DELIM.join(_scalar(x).replace(ARRAY_DELIM, ",") for x in vals)
    if typ == "boolean":
        return "true" if v else "false"
    return _scalar(v)


def _write(path: Path, header: list[str], rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def write_import(nodes_path: Path, edges_path: Path, outdir: Path, tk) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    lm = LabelMaker(tk)
    groups: dict[str, list] = defaultdict(list)
    sds = []
    for n in read_jsonl(nodes_path):
        labs = lm.labels(n["category"])
        groups[labs[0]].append((labs, n))
        if "biolink:StudyResult" in n["category"]:
            for d in n.get("supporting_data_set") or []:
                sds.append((n["id"], d))
    files = {"nodes": [], "relationships": []}
    counts = {"nodes": {}, "relationships": {}}
    for lab in sorted(groups):
        recs = groups[lab]
        keys = sorted({k for _, n in recs for k in n if k not in ("id",)})
        types = {k: _type(n.get(k) for _, n in recs) for k in keys}
        types["category"] = "string[]"
        header = ["id:ID", ":LABEL"] + [f"{k}:{types[k]}" for k in keys]
        rows = (
            [n["id"], ARRAY_DELIM.join(labs)] + [_cell(n.get(k), types[k]) for k in keys]
            for labs, n in sorted(recs, key=lambda x: x[1]["id"])
        )
        fn = f"nodes_{lab}.csv"
        _write(outdir / fn, header, rows)
        files["nodes"].append(fn)
        counts["nodes"][lab] = len(recs)
    del groups
    rels: dict[str, list] = defaultdict(list)
    for e in read_jsonl(edges_path):
        rels[local(e["predicate"])].append(e)
    for t, recs in sorted(rels.items()):
        keys = sorted({k for e in recs for k in e if k not in ("subject", "object", "id")})
        types = {k: _type(e.get(k) for e in recs) for k in keys}
        header = [":START_ID", ":END_ID", ":TYPE", "id"] + [f"{k}:{types[k]}" for k in keys]
        rows = (
            [e["subject"], e["object"], t, e["id"]] + [_cell(e.get(k), types[k]) for k in keys]
            for e in sorted(recs, key=lambda x: x["id"])
        )
        fn = f"rels_{t}.csv"
        _write(outdir / fn, header, rows)
        files["relationships"].append(fn)
        counts["relationships"][t] = len(recs)
    if sds:
        fn = "rels_supporting_data_set.csv"
        _write(
            outdir / fn,
            [":START_ID", ":END_ID", ":TYPE", "predicate"],
            ([s, d, "supporting_data_set", "biolink:supporting_data_set"] for s, d in sorted(sds)),
        )
        files["relationships"].append(fn)
        counts["relationships"]["supporting_data_set"] = len(sds)
    args = ["--id-type=string", f"--array-delimiter={ARRAY_DELIM}", "--skip-duplicate-nodes=false", "--bad-tolerance=0"]
    # "./" prefix: neo4j-admin 2026.10 rejects bare file names ("Unable to find the parent of the path"); run from neo4j-import/
    args += [f"--nodes=./{f}" for f in files["nodes"]] + [f"--relationships=./{f}" for f in files["relationships"]]
    (outdir / "import.args").write_text("\n".join(args) + "\n")
    labels_ft = sorted(set(counts["nodes"]) - FULLTEXT_EXCLUDE)
    cy = [
        "// VisionKB post-import schema (contracts §5). Run once after neo4j-admin database import full.",
        "CREATE CONSTRAINT namedthing_id IF NOT EXISTS FOR (n:NamedThing) REQUIRE n.id IS UNIQUE;",
        "CREATE INDEX studyresult_edge_id IF NOT EXISTS FOR (n:StudyResult) ON (n.edge_id);",
        "CREATE INDEX studyresult_dataset_id IF NOT EXISTS FOR (n:StudyResult) ON (n.dataset_id);",
        "CREATE INDEX gene_symbol IF NOT EXISTS FOR (n:Gene) ON (n.symbol);",
        f"CREATE FULLTEXT INDEX entity_search IF NOT EXISTS FOR (n:{'|'.join(labels_ft)}) ON EACH [n.name, n.symbol, n.synonym];",
    ]
    for t in sorted(counts["relationships"]):
        if t != "supporting_data_set":
            cy.append(f"CREATE INDEX rel_id_{t} IF NOT EXISTS FOR ()-[r:{t}]-() ON (r.id);")
    cy.append("CALL db.awaitIndexes(3600);")
    (outdir / "post_import.cypher").write_text("\n".join(cy) + "\n")
    return {"files": files, "counts": counts}
