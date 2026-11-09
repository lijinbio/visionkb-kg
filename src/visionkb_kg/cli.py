"""visionkb-kg command line: ingest <source>, build, pack, verify, all."""

from __future__ import annotations

import argparse
import datetime
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import yaml

from . import BIOLINK_VERSION, SCHEMA_VERSION, __version__
from .kgx import read_jsonl, write_source
from .sources import INFORES, Sources, build_dir, data_root, dirname_of, release_dir

INGESTS = ("hgnc", "heca", "retigene", "datasets", "ontologies")
ORDER = ("ontologies", "hgnc", "heca", "retigene", "datasets")  # merge priority: first source wins for scalars


def log(msg: str):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def wjson(obj, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str) + "\n")


def _atlases(src: Sources) -> tuple[list[dict], dict]:
    """Registry records of the atlases in this release, and per-dataset context from atlases.yaml."""
    cfg = yaml.safe_load(Path(src.get("atlases_config")).read_text())
    version = cfg.get("atlas_version")
    root = Path(src.get("atlas_root"))
    regs, ctx = [], {}
    for a in cfg["atlases"]:
        v = a.get("atlas_version", version)
        d = root / dirname_of(a["dataset_id"]) / v
        reg = json.loads((d / "registry.json").read_text())
        reg["_dir"] = str(d)
        regs.append(reg)
        ctx[a["dataset_id"]] = {
            "donor_col": a.get("donor_col"),
            "tissue": a["tissue"],
            "region_uberon": a.get("region_uberon") or {},
            "access_level": a.get("access_level", "lab"),
        }
    return regs, ctx


def _hgnc_df(src: Sources):
    from .ingest import hgnc

    return hgnc.load(src.get("hgnc", "path"))


def ingest(name: str, release: str, src: Sources) -> dict:
    out = build_dir(release) / "kgx" / name
    rep_dir = build_dir(release) / "reports"
    if name == "heca":
        import pandas as pd

        from .ingest import heca

        regs, ctx = _atlases(src)
        files = [Path(r["_dir"]) / "summaries" / "expr_summary.parquet" for r in regs]
        groups = pd.concat(
            [pd.read_parquet(Path(r["_dir"]) / "summaries" / "groups.parquet") for r in regs], ignore_index=True
        )
        allrows = heca.load_levels(files)
        donor_ok = {d: bool(a.get("donor_col")) for d, a in ctx.items()}
        dcp = src.get("atlas_group_donor_cells")
        donor_cells = pd.read_parquet(dcp) if dcp else None
        rows, d18 = heca.select_d18(allrows, donor_ok, donor_cells=donor_cells)
        del allrows
        log(
            f"heca: {d18['totals']['kept']:,} of {d18['totals']['candidate_rows']:,} gene x group rows kept by {heca.RULE['decision']} ({len(files)} atlases)"
        )
        nodes, edges, stats = heca.build(rows, groups, ctx, release, heca.method_string())
        stats["edge_rule"] = d18
    elif name == "retigene":
        import pandas as pd

        from .ingest import hgnc, retigene

        r = src.get("retigene")
        d = Path(r["dir"])
        genes = pd.read_csv(d / r["genes"])
        ev = pd.read_csv(d / r["evidence"])
        queries = yaml.safe_load(Path(r["phenotype_queries"]).read_text())["genes"]
        mapping = retigene.map_phenotypes(genes, queries)
        rep_dir.mkdir(parents=True, exist_ok=True)
        mapping.to_csv(rep_dir / "retigene_phenotype_mapping.tsv", sep="\t", index=False)
        nodes, edges, stats = retigene.build(genes, ev, mapping, hgnc.symbol_index(_hgnc_df(src)), release)
        stats["mapping"] = mapping["status"].value_counts().to_dict()
    elif name == "datasets":
        from .ingest import datasets

        regs, _ = _atlases(src)
        nodes = [datasets.node({k: v for k, v in r.items() if k != "_dir"}, release) for r in regs]
        edges = []
        stats = {"datasets": len(nodes)}
    elif name == "hgnc":
        import pyarrow.dataset as pds

        from .ingest import hgnc

        extra = set()
        for s in ("heca", "retigene"):
            for e in read_jsonl(build_dir(release) / "kgx" / s / "edges.jsonl"):
                for v in (e["subject"], e["object"]):
                    if v.startswith("HGNC:"):
                        extra.add(v)
        regs, _ = _atlases(src)
        for r in regs:
            t = pds.dataset(str(Path(r["_dir"]) / "summaries" / "expr_summary.parquet")).to_table(columns=["gene_id"])
            extra |= {g for g in set(t.column("gene_id").to_pylist()) if g.startswith("HGNC:")}
        nodes, stats = hgnc.nodes(_hgnc_df(src), extra)
        edges = []
    elif name == "ontologies":
        from .ingest import ontologies

        ref = set()
        for s in ("heca", "retigene", "datasets", "hgnc"):
            p = build_dir(release) / "kgx" / s
            for n in read_jsonl(p / "nodes.jsonl"):
                for k in ("parent_cell_type", "in_taxon", "species", "assay"):
                    if isinstance(n.get(k), str):
                        ref.add(n[k])
                for k in ("tissue", "cell_types"):
                    ref |= set(n.get(k) or [])
            for e in read_jsonl(p / "edges.jsonl"):
                ref |= {e["subject"], e["object"], *(e.get("anatomical_context_qualifier") or [])}
                if e.get("species_context_qualifier"):
                    ref.add(e["species_context_qualifier"])
        ref = {x for x in ref if x.split(":")[0] in ("CL", "UBERON", "MONDO", "HP", "EFO", "NCBITaxon")}
        obo = ontologies.load_obo_sets({k: src.get("ontologies", k, "path") for k in ("cl", "uberon", "mondo", "hp")})
        nodes, edges, stats = ontologies.slim(ref, obo, release)
    else:
        raise SystemExit(f"unknown source {name}")
    stats.update(write_source(out, nodes, edges))
    wjson(stats, rep_dir / f"ingest_{name}.json")
    log(f"{name}: {stats['nodes']:,} nodes, {stats['edges']:,} edges")
    return stats


def _git(path: str | None) -> str | None:
    if not path:
        return None
    try:
        return subprocess.run(
            ["git", "-C", path, "rev-parse", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def term_sets(src: Sources):
    from .ingest import hgnc as hg
    from .obo import read_obo
    from .validate import TermSets

    df = _hgnc_df(src)
    a = df[df.status == "Approved"]
    obo = {}
    for key, pre in (("cl", "CL"), ("uberon", "UBERON"), ("mondo", "MONDO"), ("hp", "HP")):
        _, terms = read_obo(src.get("ontologies", key, "path"), prefixes=(pre,))
        obo[pre] = {t for t, r in terms.items() if not r["obsolete"]}
    del hg
    return TermSets(dict(zip(a.hgnc_id, a.symbol, strict=True)), obo)


def build(release: str, src: Sources, *, run_kgx: bool = True) -> Path:
    import pyarrow as pa
    import pyarrow.parquet as pq

    from .biolink import biolink_sha256, toolkit
    from .neo4j_csv import write_import
    from .release import file_list, merge, release_counts, sha256_file
    from .validate import apln_self_test, validate_release

    final = release_dir(release)
    if final.exists():
        raise SystemExit(f"{final} exists; releases are immutable")
    b = build_dir(release)
    stage = b / "release" / release
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    (stage / "RELEASE").write_text(release + "\n")
    log("merge")
    mstats = merge([b / "kgx" / s for s in ORDER], stage / "kgx")
    wjson(mstats, b / "reports" / "merge.json")
    log("validate")
    ts = term_sets(src)
    rep = validate_release(stage / "kgx" / "nodes.jsonl", stage / "kgx" / "edges.jsonl", ts, run_kgx=run_kgx, log=log)
    rep["self_tests"] = {"apln_wrong_id": apln_self_test(ts)}
    rep["ingest"] = {s: json.loads((b / "reports" / f"ingest_{s}.json").read_text()) for s in INGESTS}
    rep["release"] = release
    if rep["self_tests"]["apln_wrong_id"]["status"] != "pass":
        rep["status"] = "fail"
    wjson(rep, stage / "validation" / "report.json")
    (stage / "validation" / "report.md").write_text(report_md(rep))
    shutil.copy(src.crosswalk, stage / "validation" / "heca_celltype_crosswalk.sssom.tsv")
    if (b / "reports" / "retigene_phenotype_mapping.tsv").exists():
        shutil.copy(
            b / "reports" / "retigene_phenotype_mapping.tsv", stage / "validation" / "retigene_phenotype_mapping.tsv"
        )
    log(f"validation {rep['status']}: {rep['n_errors']} errors, {rep['n_warnings']} warnings")
    log("neo4j-import")
    neo = write_import(stage / "kgx" / "nodes.jsonl", stage / "kgx" / "edges.jsonl", stage / "neo4j-import", toolkit())
    log("summaries + registry")
    regs, _ = _atlases(src)
    tables = []
    atl = []
    for r in sorted(regs, key=lambda x: x["id"]):
        f = Path(r["_dir"]) / "summaries" / "expr_summary.parquet"
        tables.append(pq.read_table(f))
        atl.append({"dataset_id": r["id"], "atlas_version": r["version"], "summary_sha256": sha256_file(f)})
    (stage / "summaries").mkdir()
    # atlases without a region column store `region` as an all-null column: promote null -> string
    merged = pa.concat_tables(tables, promote_options="permissive")
    pq.write_table(merged, stage / "summaries" / "expr_summary.parquet", compression="zstd")
    del tables
    log("known issues")
    ki_entries = known_issues_step(stage, regs, rep)
    from .ingest.datasets import registry_line

    (stage / "registry").mkdir()
    with open(stage / "registry" / "datasets.jsonl", "w") as f:
        for r in sorted(regs, key=lambda x: x["id"]):
            f.write(
                json.dumps(registry_line({k: v for k, v in r.items() if k != "_dir"}, release), sort_keys=True) + "\n"
            )
    srcs = []
    h = src.get("hgnc")
    srcs.append(
        {
            "name": "HGNC complete set",
            "version": str(h.get("date")),
            "path_or_url": h["url"],
            "sha256": sha256_file(Path(h["path"])),
            "license": h.get("license"),
            "infores": INFORES["hgnc"],
        }
    )
    vers = rep["ingest"]["ontologies"]["versions"]
    for key, pre in (
        ("cl", "CL"),
        ("uberon", "UBERON"),
        ("mondo", "MONDO"),
        ("hp", "HP"),
        ("efo", "EFO"),
        ("ncbitaxon", "NCBITaxon"),
    ):
        o = src.get("ontologies", key) or {}
        rec = {
            "name": pre,
            "version": vers.get(pre),
            "path_or_url": o.get("url"),
            "license": o.get("license"),
            "infores": INFORES[key],
        }
        if o.get("path"):
            rec["sha256"] = sha256_file(Path(o["path"]))
        srcs.append(rec)
    rg = src.get("retigene")
    srcs.append(
        {
            "name": "retigene curated haploinsufficiency review",
            "version": rg.get("version"),
            "path_or_url": "lab-internal",
            "sha256": sha256_file(Path(rg["dir"]) / rg["genes"]),
            "license": rg.get("license"),
            "infores": INFORES["retigene"],
        }
    )
    srcs.append(
        {
            "name": "HECA atlas summaries",
            "version": sorted({r["version"] for r in regs}),
            "path_or_url": "summaries/expr_summary.parquet",
            "license": "lab-internal",
            "infores": INFORES["heca"],
        }
    )
    srcs.append(
        {
            "name": "HECA label to CL crosswalk",
            "version": "0.1",
            "path_or_url": "validation/heca_celltype_crosswalk.sssom.tsv",
            "sha256": sha256_file(src.crosswalk),
            "license": "CC-BY-4.0",
            "infores": INFORES["heca"],
        }
    )
    manifest = {
        "kg_release": release,
        "created_utc": datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "biolink_version": BIOLINK_VERSION,
        "biolink_schema_sha256": biolink_sha256(),
        "schema_version": SCHEMA_VERSION,
        "code": {
            "visionkb-kg": {"version": __version__, "git_sha": _git(src.get("code", "visionkb-kg"))},
            "visionkb-atlas": {
                "version": src.get("code", "visionkb-atlas-version"),
                "git_sha": _git(src.get("code", "visionkb-atlas")),
            },
        },
        "sources": srcs,
        "atlases": atl,
        "counts": {**release_counts(stage), "neo4j": neo["counts"]},
        "validation": {"status": rep["status"], "n_errors": rep["n_errors"], "n_warnings": rep["n_warnings"]},
        "expression_edge_rule": {
            **rep["ingest"]["heca"]["edge_rule"]["rule"],
            "donor_support_checked": rep["ingest"]["heca"]["edge_rule"]["donor_support_checked"],
            "robustness_specificity_fraction": rep["ingest"]["heca"]["edge_rule"]["robustness"],
        },
        "known_issues": ki_entries,
        "neo4j_import": {
            "dir": "neo4j-import",
            "args_file": "neo4j-import/import.args",
            "post_import": "neo4j-import/post_import.cypher",
            "run_from": "neo4j-import",
        },
        "files": None,
    }
    manifest["files"] = file_list(stage)
    wjson(manifest, stage / "manifest.json")
    final.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(stage), str(final))
    log(f"release written: {final}")
    return final


def known_issues_step(stage: Path, regs: list[dict], rep: dict) -> list[dict]:
    """Compute the known issues (report only) from the release summaries; write validation/known_issues.json, append the
    'Known issues' section to validation/report.md, store a summary in report.json and return the manifest entries."""
    import pandas as pd

    from . import known_issues as KI

    groups = pd.concat(
        [pd.read_parquet(Path(r["_dir"]) / "summaries" / "groups.parquet") for r in regs], ignore_index=True
    )
    donor_ok = rep["ingest"]["heca"]["edge_rule"]["donor_support_checked"]
    ki = KI.cross_lineage(stage / "summaries" / "expr_summary.parquet", groups, donor_ok)
    wjson(ki, stage / "validation" / "known_issues.json")
    md = stage / "validation" / "report.md"
    md.write_text(md.read_text().split("\n## Known issues\n")[0].rstrip("\n") + "\n" + KI.markdown(ki))
    rp = stage / "validation" / "report.json"
    r = json.loads(rp.read_text())
    r["known_issues"] = {k: ki[k] for k in ("id", "title", "n_pairs", "n_groups", "n_study_results")}
    wjson(r, rp)
    rod = [p for p in ki["pairs"] if p["program_of"] == "Rod"]
    return [
        {
            "id": ki["id"],
            "title": ki["title"],
            "status": "open; report only, no edges flagged or removed (decision D19); contamination QC planned",
            "section": "validation/report.md#known-issues",
            "details": "validation/known_issues.json",
            "n_pairs": ki["n_pairs"],
            "n_groups": ki["n_groups"],
            "n_study_results": ki["n_study_results"],
            "rod_program_groups": sorted({f"{p['dataset_id']}|{p['group_label']}" for p in rod}),
            "rod_program_study_results": int(sum(p["d18_study_results_to_those_markers"] for p in rod)),
        }
    ]


def post_build(release: str, *, neo4j: bool = False) -> dict:
    """Post-build step on an unpacked, not yet bundled release; the graph (KGX), summaries and registry are not touched.
    Rewrites validation/report.md from report.json + known_issues.json (small-group tables, Known issues section),
    optionally regenerates neo4j-import/ from the release KGX, and refreshes the manifest (counts, file checksums, code)."""
    from . import known_issues as KI
    from .biolink import toolkit
    from .neo4j_csv import write_import
    from .release import file_list, release_counts

    root = release_dir(release)
    if (data_root() / "kg" / "bundles" / f"{release}.tar").exists():
        raise SystemExit(f"{release} is already bundled; remove the bundle first (post-build regenerates it)")
    m = json.loads((root / "manifest.json").read_text())
    rep = json.loads((root / "validation" / "report.json").read_text())
    ki = json.loads((root / "validation" / "known_issues.json").read_text())
    md = (root / "validation" / "report.md").read_text().split("\n## Known issues\n")[0].rstrip("\n") + "\n"
    sg = small_groups_md(rep["ingest"]["heca"].get("edge_rule", {}).get("small_groups"))
    if sg and sg[1] not in md:
        head, tail = md.split("\n## Ingests\n", 1)
        md = head.rstrip("\n") + "\n" + "\n".join(sg) + "\n\n## Ingests\n" + tail
    (root / "validation" / "report.md").write_text(md + KI.markdown(ki))
    neo_counts = m.get("counts", {}).get("neo4j")
    if neo4j:
        tmp = root / "neo4j-import.new"
        if tmp.exists():
            shutil.rmtree(tmp)
        neo = write_import(root / "kgx" / "nodes.jsonl", root / "kgx" / "edges.jsonl", tmp, toolkit())
        if (root / "neo4j-import").exists():
            shutil.rmtree(root / "neo4j-import")
        tmp.rename(root / "neo4j-import")
        neo_counts = neo["counts"]
    code = m["code"]["visionkb-kg"]
    sha = _git(str(Path(__file__).resolve().parents[2]))
    if sha and sha != code.get("git_sha"):
        code.setdefault("graph_build_git_sha", code.get("git_sha"))
        code["git_sha"] = sha
    steps = m.setdefault("post_build", [])
    steps.append(
        {
            "git_sha": sha,
            "utc": datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "regenerated": ["validation/report.md"] + (["neo4j-import/"] if neo4j else []),
            "unchanged": ["kgx/", "summaries/", "registry/", "validation/report.json", "validation/known_issues.json"],
        }
    )
    m["counts"] = {**release_counts(root), "neo4j": neo_counts}
    m["files"] = file_list(root)
    wjson(m, root / "manifest.json")
    return {
        "code": m["code"],
        "post_build": steps[-1],
        "counts": {k: v for k, v in m["counts"].items() if k != "neo4j"},
    }


def refresh_manifest(release: str) -> dict:
    """Recompute manifest counts (contracts v0.7) and the file list of an unpacked, not yet bundled release."""
    from .release import file_list, release_counts

    root = release_dir(release)
    if (data_root() / "kg" / "bundles" / f"{release}.tar").exists():
        raise SystemExit(f"{release} is already bundled; releases are immutable after pack")
    m = json.loads((root / "manifest.json").read_text())
    m["counts"] = {**release_counts(root), "neo4j": m.get("counts", {}).get("neo4j")}
    m["files"] = file_list(root)
    wjson(m, root / "manifest.json")
    return m["counts"]


def small_groups_md(sg: dict | None) -> list[str]:
    """Markdown lines of the D18 small-group report (groups below the cell threshold)."""
    if not sg:
        return []
    L: list[str] = []
    L += [
        "",
        f"### Small groups (< {sg['threshold_n_cells']} cells)",
        "",
        f"{len(sg['groups'])} groups below {sg['threshold_n_cells']} cells. Rows removed only by (c) over all groups: "
        f"{sg['only_c_rows_total']:,}; {sg['only_c_rows_in_groups_with_median_cells_per_donor_lt_10']:,} of them in "
        f"{sg['groups_affected_by_only_c_with_median_lt_10']} groups with median < 10 cells per donor "
        f"({sg['groups_with_median_cells_per_donor_lt_10']} groups have median < 10 cells per donor).",
        "",
        "| dataset | level | group | cells | donors | cells/donor median (min) | candidates | fail a | fail b | fail c "
        "| kept | only a | only c | top kept | top removed only by (c) | top removed only by (a) |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for g in sg["groups"]:
        L.append(
            f"| {g['dataset_id'].removeprefix('visionkb:ds-')} | {g['level']} | {g['group_label']} | {g['n_cells']} | "
            f"{g['n_donors']} | {g['median_cells_per_donor']:g} ({g['min_cells_per_donor']}) | {g['candidate_rows']:,} | "
            f"{g['fail_a']:,} | {g['fail_b']:,} | {g['fail_c']:,} | {g['kept']:,} | {g['only_a']:,} | {g['only_c']:,} | "
            f"{', '.join(g['top_kept']) or '-'} | {', '.join(g['top_removed_only_c']) or '-'} | "
            f"{', '.join(g['top_removed_only_a']) or '-'} |"
        )
    return L


def report_md(rep: dict) -> str:
    L = [
        f"# Validation report: {rep['release']}",
        "",
        f"Status: **{rep['status']}**; errors {rep['n_errors']}, warnings {rep['n_warnings']}.",
        f"Biolink {rep['biolink_version']} (bmt, vendored schema).",
        "",
        "## Nodes by category",
        "",
        "| category | n |",
        "|---|---|",
    ]
    L += [f"| {k} | {v:,} |" for k, v in rep["nodes_by_category"].items()]
    L += ["", "## Edges by predicate", "", "| predicate | n |", "|---|---|"]
    L += [f"| {k} | {v:,} |" for k, v in rep["edges_by_predicate"].items()]
    L += ["", "## Category-predicate triples (bmt domain/range)", "", "| triple | n | bmt valid |", "|---|---|---|"]
    L += [f"| {k} | {v['count']:,} | {v['bmt_valid']} |" for k, v in rep["category_predicate_triples"].items()]
    for lvl in ("errors", "warnings"):
        L += ["", f"## {lvl.capitalize()}", ""]
        if not rep[lvl]:
            L.append("none")
        for k, v in rep[lvl].items():
            L.append(f"- `{k}`: {v['count']:,} (e.g. {json.dumps(v['examples'][:3], default=str)[:400]})")
    st = rep["self_tests"]["apln_wrong_id"]
    L += [
        "",
        "## Self-test: wrong APLN id",
        "",
        f"Probe gene node `HGNC:16693` with symbol APLN (the design proposal's id): flagged = {st['flagged']} "
        f"(HGNC:16693 is {st['hgnc_symbol_of_HGNC:16693']}); correct id HGNC:16665 passes = {st['correct_id_passes']}. Status {st['status']}.",
    ]
    er = rep["ingest"].get("heca", {}).get("edge_rule")
    if er:
        r = er["rule"]
        L += [
            "",
            f"## Expression edge rule {r['decision']}",
            "",
            f"(a) pct_expressing >= {r['min_pct_expressing']} and n_cells >= {r['min_n_cells']}; (b) pseudobulk_cp10k >= "
            f"{r['specificity_fraction']} x max over groups of the same dataset and level; (c) n_donors_expressing >= "
            f"{r['min_n_donors_expressing']}. Levels: {', '.join(r['levels'])} (region_majorclass rows stay in the Parquet only).",
            "",
            "Donor support checked: "
            + ", ".join(f"{k.removeprefix('visionkb:ds-')}={v}" for k, v in er["donor_support_checked"].items())
            + (
                ". Atlases without a usable donor column: "
                + ", ".join(k for k, v in er["donor_support_checked"].items() if not v)
                if not all(er["donor_support_checked"].values())
                else ". All atlases have a donor column; (c) applied everywhere."
            ),
            "",
            "Rows removed per rule (a, b, c counted separately; 'any' = removed by at least one); candidates = all gene x group rows of the level with an HGNC gene.",
            "",
            "| dataset | level | candidates | non-HGNC rows | fail a | fail b | fail c | fail any | kept D18 | kept D8 (old) |",
            "|---|---|---|---|---|---|---|---|---|---|",
        ]
        for x in er["by_dataset_level"]:
            L.append(
                f"| {x['dataset_id'].removeprefix('visionkb:ds-')} | {x['level']} | {x['candidates_hgnc']:,} | {x['non_hgnc_gene_rows']:,} | {x['fail_a_detection']:,} | "
                f"{x['fail_b_specificity']:,} | {x['fail_c_donor_support']:,} | {x['fail_any']:,} | {x['kept']:,} | {x['kept_old_rule_D8']:,} |"
            )
        t = er["totals"]
        L.append(
            f"| **total** | | {t['candidates_hgnc']:,} | {t['non_hgnc_gene_rows']:,} | {t['fail_a_detection']:,} | {t['fail_b_specificity']:,} | "
            f"{t['fail_c_donor_support']:,} | {t['fail_any']:,} | {t['kept']:,} | {t['kept_old_rule_D8']:,} |"
        )
        L += [
            "",
            "Robustness of the specificity fraction (report only; the release uses 0.10):",
            "",
            "| fraction | study results | gene x cell type x dataset |",
            "|---|---|---|",
        ]
        L += [
            f"| {k} | {v['study_results']:,} | {v['edges_gene_celltype_dataset_tissue']:,} |"
            for k, v in er["robustness"].items()
        ]
        L += ["", "Sanity checks: majorclass groups with an edge, old rule D8 -> D18 (celltype-level counts after)", ""]
        for g, v in er["sanity"].items():
            L.append(
                f"- **{g}**: "
                + "; ".join(
                    f"{ds}: {', '.join(d.get('D8', [])) or '-'} -> {', '.join(d.get('D18', [])) or '-'}"
                    for ds, d in sorted(v["majorclass"].items())
                )
                + f" | celltype groups D8 {v['celltype_groups']['D8']} -> D18 {v['celltype_groups']['D18']}"
            )
        L += small_groups_md(er.get("small_groups"))
    L += ["", "## Ingests", ""]
    for s, v in rep["ingest"].items():
        brief = {
            k: v[k]
            for k in v
            if not isinstance(v[k], list | dict) or k in ("versions", "mapping", "edges_by_predicate", "rows_by_level")
        }
        L.append(f"- **{s}**: {json.dumps(brief, default=str)[:1200]}")
        for k in (
            "flagged",
            "extra_ids_not_approved",
            "referenced_not_found_or_obsolete",
            "cell_states_without_parent",
        ):
            if v.get(k):
                L.append(f"  - {k}: {json.dumps(v[k], default=str)[:1500]}")
    return "\n".join(L) + "\n"


def main(argv=None):
    p = argparse.ArgumentParser(prog="visionkb-kg", description=__doc__)
    p.add_argument("--version", action="version", version=__version__)
    p.add_argument("--sources", help="sources.yaml (default $VISIONKB_DATA_ROOT/config/sources.yaml)")
    p.add_argument(
        "--work-dir", help="build working root (default $VISIONKB_KG_WORK or $VISIONKB_WORK_DIR/visionkb-kg)"
    )
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("ingest")
    s.add_argument("source", choices=INGESTS)
    s.add_argument("--release", required=True)
    s = sub.add_parser("build")
    s.add_argument("release")
    s.add_argument("--no-kgx", action="store_true", help="skip the KGX validator pass")
    s = sub.add_parser("refresh-manifest", help="recompute manifest counts and file checksums before pack")
    s.add_argument("release")
    s = sub.add_parser(
        "post-build", help="rewrite report.md, optionally regenerate neo4j-import/, refresh manifest (before pack)"
    )
    s.add_argument("release")
    s.add_argument("--neo4j", action="store_true", help="regenerate neo4j-import/ from the release KGX")
    s = sub.add_parser("pack")
    s.add_argument("release")
    s.add_argument("--out", help="bundle dir (default $VISIONKB_DATA_ROOT/kg/bundles)")
    s = sub.add_parser("verify")
    s.add_argument("bundle")
    s.add_argument("--extract-to")
    s = sub.add_parser("all")
    s.add_argument("release")
    args = p.parse_args(argv)
    import re

    rel = getattr(args, "release", None)
    if rel and not re.fullmatch(r"visionkb-\d{4}\.\d{2}\.\d+", rel):
        raise SystemExit(f"release id must look like visionkb-YYYY.MM.N: {rel}")
    if args.work_dir:
        import os

        os.environ["VISIONKB_KG_WORK"] = str(Path(args.work_dir).resolve())
    if args.cmd == "verify":
        from .release import verify

        res = verify(Path(args.bundle), Path(args.extract_to) if args.extract_to else None)
        print(json.dumps(res, indent=2))
        sys.exit(0 if res["status"] == "pass" else 1)
    if args.cmd == "refresh-manifest":
        print(json.dumps({k: v for k, v in refresh_manifest(args.release).items() if k != "neo4j"}, indent=2))
        return
    if args.cmd == "post-build":
        print(json.dumps(post_build(args.release, neo4j=args.neo4j), indent=2))
        return
    if args.cmd == "pack":
        from .release import pack

        out = Path(args.out) if args.out else data_root() / "kg" / "bundles"
        print(json.dumps(pack(release_dir(args.release), out), indent=2))
        return
    src = Sources.load(args.sources)
    if args.cmd == "ingest":
        ingest(args.source, args.release, src)
    elif args.cmd == "build":
        build(args.release, src, run_kgx=not args.no_kgx)
    elif args.cmd == "all":
        for s in ("heca", "retigene", "datasets", "hgnc", "ontologies"):
            ingest(s, args.release, src)
        build(args.release, src)


if __name__ == "__main__":
    main()
