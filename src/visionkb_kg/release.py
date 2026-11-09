"""Merge per-source KGX, assemble the release folder (contracts §1, §6), pack and verify the portable bundle (D15-D17)."""

from __future__ import annotations

import hashlib
import json
import os
import tarfile
from collections import Counter
from pathlib import Path

from .kgx import dumps, read_jsonl

LIST_KEYS = {"provided_by", "category", "synonym", "xref"}


def sha256_file(path: Path, bufsize: int = 16 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(bufsize):
            h.update(chunk)
    return h.hexdigest()


def merge(source_dirs: list[Path], out_dir: Path) -> dict:
    """Nodes merged by id (first source wins for scalars; list fields united); edges deduplicated by id."""
    nodes: dict[str, dict] = {}
    edges: dict[str, dict] = {}
    stats = {"sources": {}, "node_merges": 0, "edge_duplicates": 0}
    for d in source_dirs:
        n_n = n_e = 0
        for n in read_jsonl(d / "nodes.jsonl"):
            n_n += 1
            cur = nodes.get(n["id"])
            if cur is None:
                nodes[n["id"]] = n
                continue
            stats["node_merges"] += 1
            for k, v in n.items():
                if k in LIST_KEYS and isinstance(v, list):
                    merged = list(cur.get(k) or [])
                    merged += [x for x in v if x not in merged]
                    cur[k] = merged if k == "category" else sorted(merged)
                elif cur.get(k) in (None, "", []):
                    cur[k] = v
        for e in read_jsonl(d / "edges.jsonl"):
            n_e += 1
            if e["id"] in edges:
                stats["edge_duplicates"] += 1
                continue
            edges[e["id"]] = e
        stats["sources"][d.name] = {"nodes": n_n, "edges": n_e}
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, recs in (("nodes.jsonl", nodes), ("edges.jsonl", edges)):
        tmp = out_dir / (name + ".tmp")
        with open(tmp, "w") as f:
            for k in sorted(recs):
                f.write(dumps(recs[k]) + "\n")
        tmp.replace(out_dir / name)
    stats["nodes"] = len(nodes)
    stats["edges"] = len(edges)
    stats["nodes_by_category"] = dict(
        sorted(
            Counter(
                n["category"][-1] if n["category"][0] == "biolink:Cell" and len(n["category"]) > 1 else n["category"][0]
                for n in nodes.values()
            ).items()
        )
    )
    stats["edges_by_predicate"] = dict(sorted(Counter(e["predicate"] for e in edges.values()).items()))
    return stats


def most_specific(categories: list[str]) -> str:
    """One category per node: a visionkb class if present (it extends a Biolink class), else the first Biolink category."""
    vk = [c for c in categories if c.startswith("visionkb:")]
    return vk[0] if vk else categories[0]


def release_counts(root: Path) -> dict:
    """Manifest counts (contracts v0.7): every node once under its most specific category, every relationship written to
    neo4j-import once under its predicate (StudyResult -> Dataset links as ``supporting_data_set``); totals are checked
    against the row counts of the neo4j-import CSVs."""
    import csv

    by_cat: Counter = Counter()
    for n in read_jsonl(root / "kgx" / "nodes.jsonl"):
        by_cat[most_specific(n["category"])] += 1
    by_pred: Counter = Counter(e["predicate"] for e in read_jsonl(root / "kgx" / "edges.jsonl"))
    imp = root / "neo4j-import"

    def rows(p: Path) -> int:
        with open(p, newline="", encoding="utf-8") as f:
            return sum(1 for _ in csv.reader(f)) - 1

    sds = imp / "rels_supporting_data_set.csv"
    if sds.exists():
        by_pred["supporting_data_set"] = rows(sds)
    csv_nodes = sum(rows(p) for p in sorted(imp.glob("nodes_*.csv")))
    csv_rels = sum(rows(p) for p in sorted(imp.glob("rels_*.csv")))
    out = {
        "nodes_total": sum(by_cat.values()),
        "edges_total": sum(by_pred.values()),
        "study_results": by_cat.get("biolink:StudyResult", 0),
        "nodes_by_category": dict(sorted(by_cat.items())),
        "edges_by_predicate": dict(sorted(by_pred.items())),
    }
    if out["nodes_total"] != csv_nodes or out["edges_total"] != csv_rels:
        raise RuntimeError(
            f"KGX counts {out['nodes_total']}/{out['edges_total']} != neo4j-import rows {csv_nodes}/{csv_rels}"
        )
    return out


def file_list(root: Path, exclude=("manifest.json",)) -> list[dict]:
    out = []
    for p in sorted(root.rglob("*")):
        if p.is_file() and p.relative_to(root).as_posix() not in exclude:
            out.append({"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size, "sha256": sha256_file(p)})
    return out


def _abs_strings(obj, where="") -> list[str]:
    """Strings in a JSON value that are absolute filesystem paths (URLs are fine)."""
    if isinstance(obj, str):
        return [f"{where}: {obj}"] if obj.startswith(("/", "~")) or obj.startswith("file:") else []
    if isinstance(obj, dict):
        return [x for k, v in obj.items() for x in _abs_strings(v, f"{where}.{k}")]
    if isinstance(obj, list):
        return [x for i, v in enumerate(obj) for x in _abs_strings(v, f"{where}[{i}]")]
    return []


def check_relative(root: Path) -> list[str]:
    """Absolute paths in the files the service reads (manifest.json, neo4j-import/import.args)."""
    bad = []
    m = root / "manifest.json"
    if m.exists():
        bad += [f"manifest.json{x}" for x in _abs_strings(json.loads(m.read_text()))]
    a = root / "neo4j-import" / "import.args"
    if a.exists():
        for line in a.read_text().splitlines():
            if not line.startswith(("--nodes", "--relationships")):
                continue
            spec = line.split("=", 1)[1]
            spec = spec.split("=", 1)[1] if "=" in spec else spec
            for v in spec.split(","):
                if os.path.isabs(v) or ".." in v.split("/") or v.startswith("~"):
                    bad.append(f"import.args: non-relative path {v}")
    return bad


def pack(release_root: Path, bundle_dir: Path) -> dict:
    release = release_root.name
    bundle_dir.mkdir(parents=True, exist_ok=True)
    bad = check_relative(release_root)
    if bad:
        raise RuntimeError(f"release has absolute paths: {bad}")
    tar_path = bundle_dir / f"{release}.tar"
    tmp = tar_path.with_suffix(".tar.tmp")

    def _norm(ti: tarfile.TarInfo) -> tarfile.TarInfo:
        ti.uid = ti.gid = 0
        ti.uname = ti.gname = ""
        ti.mtime = 0
        ti.mode = 0o755 if ti.isdir() else 0o644
        return ti

    with tarfile.open(tmp, "w", format=tarfile.PAX_FORMAT) as tf:
        tf.add(release_root, arcname=release, recursive=False, filter=_norm)
        for p in sorted(release_root.rglob("*")):
            tf.add(p, arcname=f"{release}/{p.relative_to(release_root).as_posix()}", recursive=False, filter=_norm)
    tmp.replace(tar_path)
    digest = sha256_file(tar_path)
    (bundle_dir / f"{release}.tar.sha256").write_text(f"{digest}  {release}.tar\n")
    return {"bundle": str(tar_path), "sha256": digest, "bytes": tar_path.stat().st_size}


def verify(bundle: Path, extract_to: Path | None = None) -> dict:
    """Check the outer sha256, member paths (relative, no '..'), every manifest checksum, and no absolute paths."""
    res = {"bundle": bundle.name, "checks": {}}
    side = bundle.with_name(bundle.name + ".sha256")
    want = side.read_text().split()[0] if side.exists() else None
    got = sha256_file(bundle)
    res["checks"]["outer_sha256"] = want == got
    with tarfile.open(bundle) as tf:
        members = tf.getmembers()
        names = [m.name for m in members]
        root = names[0].split("/")[0]
        res["release"] = root
        res["checks"]["relative_member_paths"] = all(
            not os.path.isabs(n) and ".." not in n.split("/") and n.split("/")[0] == root for n in names
        )
        res["checks"]["no_links"] = all(not (m.issym() or m.islnk()) for m in members)
        manifest = json.load(tf.extractfile(f"{root}/manifest.json"))
        listed = {f["path"]: f for f in manifest["files"]}
        bad = []
        for m in members:
            if not m.isfile():
                continue
            rel = m.name[len(root) + 1 :]
            if rel == "manifest.json":
                continue
            f = listed.get(rel)
            if f is None:
                bad.append(f"{rel}: not in manifest")
                continue
            h = hashlib.sha256()
            fh = tf.extractfile(m)
            while chunk := fh.read(16 << 20):
                h.update(chunk)
            if h.hexdigest() != f["sha256"] or m.size != f["bytes"]:
                bad.append(f"{rel}: checksum/size mismatch")
        in_tar = {m.name[len(root) + 1 :] for m in members if m.isfile()}
        bad += [f"{p}: listed but missing" for p in listed if p not in in_tar]
        res["checks"]["manifest_checksums"] = not bad
        res["problems"] = bad[:50]
        res["n_files"] = len(listed)
        if extract_to is not None:
            extract_to.mkdir(parents=True, exist_ok=True)
            tf.extractall(extract_to, filter="data")
            res["extracted_to"] = str(extract_to / root)
            res["checks"]["no_absolute_paths_in_service_files"] = not check_relative(extract_to / root)
    res["status"] = "pass" if all(res["checks"].values()) else "fail"
    return res
