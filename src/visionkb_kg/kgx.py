"""KGX JSONL read/write, deterministic (sorted records and keys)."""

from __future__ import annotations

import json
from collections.abc import Iterable, Iterator
from pathlib import Path


def _clean(rec: dict) -> dict:
    return {k: v for k, v in rec.items() if v is not None and v != [] and v != ""}


def dumps(rec: dict) -> str:
    return json.dumps(_clean(rec), sort_keys=True, ensure_ascii=False, separators=(",", ":"))


def write_jsonl(records: Iterable[dict], path: Path, key: str = "id") -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    recs = sorted((_clean(r) for r in records), key=lambda r: r[key])
    ids = [r[key] for r in recs]
    if len(set(ids)) != len(ids):
        raise ValueError(f"{path}: duplicate {key}")
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as f:
        for r in recs:
            f.write(dumps(r) + "\n")
    tmp.replace(path)
    return len(recs)


def read_jsonl(path: Path) -> Iterator[dict]:
    with open(path) as f:
        for line in f:
            if line.strip():
                yield json.loads(line)


def write_source(outdir: Path, nodes: Iterable[dict], edges: Iterable[dict]) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    n = write_jsonl(nodes, outdir / "nodes.jsonl")
    e = write_jsonl(edges, outdir / "edges.jsonl")
    return {"nodes": n, "edges": e}
