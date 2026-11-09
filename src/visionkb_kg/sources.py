"""Build inputs (sources.yaml, contracts §1) and environment roots (contracts §9)."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import yaml

PACKAGE_CROSSWALK = Path(__file__).parent / "mappings" / "heca-celltypes.sssom.tsv"
INFORES = {
    "hgnc": "infores:hgnc",
    "cl": "infores:cl",
    "uberon": "infores:uberon",
    "mondo": "infores:mondo",
    "hp": "infores:hpo",
    "efo": "infores:efo",
    "ncbitaxon": "infores:ncbi-taxon",
    "heca": "infores:visionkb-heca",
    "retigene": "infores:visionkb-retigene",
}


def dirname_of(curie: str) -> str:
    """File/directory name for a CURIE: the local part without prefix (no ':' in paths; contracts v0.4)."""
    if ":" not in curie:
        raise ValueError(f"not a CURIE: {curie}")
    name = curie.split(":", 1)[1]
    if not name or "/" in name or ":" in name or name in (".", ".."):
        raise ValueError(f"CURIE {curie!r} has no safe file name")
    return name


def data_root() -> Path:
    v = os.environ.get("VISIONKB_DATA_ROOT")
    if not v:
        raise RuntimeError("VISIONKB_DATA_ROOT is not set")
    return Path(v)


def work_root() -> Path:
    """Build working files: $VISIONKB_KG_WORK if set (e.g. a workflow step folder), else $VISIONKB_WORK_DIR/visionkb-kg."""
    if v := os.environ.get("VISIONKB_KG_WORK"):
        return Path(v)
    v = os.environ.get("VISIONKB_WORK_DIR")
    if not v:
        raise RuntimeError("VISIONKB_WORK_DIR (or VISIONKB_KG_WORK) is not set")
    return Path(v) / "visionkb-kg"


def build_dir(release: str) -> Path:
    return work_root() / "build" / release


def release_dir(release: str) -> Path:
    return data_root() / "kg" / "releases" / release


def _expand(v):
    if isinstance(v, str):
        return os.path.expandvars(v)
    if isinstance(v, dict):
        return {k: _expand(x) for k, x in v.items()}
    if isinstance(v, list):
        return [_expand(x) for x in v]
    return v


@dataclass
class Sources:
    raw: dict

    @classmethod
    def load(cls, path: str | Path | None = None) -> Sources:
        p = Path(path) if path else data_root() / "config" / "sources.yaml"
        return cls(_expand(yaml.safe_load(p.read_text())))

    def get(self, *keys, default=None):
        x = self.raw
        for k in keys:
            if not isinstance(x, dict) or k not in x:
                return default
            x = x[k]
        return x

    @property
    def crosswalk(self) -> Path:
        return Path(self.get("crosswalk", default=str(PACKAGE_CROSSWALK)))
