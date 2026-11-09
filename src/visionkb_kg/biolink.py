"""Biolink Model 4.4.6 pinned: the vendored schema and a bmt Toolkit built on it (bmt's default bundles another version)."""

from __future__ import annotations

import hashlib
from functools import lru_cache
from pathlib import Path

SCHEMA_DIR = Path(__file__).parent / "schema"
BIOLINK_YAML = SCHEMA_DIR / "biolink" / "biolink-model_v4.4.6.yaml"
BIOLINK_SHA256 = "5837907e1eb579812856579feea5a099665252f2454bd6f950d63fbb00f4c9f0"
BIOLINK_URL = "https://raw.githubusercontent.com/biolink/biolink-model/v4.4.6/biolink-model.yaml"
VISIONKB_YAML = SCHEMA_DIR / "visionkb.yaml"


def biolink_sha256() -> str:
    return hashlib.sha256(BIOLINK_YAML.read_bytes()).hexdigest()


@lru_cache(maxsize=1)
def toolkit():
    from bmt import Toolkit

    t = Toolkit(schema=str(BIOLINK_YAML))
    v = t.get_model_version()
    if v != "4.4.6":
        raise RuntimeError(f"bmt loaded Biolink {v}, expected 4.4.6")
    return t


def visionkb_classes() -> dict[str, str]:
    """visionkb class name -> Biolink parent class name (from the LinkML schema, no imports)."""
    import yaml

    s = yaml.safe_load(VISIONKB_YAML.read_text())
    return {k: v.get("is_a") for k, v in s["classes"].items()}
