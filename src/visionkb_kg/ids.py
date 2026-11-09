"""Identifiers: CURIE checks, deterministic edge and study-result ids (contracts §5)."""

from __future__ import annotations

import hashlib
import json
import re

CURIE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*:[^\s]+$")
PATTERNS = {
    "HGNC": re.compile(r"^HGNC:\d+$"),
    "CL": re.compile(r"^CL:\d{7}$"),
    "UBERON": re.compile(r"^UBERON:\d{7}$"),
    "MONDO": re.compile(r"^MONDO:\d{7}$"),
    "HP": re.compile(r"^HP:\d{7}$"),
    "EFO": re.compile(r"^EFO:\d{7}$"),
    "NCBITaxon": re.compile(r"^NCBITaxon:\d+$"),
    "PMID": re.compile(r"^PMID:\d+$"),
    "infores": re.compile(r"^infores:[a-z0-9][a-z0-9.-]*$"),
    "visionkb": re.compile(r"^visionkb:[a-z]+-[A-Za-z0-9._-]+$"),
    "ENSEMBL": re.compile(r"^ENSEMBL:ENSG\d{11}$"),
}
QUALIFIER_KEYS = (
    "anatomical_context_qualifier",
    "species_context_qualifier",
    "stage_qualifier",
    "object_aspect_qualifier",
    "object_direction_qualifier",
    "qualified_predicate",
    "subject_form_or_variant_qualifier",
    "subject_aspect_qualifier",
)


def prefix(curie: str) -> str:
    return curie.split(":", 1)[0]


def check_curie(curie: str) -> str | None:
    """None if the CURIE is well-formed for its prefix, else a reason."""
    if not isinstance(curie, str) or not CURIE.match(curie):
        return "not a CURIE"
    p = PATTERNS.get(prefix(curie))
    if p and not p.match(curie):
        return f"malformed {prefix(curie)} id"
    return None


def _canon(v):
    if isinstance(v, list | tuple | set):
        return sorted(str(x) for x in v)
    return str(v)


def qualifiers_of(edge: dict) -> dict:
    return {k: _canon(edge[k]) for k in QUALIFIER_KEYS if edge.get(k) not in (None, [], "")}


def edge_id(subject: str, predicate: str, obj: str, qualifiers: dict | None = None) -> str:
    """sha1 of subject|predicate|object|sorted qualifiers, prefixed visionkb:e- (contracts §5)."""
    q = json.dumps(
        {k: _canon(v) for k, v in sorted((qualifiers or {}).items()) if v not in (None, [], "")},
        sort_keys=True,
        separators=(",", ":"),
    )
    return "visionkb:e-" + hashlib.sha1(f"{subject}|{predicate}|{obj}|{q}".encode()).hexdigest()  # noqa: S324 - id, not security


def short_id(kind: str, *parts) -> str:
    h = hashlib.sha1("|".join(str(p) for p in parts).encode()).hexdigest()[:20]  # noqa: S324 - id, not security
    return f"visionkb:{kind}-{h}"
