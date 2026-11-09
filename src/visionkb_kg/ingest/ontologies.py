"""Ontology slims: CL, UBERON, MONDO, HP from OBO releases; EFO and NCBITaxon terms from OLS.

Slim = terms referenced by any ingest + their is_a/part_of ancestors + the eye branch of UBERON (is_a/part_of
descendants of UBERON:0000019 'camera-type eye') + CL cell types of the retina/eye (descendants of CL:0009004
'retinal cell' and CL classes part_of an eye-branch structure). Edges: subclass_of (is_a) and part_of between kept terms.
"""

from __future__ import annotations

import time
import urllib.parse

import requests

from ..obo import ancestors, descendants, read_obo, version_of
from ..sources import INFORES
from .common import make_edge

EYE_ROOT = "UBERON:0000019"
RETINAL_CELL = "CL:0009004"
CATEGORY = {
    "CL": "biolink:Cell",
    "UBERON": "biolink:AnatomicalEntity",
    "MONDO": "biolink:Disease",
    "HP": "biolink:PhenotypicFeature",
    "EFO": "biolink:NamedThing",
    "NCBITaxon": "biolink:OrganismTaxon",
}
OBO_SOURCES = {"CL": "cl", "UBERON": "uberon", "MONDO": "mondo", "HP": "hp"}
OLS_SOURCES = {"EFO": "efo", "NCBITaxon": "ncbitaxon"}
OLS = "https://www.ebi.ac.uk/ols4/api"


def load_obo_sets(paths: dict[str, str]) -> dict[str, dict]:
    """{'CL': {'header':..., 'terms':...}, ...} for the four OBO ontologies (terms filtered to the ontology's prefix)."""
    out = {}
    for pre, key in OBO_SOURCES.items():
        header, terms = read_obo(paths[key], prefixes=(pre,))
        out[pre] = {"header": header, "terms": terms, "version": version_of(header), "path": paths[key]}
    # CL also carries cross-ontology part_of to UBERON; keep the raw CL file's relations
    return out


def _iri(curie: str) -> str:
    pre, loc = curie.split(":", 1)
    if pre == "EFO":
        return f"http://www.ebi.ac.uk/efo/EFO_{loc}"
    return f"http://purl.obolibrary.org/obo/{pre}_{loc}"


def _get(url, params=None):
    for attempt in range(4):
        try:
            r = requests.get(url, params=params, timeout=60)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.json()
        except requests.RequestException:
            if attempt == 3:
                raise
            time.sleep(2 * (attempt + 1))
    return None


def ols_term(curie: str) -> dict | None:
    ont = OLS_SOURCES[curie.split(":")[0]]
    iri2 = urllib.parse.quote(urllib.parse.quote(_iri(curie), safe=""), safe="")
    t = _get(f"{OLS}/ontologies/{ont}/terms/{iri2}")
    if not t:
        return None
    parents = _get(f"{OLS}/ontologies/{ont}/terms/{iri2}/parents", {"size": 100}) or {}
    pids = [p.get("obo_id") for p in parents.get("_embedded", {}).get("terms", []) if p.get("obo_id")]
    return {
        "id": curie,
        "name": t.get("label"),
        "obsolete": bool(t.get("is_obsolete")),
        "synonyms": t.get("synonyms") or [],
        "def": (t.get("description") or [None])[0],
        "is_a": pids,
    }


def ols_version(ont: str) -> str:
    d = _get(f"{OLS}/ontologies/{ont}") or {}
    cfg = d.get("config", {})
    return cfg.get("version") or cfg.get("versionIri") or d.get("updated") or "unknown"


def ols_slim(seeds: set[str], max_terms: int = 400) -> tuple[dict, dict]:
    """Seed terms + is_a ancestors (same prefix) via OLS."""
    terms, todo = {}, sorted(seeds)
    while todo and len(terms) < max_terms:
        c = todo.pop()
        if c in terms:
            continue
        t = ols_term(c)
        if t is None:
            terms[c] = None
            continue
        t["is_a"] = [p for p in t["is_a"] if p.split(":")[0] == c.split(":")[0]]
        terms[c] = t
        todo.extend(p for p in t["is_a"] if p not in terms)
    versions = {
        pre: ols_version(ont) for pre, ont in OLS_SOURCES.items() if any(s.startswith(pre + ":") for s in seeds)
    }
    return terms, versions


def slim(referenced: set[str], obo: dict, release: str) -> tuple[list[dict], list[dict], dict]:
    cl, ub = obo["CL"]["terms"], obo["UBERON"]["terms"]
    eye = descendants(ub, EYE_ROOT)
    eye_cells = descendants(cl, RETINAL_CELL, rels=("is_a",))
    # CL classes with part_of to an eye-branch structure (cross-ontology relations in cl.obo)
    _, cl_all = read_obo(obo["CL"]["path"])
    for tid, rec in cl_all.items():
        if tid.startswith("CL:") and not rec["obsolete"] and any(p in eye for p in rec["part_of"]):
            eye_cells.add(tid)
            cl.setdefault(tid, rec)
    for tid, rec in cl_all.items():
        if tid in cl:
            cl[tid]["part_of"] = sorted(set(cl[tid]["part_of"]) | set(rec["part_of"]))
    seeds = {c for c in referenced if c.split(":")[0] in OBO_SOURCES} | eye | eye_cells
    keep: dict[str, dict] = {}
    unknown = []
    all_terms = {pre: obo[pre]["terms"] for pre in OBO_SOURCES}
    closure = set()
    for s in seeds:
        pre = s.split(":")[0]
        terms = all_terms[pre]
        if s not in terms or terms[s]["obsolete"]:
            unknown.append(s)
            continue
        closure |= ancestors(terms, [s]) if pre != "CL" else ancestors({**cl, **ub}, [s])
    for c in sorted(closure):
        pre = c.split(":")[0]
        rec = all_terms.get(pre, {}).get(c)
        if rec and not rec["obsolete"]:
            keep[c] = rec
    nodes, edges = [], []
    for c, rec in sorted(keep.items()):
        pre = c.split(":")[0]
        nodes.append(
            {
                "id": c,
                "category": [CATEGORY[pre]],
                "name": rec.get("name"),
                "description": rec.get("def"),
                "synonym": sorted(set(rec["synonyms"])),
                "xref": sorted(set(rec["xrefs"])),
                "provided_by": [INFORES[OBO_SOURCES[pre]]],
                "in_eye_branch": c in eye or c in eye_cells or None,
            }
        )
        src = INFORES[OBO_SOURCES[pre]]
        for p in rec["is_a"]:
            if p in keep:
                edges.append(
                    make_edge(
                        c,
                        "biolink:subclass_of",
                        p,
                        source=src,
                        knowledge_level="knowledge_assertion",
                        agent_type="manual_agent",
                        release=release,
                    )
                )
        for p in rec["part_of"]:
            if p in keep:
                edges.append(
                    make_edge(
                        c,
                        "biolink:part_of",
                        p,
                        source=src,
                        knowledge_level="knowledge_assertion",
                        agent_type="manual_agent",
                        release=release,
                    )
                )
    # OLS-backed ontologies
    ols_seeds = {c for c in referenced if c.split(":")[0] in OLS_SOURCES}
    ols_terms, ols_versions = ols_slim(ols_seeds) if ols_seeds else ({}, {})
    for c, t in sorted(ols_terms.items()):
        if t is None or t["obsolete"]:
            unknown.append(c)
            continue
        pre = c.split(":")[0]
        nodes.append(
            {
                "id": c,
                "category": [CATEGORY[pre]],
                "name": t["name"],
                "description": t["def"],
                "synonym": sorted(set(t["synonyms"])),
                "provided_by": [INFORES[OLS_SOURCES[pre]]],
                "ols_verified": True,
            }
        )
        for p in t["is_a"]:
            if ols_terms.get(p):
                edges.append(
                    make_edge(
                        c,
                        "biolink:subclass_of",
                        p,
                        source=INFORES[OLS_SOURCES[pre]],
                        knowledge_level="knowledge_assertion",
                        agent_type="manual_agent",
                        release=release,
                    )
                )
    ids = [e["id"] for e in edges]
    edges = list({e["id"]: e for e in edges}.values())
    stats = {
        "versions": {**{pre: obo[pre]["version"] for pre in OBO_SOURCES}, **ols_versions},
        "referenced": len(referenced),
        "eye_branch_uberon": len(eye),
        "eye_cell_types": len(eye_cells),
        "nodes_by_prefix": {pre: sum(n["id"].startswith(pre + ":") for n in nodes) for pre in CATEGORY},
        "edges": len(edges),
        "duplicate_edges_dropped": len(ids) - len(edges),
        "referenced_not_found_or_obsolete": sorted(set(unknown)),
    }
    return nodes, edges, stats
