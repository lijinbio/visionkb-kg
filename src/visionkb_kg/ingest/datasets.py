"""Dataset nodes from atlas registry records (contracts §4)."""

from __future__ import annotations

FLAT = (
    "version",
    "title",
    "species",
    "tissue",
    "assay",
    "disease",
    "stage",
    "n_cells",
    "n_cells_annotated",
    "n_donors",
    "n_samples",
    "n_genes",
    "levels",
    "cell_types",
    "format",
    "h5ad_path",
    "soma_uri",
    "checksum",
    "access_level",
    "license",
    "citation",
)


def node(rec: dict, release: str) -> dict:
    n = {
        "id": rec["id"],
        "category": ["biolink:Dataset"],
        "name": rec["title"],
        "provided_by": ["infores:visionkb-heca"],
        "kg_release": release,
    }
    for k in FLAT:
        if rec.get(k) not in (None, [], ""):
            n["disease_context" if k == "disease" else k] = rec[k]
    n["atlas_version"] = rec["version"]
    b = rec.get("bind") or {}
    n.update(bind_path=b.get("path"), bind_name=b.get("name"), bind_role=b.get("role"))
    return n


def registry_line(rec: dict, release: str) -> dict:
    out = dict(rec)
    out["kg_release"] = release
    return out
