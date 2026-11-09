"""HGNC genes: all approved protein-coding genes plus any HGNC id another source uses."""

from __future__ import annotations

import pandas as pd

from .common import HUMAN

XREF_COLS = {"ensembl_gene_id": "ENSEMBL", "entrez_id": "NCBIGene", "omim_id": "OMIM", "uniprot_ids": "UniProtKB"}


def load(path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False).fillna("")
    return df


def split(v: str) -> list[str]:
    return [x.strip().strip('"') for x in v.split("|") if x.strip()] if v else []


def symbol_index(df: pd.DataFrame) -> dict[str, str]:
    a = df[df.status == "Approved"]
    return dict(zip(a.symbol, a.hgnc_id, strict=True))


def nodes(df: pd.DataFrame, extra_ids: set[str]) -> tuple[list[dict], dict]:
    a = df[df.status == "Approved"]
    keep = a[(a.locus_group == "protein-coding gene") | a.hgnc_id.isin(extra_ids)]
    missing = sorted(extra_ids - set(a.hgnc_id))
    out = []
    for r in keep.itertuples(index=False):
        xref = []
        for col, pre in XREF_COLS.items():
            xref += [f"{pre}:{x}" for x in split(getattr(r, col))]
        out.append(
            {
                "id": r.hgnc_id,
                "category": ["biolink:Gene"],
                "name": r.symbol,
                "symbol": r.symbol,
                "full_name": r.name,
                "synonym": sorted(set(split(r.alias_symbol) + split(r.prev_symbol))),
                "xref": sorted(set(xref)),
                "in_taxon": HUMAN,
                "locus_group": r.locus_group,
                "locus_type": r.locus_type,
                "provided_by": ["infores:hgnc"],
            }
        )
    stats = {
        "approved": len(a),
        "protein_coding": int((a.locus_group == "protein-coding gene").sum()),
        "nodes": len(out),
        "extra_ids_requested": len(extra_ids),
        "extra_ids_not_approved": missing,
    }
    return out, stats
