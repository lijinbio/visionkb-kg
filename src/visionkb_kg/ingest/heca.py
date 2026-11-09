"""HECA atlas summaries -> expressed_in edges (gene -> CL class or CellState) with ExpressionStudyResult nodes.

One edge per (gene, cell type, anatomical context, species); each passing atlas group (D8: pct_expressing >= 0.10 and
n_cells >= 20) becomes one StudyResult linked to its Dataset. Level majorclass/celltype rows take the dataset tissues
as anatomical context; region_majorclass rows take the region's UBERON id (rows whose region has no id are skipped).
CellState nodes come from groups.parquet and are subclass_of their CL class.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as pds

from ..ids import edge_id, short_id
from .common import HUMAN, make_edge

SOURCE = "infores:visionkb-heca"
MIN_PCT = 0.10
MIN_CELLS = 20


def load_passing(summary_files: list[Path], min_pct=MIN_PCT, min_cells=MIN_CELLS) -> pd.DataFrame:
    parts = []
    for f in summary_files:
        d = pds.dataset(str(f))
        t = d.to_table(filter=(pds.field("pct_expressing") >= min_pct) & (pds.field("n_cells") >= min_cells))
        parts.append(t.to_pandas())
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def build(
    rows: pd.DataFrame, groups: pd.DataFrame, atlases: dict, release: str, method: str
) -> tuple[list, list, dict]:
    """rows: passing expr_summary rows; groups: concatenated groups.parquet; atlases: dataset_id -> {tissue, region_uberon, access_level}."""
    stats = {"rows_passing": len(rows)}
    rows = rows.copy()
    is_hgnc = rows["gene_id"].str.startswith("HGNC:")
    stats["rows_non_hgnc_gene_skipped"] = int((~is_hgnc).sum())
    rows = rows[is_hgnc]
    ctx = []
    for ds, level, region in zip(rows["dataset_id"], rows["level"], rows["region"], strict=True):
        a = atlases[ds]
        if level == "region_majorclass":
            u = a.get("region_uberon", {}).get(region)
            ctx.append((u,) if u else None)
        else:
            ctx.append(tuple(sorted(a["tissue"])))
    rows["ctx"] = ctx
    stats["rows_region_without_uberon_skipped"] = int(rows["ctx"].isna().sum())
    rows = rows[rows["ctx"].notna()].reset_index(drop=True)
    rows["ctx_key"] = rows["ctx"].map(lambda t: "|".join(t))
    # deterministic edge ids
    keys = (
        rows[["gene_id", "cell_type_id", "ctx_key"]]
        .drop_duplicates()
        .sort_values(["gene_id", "cell_type_id", "ctx_key"])
    )
    keys["edge_id"] = [
        edge_id(
            g,
            "biolink:expressed_in",
            c,
            {"anatomical_context_qualifier": k.split("|"), "species_context_qualifier": HUMAN},
        )
        for g, c, k in zip(keys.gene_id, keys.cell_type_id, keys.ctx_key, strict=True)
    ]
    rows = rows.merge(keys, on=["gene_id", "cell_type_id", "ctx_key"], how="left", validate="m:1")
    rows["sr_id"] = [
        short_id("sr", d, lv, g, gid)
        for d, lv, g, gid in zip(rows.dataset_id, rows.level, rows.group_label, rows.gene_id, strict=True)
    ]
    if rows["sr_id"].duplicated().any():
        raise RuntimeError("duplicate study result ids")
    srs = []
    for r in rows.itertuples(index=False):
        srs.append(
            {
                "id": r.sr_id,
                "category": ["biolink:StudyResult"],
                "name": f"{r.gene_symbol} in {r.group_label} ({r.dataset_id.removeprefix('visionkb:ds-')}, {r.level})",
                "edge_id": r.edge_id,
                "dataset_id": r.dataset_id,
                "atlas_version": r.atlas_version,
                "level": r.level,
                "group_label": r.group_label,
                "region": r.region if isinstance(r.region, str) else None,
                "n_cells": int(r.n_cells),
                "n_donors": int(r.n_donors),
                "pct_expressing": round(float(r.pct_expressing), 6),
                "mean_log1p_cp10k": round(float(r.mean_log1p_cp10k), 6),
                "pseudobulk_cp10k": round(float(r.pseudobulk_cp10k), 6),
                "n_donors_expressing": int(r.n_donors_expressing),
                "method": method,
                "supporting_data_set": [r.dataset_id],
                "provided_by": [SOURCE],
            }
        )
    vis = {ds: ("public" if a.get("access_level") == "public" else "lab") for ds, a in atlases.items()}
    # aggregate per edge with sorted columns (one pass; no per-group pandas objects)
    rows = rows.sort_values(["edge_id", "sr_id"], kind="mergesort").reset_index(drop=True)
    rows["public"] = rows["dataset_id"].map(lambda d: vis[d] == "public")
    agg = rows.groupby("edge_id", sort=True).agg(
        gene_id=("gene_id", "first"),
        cell_type_id=("cell_type_id", "first"),
        ctx=("ctx", "first"),
        sr_ids=("sr_id", list),
        datasets=("dataset_id", lambda x: sorted(set(x))),
        max_pct=("pct_expressing", "max"),
        max_mean=("mean_log1p_cp10k", "max"),
        all_public=("public", "all"),
    )
    edges = []
    for eid, r in zip(agg.index, agg.itertuples(index=False), strict=True):
        e = make_edge(
            r.gene_id,
            "biolink:expressed_in",
            r.cell_type_id,
            source=SOURCE,
            knowledge_level="observation",
            agent_type="data_analysis_pipeline",
            release=release,
            anatomical_context_qualifier=list(r.ctx),
            species_context_qualifier=HUMAN,
            has_study_results=r.sr_ids,
            evidence_count=len(r.sr_ids),
            n_datasets=len(r.datasets),
            max_pct_expressing=round(float(r.max_pct), 6),
            max_mean_log1p_cp10k=round(float(r.max_mean), 6),
            supporting_data_set=r.datasets,
            visibility="public" if r.all_public else "lab",
        )
        if e["id"] != eid:
            raise RuntimeError(f"edge id mismatch for {eid}")
        edges.append(e)

    # cell states
    cs = groups[groups["cell_type_id"].str.startswith("visionkb:cs-")]
    nodes = []
    cs_edges = []
    for cid, g in cs.groupby("cell_type_id", sort=True):
        f = g.iloc[0]
        parent = next((p for p in g["parent_cl"] if isinstance(p, str) and p), None)
        nodes.append(
            {
                "id": cid,
                "category": ["biolink:Cell", "visionkb:CellState"],
                "name": f"{f.group_label} ({f.dataset_id.removeprefix('visionkb:ds-')})",
                "group_label": f.group_label,
                "in_dataset": f.dataset_id,
                "dataset_id": f.dataset_id,
                "summary_level": sorted(set(g.level)),
                "parent_cell_type": parent,
                "mapping_predicate": next((p for p in g["mapping_predicate"] if p), None),
                "mapping_source": sorted(set(g.mapping_source)),
                "in_taxon": HUMAN,
                "provided_by": [SOURCE],
            }
        )
        if parent:
            how = set(g.mapping_source)
            cs_edges.append(
                make_edge(
                    cid,
                    "biolink:subclass_of",
                    parent,
                    source=SOURCE,
                    knowledge_level="knowledge_assertion",
                    agent_type="manual_agent" if "crosswalk" in how else "automated_agent",
                    release=release,
                    visibility=vis[f.dataset_id],
                )
            )
    stats.update(
        study_results=len(srs),
        expressed_in_edges=len(edges),
        cell_states=len(nodes),
        cell_state_subclass_edges=len(cs_edges),
        cell_states_without_parent=sorted(n["id"] for n in nodes if not n.get("parent_cell_type")),
        rows_by_level=rows["level"].value_counts().sort_index().to_dict(),
        edges_to_cellstate=int(sum(e["object"].startswith("visionkb:cs-") for e in edges)),
        genes=int(rows["gene_id"].nunique()),
        max_evidence_count=int(max((e["evidence_count"] for e in edges), default=0)),
        median_evidence_count=float(np.median([e["evidence_count"] for e in edges])) if edges else 0.0,
    )
    return nodes + srs, edges + cs_edges, stats
