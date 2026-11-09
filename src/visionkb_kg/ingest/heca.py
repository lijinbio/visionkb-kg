"""HECA atlas summaries -> expressed_in edges (gene -> CL class or CellState) with ExpressionStudyResult nodes.

Edge rule D18 (replaces D8). A (gene, dataset, group) at level majorclass or celltype becomes a StudyResult (and supports an
expressed_in edge) only if all hold: (a) detection: pct_expressing >= 0.10 and n_cells >= 20; (b) specificity:
pseudobulk_cp10k >= 0.10 x the gene's maximum pseudobulk_cp10k over the groups of the same dataset and level; (c) donor
support: n_donors_expressing >= 2 (only for atlases with a usable donor column; otherwise (a)+(b) and the edge is flagged
donor_support_checked=false). region_majorclass rows stay in the Parquet only. One edge per (gene, cell type, dataset
tissues, species); edge summaries come from the surviving StudyResults. CellState nodes come from groups.parquet and are
subclass_of their CL class.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as pds

from ..ids import edge_id, short_id
from .common import HUMAN, make_edge

SOURCE = "infores:visionkb-heca"
LEVELS = ("majorclass", "celltype")
RULE = {
    "min_pct_expressing": 0.10,
    "min_n_cells": 20,
    "specificity_fraction": 0.10,
    "min_n_donors_expressing": 2,
    "levels": list(LEVELS),
    "decision": "D18",
}
ROBUSTNESS_FRACTIONS = (0.05, 0.10, 0.20)
SANITY_GENES = ("RHO", "PDE6A", "OPN1SW", "RLBP1", "GLUL", "RPE65", "PAX6", "ACTB")
COLUMNS = [
    "dataset_id",
    "atlas_version",
    "level",
    "group_label",
    "cell_type_id",
    "gene_id",
    "gene_symbol",
    "var_symbol",
    "n_cells",
    "n_donors",
    "pct_expressing",
    "mean_log1p_cp10k",
    "pseudobulk_cp10k",
    "n_donors_expressing",
]


def method_string(rule: dict = RULE) -> str:
    return (
        f"visionkb-atlas streamed summary; {rule['decision']}: pct_expressing>={rule['min_pct_expressing']}, "
        f"n_cells>={rule['min_n_cells']}, pseudobulk_cp10k>={rule['specificity_fraction']}*max over groups of the same dataset and level, "
        f"n_donors_expressing>={rule['min_n_donors_expressing']}; levels {'/'.join(rule['levels'])}"
    )


def small_group_report(df, a, b, c, hgnc, donor_cells: pd.DataFrame, small_n: int = 500, top: int = 10) -> dict:
    """Groups with n_cells < small_n: cells per donor, rows removed per rule, top retained genes and genes removed only by (a)
    or only by (c); plus, over all groups, rows removed only by (c) where the group's median cells per donor is < 10."""
    dc = donor_cells[donor_cells["level"].isin(LEVELS)]
    per = (
        dc.groupby(["dataset_id", "level", "group_label"])["n_cells"]
        .agg(["sum", "count", "min", "median"])
        .reset_index()
    )
    per.columns = [
        "dataset_id",
        "level",
        "group_label",
        "n_cells",
        "n_donors",
        "min_cells_per_donor",
        "median_cells_per_donor",
    ]
    key = ["dataset_id", "level", "group_label"]
    f = pd.DataFrame({k: df[k].astype(str) for k in key})
    f["gene"] = df["gene_symbol"].astype(str)
    f["pb"] = df["pseudobulk_cp10k"].to_numpy()
    f["a"], f["b"], f["c"], f["hgnc"] = a.to_numpy(), b, c.to_numpy(), hgnc.to_numpy()
    f = f[f.hgnc]
    f["keep"] = f.a & f.b & f.c
    f["only_c"] = f.a & f.b & ~f.c
    f["only_a"] = ~f.a & f.b & f.c
    m = f.merge(per[[*key, "median_cells_per_donor"]], on=key, how="left")
    low = m[m.only_c & (m.median_cells_per_donor < 10)]
    out = {
        "threshold_n_cells": small_n,
        "only_c_rows_in_groups_with_median_cells_per_donor_lt_10": len(low),
        "groups_with_median_cells_per_donor_lt_10": int((per.median_cells_per_donor < 10).sum()),
        "groups_affected_by_only_c_with_median_lt_10": int(low[key].drop_duplicates().shape[0]),
        "only_c_rows_total": int(f.only_c.sum()),
        "groups": [],
    }
    small = per[per.n_cells < small_n].sort_values(["n_cells", "dataset_id", "group_label"])
    for g in small.itertuples(index=False):
        sub = f[(f.dataset_id == g.dataset_id) & (f.level == g.level) & (f.group_label == g.group_label)]
        out["groups"].append(
            {
                "dataset_id": g.dataset_id,
                "level": g.level,
                "group_label": g.group_label,
                "n_cells": int(g.n_cells),
                "n_donors": int(g.n_donors),
                "min_cells_per_donor": int(g.min_cells_per_donor),
                "median_cells_per_donor": float(g.median_cells_per_donor),
                "candidate_rows": len(sub),
                "kept": int(sub.keep.sum()),
                "fail_a": int((~sub.a).sum()),
                "fail_b": int((~sub.b).sum()),
                "fail_c": int((~sub.c).sum()),
                "only_a": int(sub.only_a.sum()),
                "only_c": int(sub.only_c.sum()),
                "top_kept": sub[sub.keep].nlargest(top, "pb")["gene"].tolist(),
                "top_removed_only_c": sub[sub.only_c].nlargest(top, "pb")["gene"].tolist(),
                "top_removed_only_a": sub[sub.only_a].nlargest(top, "pb")["gene"].tolist(),
            }
        )
    return out


def load_levels(summary_files: list[Path], levels=LEVELS) -> pd.DataFrame:
    """All gene x group rows of the graph levels (no threshold), strings as categoricals."""
    parts = []
    for f in summary_files:
        t = pds.dataset(str(f)).to_table(columns=COLUMNS, filter=pds.field("level").isin(list(levels)))
        parts.append(t.to_pandas(strings_to_categorical=True))
    df = pd.concat(parts, ignore_index=True)
    for c in (
        "dataset_id",
        "atlas_version",
        "level",
        "group_label",
        "cell_type_id",
        "gene_id",
        "gene_symbol",
        "var_symbol",
    ):
        df[c] = df[c].astype("category")
    return df


def select_d18(
    df: pd.DataFrame,
    donor_ok: dict[str, bool],
    rule: dict = RULE,
    donor_cells: pd.DataFrame | None = None,
    small_n: int = 500,
) -> tuple[pd.DataFrame, dict]:
    """Flag rules a/b/c on every candidate row; return the kept rows and per dataset x level counts."""
    a = (df["pct_expressing"] >= rule["min_pct_expressing"]) & (df["n_cells"] >= rule["min_n_cells"])
    gmax = df.groupby(["dataset_id", "level", "var_symbol"], observed=True)["pseudobulk_cp10k"].transform("max")
    ratio = np.where(gmax > 0, df["pseudobulk_cp10k"] / gmax.where(gmax > 0, 1.0), 0.0)
    checked = df["dataset_id"].map(lambda d: bool(donor_ok.get(d, False))).astype(bool)
    c = (~checked) | (df["n_donors_expressing"] >= rule["min_n_donors_expressing"])
    hgnc = df["gene_id"].astype(str).str.startswith("HGNC:")
    out = {
        "rule": rule,
        "donor_support_checked": {d: bool(v) for d, v in sorted(donor_ok.items())},
        "by_dataset_level": [],
        "robustness": {},
    }
    b0 = ratio >= rule["specificity_fraction"]
    keys = df[["dataset_id", "level"]].astype(str).agg("|".join, axis=1)
    frame = pd.DataFrame({"key": keys, "a": a, "b": b0, "c": c, "hgnc": hgnc})
    for k, g in frame.groupby("key", sort=True):
        ds, lv = k.split("|")
        cand = g[g.hgnc]
        out["by_dataset_level"].append(
            {
                "dataset_id": ds,
                "level": lv,
                "candidate_rows": len(g),
                "non_hgnc_gene_rows": int((~g.hgnc).sum()),
                "candidates_hgnc": len(cand),
                "fail_a_detection": int((~cand.a).sum()),
                "fail_b_specificity": int((~cand.b).sum()),
                "fail_c_donor_support": int((~cand.c).sum()),
                "fail_any": int((~(cand.a & cand.b & cand.c)).sum()),
                "kept": int((cand.a & cand.b & cand.c).sum()),
                "kept_old_rule_D8": int(cand.a.sum()),
                "donor_support_checked": bool(donor_ok.get(ds, False)),
            }
        )
    for frac in ROBUSTNESS_FRACTIONS:
        keep = a & (ratio >= frac) & c & hgnc
        sub = df.loc[keep, ["gene_id", "cell_type_id", "dataset_id"]]
        out["robustness"][str(frac)] = {
            "study_results": int(keep.sum()),
            "edges_gene_celltype_dataset_tissue": int(
                sub.drop_duplicates(["gene_id", "cell_type_id", "dataset_id"]).shape[0]
            ),
        }
    keep = a & b0 & c & hgnc
    sel = keep & df["gene_symbol"].isin(list(SANITY_GENES)) & (df["level"] == "majorclass")
    det = df.loc[
        sel, ["gene_symbol", "dataset_id", "group_label", "pct_expressing", "pseudobulk_cp10k", "n_donors_expressing"]
    ].astype({"gene_symbol": str, "dataset_id": str, "group_label": str})
    det["ratio_to_max"] = ratio[sel.to_numpy()]
    sanity_detail = {
        g: [
            f"{r.dataset_id.removeprefix('visionkb:ds-')}:{r.group_label} pct={r.pct_expressing:.2f} pb={r.pseudobulk_cp10k:.3g} "
            f"ratio={r.ratio_to_max:.2f} donors={int(r.n_donors_expressing)}"
            for r in sub.sort_values(["dataset_id", "ratio_to_max"], ascending=[True, False]).itertuples()
        ]
        for g, sub in det.groupby("gene_symbol")
    }
    kept = df.loc[keep].copy()
    kept["donor_support_checked"] = checked[keep].to_numpy()
    for col in (
        "dataset_id",
        "atlas_version",
        "level",
        "group_label",
        "cell_type_id",
        "gene_id",
        "gene_symbol",
        "var_symbol",
    ):
        kept[col] = kept[col].astype(str)
    old = df.loc[a & hgnc, ["dataset_id", "level", "group_label", "gene_symbol"]].astype(str)
    new = kept[["dataset_id", "level", "group_label", "gene_symbol"]]
    sanity = {}
    for gname in SANITY_GENES:
        rec = {}
        for tag, t in (("D8", old), ("D18", new)):
            t = t[(t.gene_symbol == gname) & (t.level == "majorclass")]
            for ds, gg in t.groupby("dataset_id"):
                rec.setdefault(ds.removeprefix("visionkb:ds-"), {})[tag] = sorted(gg.group_label)
        ct_old = old[(old.gene_symbol == gname) & (old.level == "celltype")]
        ct_new = new[(new.gene_symbol == gname) & (new.level == "celltype")]
        sanity[gname] = {
            "majorclass": rec,
            "celltype_groups": {
                "D8": len(ct_old),
                "D18": len(ct_new),
                "D18_labels": sorted(
                    f"{d.removeprefix('visionkb:ds-')}:{x}"
                    for d, x in zip(ct_new.dataset_id, ct_new.group_label, strict=True)
                )[:40],
            },
        }
    for g, v in sanity.items():
        v["majorclass_D18_detail"] = sanity_detail.get(g, [])
    out["sanity"] = sanity
    if donor_cells is not None:
        out["small_groups"] = small_group_report(df, a, b0, c, hgnc, donor_cells, small_n)
    tot = pd.DataFrame(out["by_dataset_level"])
    out["totals"] = {
        k: int(tot[k].sum())
        for k in (
            "candidate_rows",
            "non_hgnc_gene_rows",
            "candidates_hgnc",
            "fail_a_detection",
            "fail_b_specificity",
            "fail_c_donor_support",
            "fail_any",
            "kept",
            "kept_old_rule_D8",
        )
    }
    return kept, out


def build(
    rows: pd.DataFrame, groups: pd.DataFrame, atlases: dict, release: str, method: str
) -> tuple[list, list, dict]:
    """rows: kept rows from select_d18; groups: concatenated groups.parquet; atlases: dataset_id -> {tissue, access_level}."""
    stats = {"study_result_rows": len(rows)}
    rows = rows.reset_index(drop=True)
    rows["ctx"] = [tuple(sorted(atlases[d]["tissue"])) for d in rows["dataset_id"]]
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
                "n_cells": int(r.n_cells),
                "n_donors": int(r.n_donors),
                "pct_expressing": round(float(r.pct_expressing), 6),
                "mean_log1p_cp10k": round(float(r.mean_log1p_cp10k), 6),
                "pseudobulk_cp10k": round(float(r.pseudobulk_cp10k), 6),
                "n_donors_expressing": int(r.n_donors_expressing),
                "donor_support_checked": bool(r.donor_support_checked),
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
        donor_checked=("donor_support_checked", "all"),
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
            donor_support_checked=bool(r.donor_checked),
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
