"""Known issues of a release, computed from the atlas summaries (report only; nothing is flagged or dropped).

KI-1 cross-lineage marker programs: per atlas, the top-50 markers of each majorclass M (HGNC genes with pct >= 0.25 and
pseudobulk_cp10k >= 1 in M, ranked by pseudobulk in M / (max pseudobulk in any other majorclass + 0.01)) are looked up in every
cell-type group G of another majorclass. elevation(G, M) = median over the markers of pseudobulk(G) / pseudobulk(M). A pair is
listed when elevation >= 0.10 and >= 5 x the median elevation of G's sibling cell types (same dominant majorclass); the number of
D18 expressed_in study results of G for those markers is counted. Typical causes: doublets or contamination with an abundant
lineage (e.g. rod transcripts in small amacrine subtypes of a retina snRNA atlas).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pyarrow.dataset as pds

from .ingest import heca

N_MARKERS = 50
MIN_PCT = 0.25
MIN_PB = 1.0
FLAG_ELEVATION = 0.10
FLAG_FOLD_OVER_SIBLINGS = 5.0
ROD_PANEL = ("PDE6A", "SAG", "RHO", "NRL")


def _dataset(summary_path, ds: str) -> pd.DataFrame:
    t = pds.dataset(str(summary_path)).to_table(
        columns=heca.COLUMNS, filter=(pds.field("dataset_id") == ds) & pds.field("level").isin(list(heca.LEVELS))
    )
    return t.to_pandas()


def cross_lineage(summary_path, groups: pd.DataFrame, donor_ok: dict[str, bool]) -> dict:
    out = {
        "id": "KI-1",
        "title": "cross-lineage marker programs in cell-type groups",
        "method": __doc__.split("\n\n")[1].replace("\n", " "),
        "parameters": {
            "n_markers": N_MARKERS,
            "min_pct": MIN_PCT,
            "min_pseudobulk_cp10k": MIN_PB,
            "flag_elevation": FLAG_ELEVATION,
            "flag_fold_over_siblings": FLAG_FOLD_OVER_SIBLINGS,
        },
        "pairs": [],
        "rod_panel": {},
    }
    for ds in sorted(groups.loc[groups.level == "celltype", "dataset_id"].unique()):
        df = _dataset(summary_path, ds)
        df = df[df["gene_id"].str.startswith("HGNC:")]
        mc = df[df.level == "majorclass"]
        pb = mc.pivot_table(
            index="gene_symbol", columns="group_label", values="pseudobulk_cp10k", aggfunc="first"
        ).fillna(0.0)
        pct = mc.pivot_table(
            index="gene_symbol", columns="group_label", values="pct_expressing", aggfunc="first"
        ).fillna(0.0)
        ct = df[df.level == "celltype"]
        cpb = ct.pivot_table(
            index="gene_symbol", columns="group_label", values="pseudobulk_cp10k", aggfunc="first"
        ).fillna(0.0)
        g = groups[(groups.dataset_id == ds) & (groups.level == "celltype")].set_index("group_label")
        dom = g["majorclass"].to_dict()
        ncells = g["n_cells"].to_dict()
        kept, _ = heca.select_d18(
            df.astype({c: "category" for c in ("dataset_id", "level", "group_label", "gene_symbol", "var_symbol")}),
            {ds: donor_ok.get(ds, True)},
        )
        kept_ct = kept[kept.level == "celltype"][["group_label", "gene_symbol"]]
        markers = {}
        for m in pb.columns:
            others = pb.drop(columns=[m]).max(axis=1) if pb.shape[1] > 1 else 0.0
            score = pb[m] / (others + 0.01)
            elig = (pct[m] >= MIN_PCT) & (pb[m] >= MIN_PB)
            markers[m] = score[elig].sort_values(ascending=False).index[:N_MARKERS].tolist()
        elev = {}
        for gl in cpb.columns:
            for m, mk in markers.items():
                if not mk or dom.get(gl) == m:
                    continue
                elev[(gl, m)] = float(np.median(cpb.loc[mk, gl].to_numpy() / pb.loc[mk, m].to_numpy()))
        for (gl, m), e in sorted(elev.items()):
            sib = [elev[(x, m)] for x in cpb.columns if dom.get(x) == dom.get(gl) and x != gl and (x, m) in elev]
            base = float(np.median(sib)) if sib else 0.0
            if e >= FLAG_ELEVATION and e >= FLAG_FOLD_OVER_SIBLINGS * max(base, 1e-9):
                mk = markers[m]
                n_edges = int(kept_ct[(kept_ct.group_label == gl) & kept_ct.gene_symbol.isin(mk)].shape[0])
                top = sorted(mk, key=lambda x: -pb.loc[x, m])[:4]
                out["pairs"].append(
                    {
                        "dataset_id": ds,
                        "group_label": gl,
                        "group_majorclass": dom.get(gl),
                        "n_cells": int(ncells.get(gl, 0)),
                        "program_of": m,
                        "elevation_median_ratio": round(e, 4),
                        "sibling_median_ratio": round(base, 4),
                        "d18_study_results_to_those_markers": n_edges,
                        "n_markers": len(mk),
                        "examples": {x: [round(float(cpb.loc[x, gl]), 3), round(float(pb.loc[x, m]), 3)] for x in top},
                    }
                )
        if "Rod" in pb.columns and all(x in pb.index for x in ROD_PANEL):
            panel = {
                "Rod (majorclass)": {
                    "n_cells": int(mc.loc[mc.group_label == "Rod", "n_cells"].iloc[0]),
                    **{x: round(float(pb.loc[x, "Rod"]), 3) for x in ROD_PANEL},
                }
            }
            ac = [x for x in cpb.columns if dom.get(x) == "AC"]
            for gl in sorted(
                x
                for x in ac
                if any(
                    p["group_label"] == x and p["program_of"] == "Rod" for p in out["pairs"] if p["dataset_id"] == ds
                )
            ) + (["HAC3"] if "HAC3" in ac else []):
                panel[gl] = {
                    "n_cells": int(ncells.get(gl, 0)),
                    **{x: round(float(cpb.loc[x, gl]), 3) for x in ROD_PANEL},
                }
            if ac:
                panel[f"median of {len(ac)} AC cell types"] = {
                    x: round(float(np.median(cpb.loc[x, ac])), 3) for x in ROD_PANEL
                }
            out["rod_panel"][ds] = panel
    out["n_pairs"] = len(out["pairs"])
    out["n_groups"] = len({(p["dataset_id"], p["group_label"]) for p in out["pairs"]})
    out["n_study_results"] = int(sum(p["d18_study_results_to_those_markers"] for p in out["pairs"]))
    return out


def markdown(ki: dict) -> str:
    L = [
        "",
        "## Known issues",
        "",
        f"### {ki['id']}: {ki['title']} (report only; nothing flagged or dropped)",
        "",
        ki["method"],
        "",
        f"{ki['n_pairs']} (group, program) pairs in {ki['n_groups']} cell-type groups; they contribute {ki['n_study_results']:,} D18 study results "
        "(expressed_in evidence) for the markers of the other lineage.",
        "",
        "| dataset | cell type (majorclass) | cells | program of | median ratio | sibling median | D18 study results to those markers | examples: group / program CP10k |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for p in sorted(ki["pairs"], key=lambda x: (x["dataset_id"], -x["elevation_median_ratio"])):
        ex = "; ".join(f"{g} {a}/{b}" for g, (a, b) in p["examples"].items())
        L.append(
            f"| {p['dataset_id'].removeprefix('visionkb:ds-')} | {p['group_label']} ({p['group_majorclass']}) | {p['n_cells']:,} | {p['program_of']} | "
            f"{p['elevation_median_ratio']:.3f} | {p['sibling_median_ratio']:.3f} | {p['d18_study_results_to_those_markers']} | {ex} |"
        )
    for ds, panel in ki["rod_panel"].items():
        L += [
            "",
            f"Rod panel, {ds.removeprefix('visionkb:ds-')} (pseudobulk CP10k):",
            "",
            "| group | cells | " + " | ".join(ROD_PANEL) + " |",
            "|---|---|" + "---|" * len(ROD_PANEL),
        ]
        for gl, v in panel.items():
            L.append(f"| {gl} | {v.get('n_cells', '')} | " + " | ".join(f"{v[x]:g}" for x in ROD_PANEL) + " |")
    return "\n".join(L) + "\n"
