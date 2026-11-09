import pandas as pd

from visionkb_kg import known_issues as KI


def _summary(tmp_path):
    ds = "visionkb:ds-t-snrna"
    base = dict(dataset_id=ds, atlas_version="t-0", n_donors=3, mean_log1p_cp10k=1.0, n_donors_expressing=3)
    rod = [f"R{i}" for i in range(5)]
    ac = [f"A{i}" for i in range(5)]
    rows = []

    def add(level, group, n, gene, pb):
        rows.append(
            dict(
                base,
                level=level,
                group_label=group,
                cell_type_id="CL:0000001",
                gene_id=f"HGNC:{abs(hash(gene)) % 10**6}",
                gene_symbol=gene,
                var_symbol=gene,
                n_cells=n,
                pct_expressing=0.9 if pb >= 1 else 0.05,
                pseudobulk_cp10k=pb,
            )
        )

    for g in rod + ac:
        add("majorclass", "Rod", 1000, g, 50.0 if g in rod else 0.05)
        add("majorclass", "AC", 500, g, 40.0 if g in ac else 0.1)
        # AC1 carries 30% of the rod program (contaminated); AC2/AC3 are clean
        add("celltype", "AC1", 60, g, 15.0 if g in rod else 40.0)
        add("celltype", "AC2", 300, g, 0.2 if g in rod else 40.0)
        add("celltype", "AC3", 200, g, 0.1 if g in rod else 40.0)
        add("celltype", "Rod1", 900, g, 50.0 if g in rod else 0.05)
    p = tmp_path / "s.parquet"
    pd.DataFrame(rows)[KI.heca.COLUMNS].to_parquet(p)
    groups = pd.DataFrame(
        dict(
            dataset_id=ds,
            level="celltype",
            group_label=["AC1", "AC2", "AC3", "Rod1"],
            majorclass=["AC", "AC", "AC", "Rod"],
            n_cells=[60, 300, 200, 900],
        )
    )
    return p, groups, ds


def test_cross_lineage_flags_contaminated_subtype(tmp_path):
    p, groups, ds = _summary(tmp_path)
    ki = KI.cross_lineage(p, groups, {ds: True})
    pairs = {(x["group_label"], x["program_of"]): x for x in ki["pairs"]}
    assert set(pairs) == {("AC1", "Rod")}
    x = pairs[("AC1", "Rod")]
    assert abs(x["elevation_median_ratio"] - 0.3) < 1e-6
    # 5 rod markers in AC1: pct 0.9, 60 cells, 15 >= 0.1 * 50, 3 donors -> all kept by D18
    assert x["d18_study_results_to_those_markers"] == 5
    assert ki["n_study_results"] == 5
    md = KI.markdown(ki)
    assert "## Known issues" in md and "AC1 (AC)" in md
