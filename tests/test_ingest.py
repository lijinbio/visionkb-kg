import pandas as pd

from visionkb_kg.ingest import heca, retigene


def _all_rows():
    base = dict(atlas_version="t-0", n_donors=3, mean_log1p_cp10k=2.0)
    r = []
    # dataset a, majorclass: RHO high in Rod, low (but detected) in Cone, absent in MG
    for g, cid, n, pct, pb, nde in (
        ("Rod", "CL:0000604", 100, 0.9, 100.0, 3),
        ("Cone", "CL:0000573", 50, 0.3, 5.0, 3),
        ("MG", "CL:0000636", 80, 0.0, 0.0, 0),
    ):
        r.append(
            dict(
                base,
                dataset_id="visionkb:ds-a-scrna",
                level="majorclass",
                group_label=g,
                cell_type_id=cid,
                gene_id="HGNC:10012",
                gene_symbol="RHO",
                var_symbol="RHO",
                n_cells=n,
                pct_expressing=pct,
                pseudobulk_cp10k=pb,
                n_donors_expressing=nde,
            )
        )
    # one donor only -> fails (c); same group in a dataset without donor column -> kept but flagged
    r.append(
        dict(
            base,
            dataset_id="visionkb:ds-a-scrna",
            level="celltype",
            group_label="Rod_1",
            cell_type_id="visionkb:cs-a-scrna-rod-1",
            gene_id="HGNC:10012",
            gene_symbol="RHO",
            var_symbol="RHO",
            n_cells=40,
            pct_expressing=0.8,
            pseudobulk_cp10k=80.0,
            n_donors_expressing=1,
        )
    )
    r.append(
        dict(
            base,
            dataset_id="visionkb:ds-b-snrna",
            level="majorclass",
            group_label="Rod",
            cell_type_id="CL:0000604",
            gene_id="HGNC:10012",
            gene_symbol="RHO",
            var_symbol="RHO",
            n_cells=50,
            pct_expressing=0.5,
            pseudobulk_cp10k=50.0,
            n_donors_expressing=1,
        )
    )
    r.append(
        dict(
            base,
            dataset_id="visionkb:ds-a-scrna",
            level="majorclass",
            group_label="Rod",
            cell_type_id="CL:0000604",
            gene_id="visionkb:gene-X",
            gene_symbol="X",
            var_symbol="X",
            n_cells=100,
            pct_expressing=0.9,
            pseudobulk_cp10k=9.0,
            n_donors_expressing=3,
        )
    )
    df = pd.DataFrame(r)
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


def test_select_d18_and_build():
    df = _all_rows()
    dcells = pd.DataFrame(
        [
            ("visionkb:ds-a-scrna", "majorclass", g, f"d{i}", n)
            for g, n in (("Rod", 50), ("Cone", 25), ("MG", 40))
            for i in (1, 2)
        ]
        + [
            ("visionkb:ds-a-scrna", "celltype", "Rod_1", "d1", 40),
            ("visionkb:ds-b-snrna", "majorclass", "Rod", "d1", 50),
        ],
        columns=["dataset_id", "level", "group_label", "donor", "n_cells"],
    )
    kept, rep = heca.select_d18(df, {"visionkb:ds-a-scrna": True, "visionkb:ds-b-snrna": False}, donor_cells=dcells)
    sg = {(x["dataset_id"], x["group_label"]): x for x in rep["small_groups"]["groups"]}
    assert sg[("visionkb:ds-a-scrna", "Rod_1")]["only_c"] == 1 and sg[("visionkb:ds-a-scrna", "Rod_1")][
        "top_removed_only_c"
    ] == ["RHO"]
    assert (
        sg[("visionkb:ds-a-scrna", "Rod")]["top_kept"] == ["RHO"]
        and sg[("visionkb:ds-a-scrna", "Cone")]["median_cells_per_donor"] == 25
    )
    got = set(zip(kept.dataset_id, kept.group_label, strict=True))
    assert got == {
        ("visionkb:ds-a-scrna", "Rod"),
        ("visionkb:ds-b-snrna", "Rod"),
    }  # Cone fails (b), MG (a), Rod_1 (c), X non-HGNC
    a_mc = next(
        x for x in rep["by_dataset_level"] if x["dataset_id"] == "visionkb:ds-a-scrna" and x["level"] == "majorclass"
    )
    assert (
        a_mc["fail_a_detection"] == 1
        and a_mc["fail_b_specificity"] == 2
        and a_mc["kept"] == 1
        and a_mc["kept_old_rule_D8"] == 2
    )
    assert a_mc["non_hgnc_gene_rows"] == 1
    assert rep["robustness"]["0.05"]["study_results"] == 3  # Cone passes (b) at 0.05
    assert rep["sanity"]["RHO"]["majorclass"]["a-scrna"] == {"D8": ["Cone", "Rod"], "D18": ["Rod"]}
    atl = {
        "visionkb:ds-a-scrna": {"tissue": ["UBERON:0000966"], "access_level": "lab"},
        "visionkb:ds-b-snrna": {"tissue": ["UBERON:0000966"], "access_level": "lab"},
    }
    groups = pd.DataFrame(
        [
            {
                "dataset_id": "visionkb:ds-a-scrna",
                "level": "celltype",
                "group_label": "Rod_1",
                "cell_type_id": "visionkb:cs-a-scrna-rod-1",
                "parent_cl": "CL:0000604",
                "mapping_predicate": "",
                "mapping_source": "majorclass_parent",
            }
        ]
    )
    nodes, edges, _ = heca.build(kept, groups, atl, "visionkb-2026.11.0", heca.method_string())
    ex = [e for e in edges if e["predicate"] == "biolink:expressed_in"]
    assert (
        len(ex) == 1
        and ex[0]["n_datasets"] == 2
        and ex[0]["evidence_count"] == 2
        and ex[0]["max_pct_expressing"] == 0.9
    )
    assert ex[0]["donor_support_checked"] is False  # one supporting dataset had no donor column
    srs = {n["id"]: n for n in nodes if n["category"] == ["biolink:StudyResult"]}
    assert {srs[s]["donor_support_checked"] for s in ex[0]["has_study_results"]} == {True, False}
    assert "D18" in next(iter(srs.values()))["method"]
    assert next(e for e in edges if e["predicate"] == "biolink:subclass_of")["agent_type"] == "automated_agent"


def test_retigene_build():
    genes = pd.DataFrame(
        [
            {
                "gene": "GENEA",
                "review_grade": "A",
                "phenotype": "retinitis pigmentosa",
                "pmid": 1,
                "review_note": "explicit",
            },
            {"gene": "GENEB", "review_grade": "B", "phenotype": "x", "pmid": 2, "review_note": "hedged"},
            {"gene": "GENEN", "review_grade": "N", "phenotype": "y", "pmid": 3, "review_note": "no"},
            {"gene": "GENEU", "review_grade": "A", "phenotype": "z", "pmid": 4, "review_note": "no map"},
        ]
    )
    ev = pd.DataFrame(
        [
            {
                "gene": "GENEA",
                "is_reference_quote": "yes",
                "dosage_wording": "haploinsufficiency",
                "quote": "w " * 400,
                "pmid": 1,
                "sentence_id": "s1",
                "page": 2,
            },
            {
                "gene": "GENEA",
                "is_reference_quote": "no",
                "dosage_wording": "",
                "quote": "ctx",
                "pmid": 9,
                "sentence_id": "s2",
                "page": 3,
            },
            {
                "gene": "GENEB",
                "is_reference_quote": "yes",
                "dosage_wording": "dosage",
                "quote": "q",
                "pmid": 2,
                "sentence_id": "s3",
                "page": 1,
            },
        ]
    )
    mapping = pd.DataFrame(
        [
            {
                "gene": "GENEA",
                "status": "mapped",
                "object_id": "MONDO:0019200",
                "confidence": "medium",
                "match": "label",
            },
            {"gene": "GENEB", "status": "mapped", "object_id": "HP:0000001", "confidence": "broad", "match": "label"},
            {"gene": "GENEU", "status": "no_exact_match", "object_id": None, "confidence": "high", "match": None},
        ]
    )
    sym = {"GENEA": "HGNC:1", "GENEB": "HGNC:2", "GENEN": "HGNC:3", "GENEU": "HGNC:4"}
    nodes, edges, st = retigene.build(genes, ev, mapping, sym, "visionkb-test.0")
    by = {e["subject"]: e for e in edges}
    assert by["HGNC:1"]["predicate"] == "biolink:causes" and by["HGNC:2"]["predicate"] == "biolink:has_phenotype"
    assert "HGNC:3" not in by and "HGNC:4" not in by
    assert by["HGNC:1"]["publications"] == ["PMID:1"]
    assert (
        by["HGNC:1"]["knowledge_level"] == "knowledge_assertion"
        and by["HGNC:1"]["agent_type"] == "manual_validation_of_automated_agent"
    )
    assert all(len(t) <= 300 for n in nodes for t in n.get("supporting_text") or [])
    assert {f["gene"] for f in st["flagged"]} == {"GENEN", "GENEU"}
