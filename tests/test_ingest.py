import pandas as pd

from visionkb_kg.ingest import heca, retigene


def _rows():
    base = dict(
        atlas_version="t-0",
        gene_symbol="RHO",
        var_symbol="RHO",
        n_donors=3,
        n_expressing=10,
        mean_log1p_cp10k=2.0,
        pseudobulk_cp10k=50.0,
        n_donors_expressing=2,
    )
    r = [
        dict(
            base,
            dataset_id="visionkb:ds-a-scrna",
            level="majorclass",
            group_label="Rod",
            cell_type_id="CL:0000604",
            region=None,
            gene_id="HGNC:10012",
            n_cells=100,
            pct_expressing=0.9,
        ),
        dict(
            base,
            dataset_id="visionkb:ds-a-snrna",
            level="majorclass",
            group_label="Rod",
            cell_type_id="CL:0000604",
            region=None,
            gene_id="HGNC:10012",
            n_cells=50,
            pct_expressing=0.5,
        ),
        dict(
            base,
            dataset_id="visionkb:ds-a-scrna",
            level="celltype",
            group_label="Rod_1",
            cell_type_id="visionkb:cs-a-scrna-rod-1",
            region=None,
            gene_id="HGNC:10012",
            n_cells=40,
            pct_expressing=0.8,
        ),
        dict(
            base,
            dataset_id="visionkb:ds-a-scrna",
            level="region_majorclass",
            group_label="fovea|Rod",
            cell_type_id="CL:0000604",
            region="fovea",
            gene_id="HGNC:10012",
            n_cells=30,
            pct_expressing=0.7,
        ),
        dict(
            base,
            dataset_id="visionkb:ds-a-scrna",
            level="region_majorclass",
            group_label="mixed|Rod",
            cell_type_id="CL:0000604",
            region="mixed",
            gene_id="HGNC:10012",
            n_cells=30,
            pct_expressing=0.7,
        ),
        dict(
            base,
            dataset_id="visionkb:ds-a-scrna",
            level="majorclass",
            group_label="Rod",
            cell_type_id="CL:0000604",
            region=None,
            gene_id="visionkb:gene-X",
            n_cells=100,
            pct_expressing=0.9,
        ),
    ]
    return pd.DataFrame(r)


def test_heca_build():
    atl = {
        "visionkb:ds-a-scrna": {
            "tissue": ["UBERON:0000966"],
            "region_uberon": {"fovea": "UBERON:0001786"},
            "access_level": "lab",
        },
        "visionkb:ds-a-snrna": {"tissue": ["UBERON:0000966"], "region_uberon": {}, "access_level": "lab"},
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
    nodes, edges, st = heca.build(_rows(), groups, atl, "visionkb-test.0", "m")
    ex = [e for e in edges if e["predicate"] == "biolink:expressed_in"]
    assert st["rows_non_hgnc_gene_skipped"] == 1 and st["rows_region_without_uberon_skipped"] == 1
    assert len(ex) == 3  # retina rod (2 datasets), CellState, fovea rod
    retina = next(
        e for e in ex if e["object"] == "CL:0000604" and e["anatomical_context_qualifier"] == ["UBERON:0000966"]
    )
    assert retina["n_datasets"] == 2 and retina["evidence_count"] == 2 and retina["max_pct_expressing"] == 0.9
    assert retina["knowledge_level"] == "observation" and retina["agent_type"] == "data_analysis_pipeline"
    srs = {n["id"]: n for n in nodes if n["category"] == ["biolink:StudyResult"]}
    assert all(srs[s]["edge_id"] == retina["id"] for s in retina["has_study_results"])
    cs = [n for n in nodes if "visionkb:CellState" in n["category"]]
    assert cs[0]["parent_cell_type"] == "CL:0000604"
    sub = [e for e in edges if e["predicate"] == "biolink:subclass_of"]
    assert sub[0]["agent_type"] == "automated_agent"


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
