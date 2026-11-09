import pytest

from visionkb_kg.validate import TermSets


@pytest.fixture
def termsets():
    # public HGNC facts: APLN is HGNC:16665; HGNC:16693 is GPR162; RHO is HGNC:10012
    return TermSets(
        {"HGNC:16665": "APLN", "HGNC:16693": "GPR162", "HGNC:10012": "RHO"},
        {"CL": {"CL:0000604", "CL:0000540"}, "UBERON": {"UBERON:0000966"}, "MONDO": {"MONDO:0019200"}, "HP": set()},
    )


@pytest.fixture
def mini_graph():
    nodes = [
        {
            "id": "HGNC:10012",
            "category": ["biolink:Gene"],
            "name": "RHO",
            "symbol": "RHO",
            "in_taxon": "NCBITaxon:9606",
        },
        {"id": "CL:0000604", "category": ["biolink:Cell"], "name": "retinal rod cell"},
        {"id": "CL:0000540", "category": ["biolink:Cell"], "name": "neuron"},
        {"id": "UBERON:0000966", "category": ["biolink:AnatomicalEntity"], "name": "retina"},
        {"id": "NCBITaxon:9606", "category": ["biolink:OrganismTaxon"], "name": "Homo sapiens", "ols_verified": True},
        {"id": "visionkb:ds-syn-scrna", "category": ["biolink:Dataset"], "name": "synthetic"},
        {
            "id": "visionkb:cs-syn-scrna-rod-1",
            "category": ["biolink:Cell", "visionkb:CellState"],
            "name": "Rod_1 (syn-scrna)",
            "parent_cell_type": "CL:0000604",
        },
    ]
    from visionkb_kg.ids import edge_id

    eid = edge_id(
        "HGNC:10012",
        "biolink:expressed_in",
        "CL:0000604",
        {"anatomical_context_qualifier": ["UBERON:0000966"], "species_context_qualifier": "NCBITaxon:9606"},
    )
    nodes.append(
        {
            "id": "visionkb:sr-1",
            "category": ["biolink:StudyResult"],
            "name": "RHO in Rod",
            "edge_id": eid,
            "supporting_data_set": ["visionkb:ds-syn-scrna"],
            "pct_expressing": 0.9,
        }
    )
    edges = [
        {
            "id": eid,
            "subject": "HGNC:10012",
            "predicate": "biolink:expressed_in",
            "object": "CL:0000604",
            "anatomical_context_qualifier": ["UBERON:0000966"],
            "species_context_qualifier": "NCBITaxon:9606",
            "has_study_results": ["visionkb:sr-1"],
            "primary_knowledge_source": "infores:visionkb-heca",
            "knowledge_level": "observation",
            "agent_type": "data_analysis_pipeline",
        },
        {
            "id": edge_id("CL:0000604", "biolink:subclass_of", "CL:0000540"),
            "subject": "CL:0000604",
            "predicate": "biolink:subclass_of",
            "object": "CL:0000540",
            "primary_knowledge_source": "infores:cl",
            "knowledge_level": "knowledge_assertion",
            "agent_type": "manual_agent",
        },
        {
            "id": edge_id("visionkb:cs-syn-scrna-rod-1", "biolink:subclass_of", "CL:0000604"),
            "subject": "visionkb:cs-syn-scrna-rod-1",
            "predicate": "biolink:subclass_of",
            "object": "CL:0000604",
            "primary_knowledge_source": "infores:visionkb-heca",
            "knowledge_level": "knowledge_assertion",
            "agent_type": "manual_agent",
        },
    ]
    return nodes, edges
