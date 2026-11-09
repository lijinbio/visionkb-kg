from visionkb_kg.ids import check_curie, edge_id, short_id
from visionkb_kg.obo import ancestors, descendants, read_obo

OBO = """format-version: 1.2
data-version: test/2026-01-01

[Term]
id: UBERON:0000001
name: eye

[Term]
id: UBERON:0000002
name: retina
synonym: "inner layer of eyeball" EXACT []
is_a: UBERON:0000009 ! tissue
relationship: part_of UBERON:0000001 ! eye

[Term]
id: UBERON:0000003
name: fovea
relationship: part_of UBERON:0000002 ! retina

[Term]
id: UBERON:0000009
name: tissue

[Term]
id: UBERON:0000004
name: old
is_obsolete: true

[Typedef]
id: part_of
"""


def test_edge_id_deterministic():
    a = edge_id("HGNC:1", "biolink:expressed_in", "CL:1", {"anatomical_context_qualifier": ["UBERON:2", "UBERON:1"]})
    b = edge_id(
        "HGNC:1",
        "biolink:expressed_in",
        "CL:1",
        {"anatomical_context_qualifier": ["UBERON:1", "UBERON:2"], "stage_qualifier": None},
    )
    assert a == b and a.startswith("visionkb:e-")
    assert a != edge_id("HGNC:1", "biolink:expressed_in", "CL:1")
    assert short_id("sr", 1, 2) == short_id("sr", 1, 2)


def test_check_curie():
    assert check_curie("HGNC:16665") is None
    assert check_curie("HGNC:APLN") == "malformed HGNC id"
    assert check_curie("CL:604") == "malformed CL id"
    assert check_curie("no colon") == "not a CURIE"


def test_obo(tmp_path):
    p = tmp_path / "t.obo"
    p.write_text(OBO)
    header, terms = read_obo(p)
    assert header["data-version"] == "test/2026-01-01"
    assert terms["UBERON:0000002"]["exact_synonyms"] == ["inner layer of eyeball"]
    assert terms["UBERON:0000004"]["obsolete"]
    assert descendants(terms, "UBERON:0000001") == {"UBERON:0000001", "UBERON:0000002", "UBERON:0000003"}
    assert ancestors(terms, ["UBERON:0000003"]) == {
        "UBERON:0000003",
        "UBERON:0000002",
        "UBERON:0000001",
        "UBERON:0000009",
    }
