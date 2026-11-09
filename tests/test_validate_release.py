import json
import tarfile

from visionkb_kg.biolink import toolkit
from visionkb_kg.kgx import write_jsonl
from visionkb_kg.neo4j_csv import write_import
from visionkb_kg.release import file_list, pack, verify
from visionkb_kg.validate import apln_self_test, validate_release


def _write(tmp_path, nodes, edges):
    write_jsonl(nodes, tmp_path / "nodes.jsonl")
    write_jsonl(edges, tmp_path / "edges.jsonl")
    return tmp_path / "nodes.jsonl", tmp_path / "edges.jsonl"


def test_validate_clean_graph(tmp_path, mini_graph, termsets):
    n, e = _write(tmp_path, *mini_graph)
    rep = validate_release(n, e, termsets, run_kgx=False, log=lambda *_: None)
    assert rep["n_errors"] == 0, rep["errors"]
    assert rep["status"] == "pass"


def test_validate_flags_problems(tmp_path, mini_graph, termsets):
    nodes, edges = mini_graph
    nodes = [
        *nodes,
        {"id": "HGNC:16693", "category": ["biolink:Gene"], "name": "APLN", "symbol": "APLN"},
        {"id": "CL:9999999", "category": ["biolink:Cell"], "name": "made up"},
    ]
    edges = [
        *edges,
        {
            "id": "visionkb:e-x",
            "subject": "HGNC:10012",
            "predicate": "biolink:expressed_in",
            "object": "CL:0000999",
            "primary_knowledge_source": "infores:x",
            "knowledge_level": "observation",
            "agent_type": "data_analysis_pipeline",
        },
    ]
    n, e = _write(tmp_path, nodes, edges)
    rep = validate_release(n, e, termsets, run_kgx=False, log=lambda *_: None)
    assert rep["status"] == "fail"
    assert set(rep["errors"]) >= {"hgnc_symbol_mismatch", "ontology_id_unknown_or_obsolete", "dangling_object"}


def test_apln_self_test(termsets):
    r = apln_self_test(termsets)
    assert r["flagged"] and r["correct_id_passes"] and r["status"] == "pass"


def test_kgx_validator_runs(tmp_path, mini_graph, termsets):
    n, e = _write(tmp_path, *mini_graph)
    rep = validate_release(n, e, termsets, run_kgx=True, log=lambda *_: None)
    assert rep["kgx"]["biolink_version"] == "4.4.6"


def test_neo4j_pack_verify(tmp_path, mini_graph):
    rel = tmp_path / "rel" / "visionkb-test.0"
    (rel / "kgx").mkdir(parents=True)
    n, e = _write(rel / "kgx", *mini_graph)
    res = write_import(n, e, rel / "neo4j-import", toolkit())
    assert res["counts"]["nodes"]["Gene"] == 1 and res["counts"]["relationships"]["supporting_data_set"] == 1
    head = (rel / "neo4j-import" / "nodes_Gene.csv").read_text().splitlines()[0]
    cols = [c.split(":")[0] if not c.startswith(":") else c for c in head.split(",")]
    assert head.startswith("id:ID,:LABEL,") and len(cols) == len(set(cols)), (
        head
    )  # neo4j-admin rejects duplicate headers
    assert "synonym" not in head or "synonym:string[]" in head
    gene = (rel / "neo4j-import" / "nodes_Gene.csv").read_text().splitlines()[1]
    assert "Gene;BiologicalEntity;NamedThing" in gene or ";NamedThing" in gene
    cs = (rel / "neo4j-import" / "nodes_CellState.csv").read_text()
    assert "CellState;Cell;" in cs
    args = (rel / "neo4j-import" / "import.args").read_text()
    assert "--nodes=./nodes_Gene.csv" in args
    assert all(x.split("=", 1)[1].startswith("./") for x in args.split() if x.startswith(("--nodes", "--rel")))
    assert "entity_search" in (rel / "neo4j-import" / "post_import.cypher").read_text()
    from visionkb_kg.release import release_counts

    cnt = release_counts(rel)
    assert cnt["nodes_total"] == len(mini_graph[0]) and sum(cnt["nodes_by_category"].values()) == cnt["nodes_total"]
    assert cnt["nodes_by_category"]["visionkb:CellState"] == 1 and cnt["edges_by_predicate"]["supporting_data_set"] == 1
    assert cnt["edges_total"] == len(mini_graph[1]) + 1 and cnt["study_results"] == 1
    (rel / "RELEASE").write_text("visionkb-test.0\n")
    (rel / "manifest.json").write_text(json.dumps({"kg_release": "visionkb-test.0", "files": file_list(rel)}))
    out = pack(rel, tmp_path / "bundles")
    v = verify(tmp_path / "bundles" / "visionkb-test.0.tar", tmp_path / "unpacked")
    assert v["status"] == "pass", v
    assert (tmp_path / "unpacked" / "visionkb-test.0" / "neo4j-import" / "import.args").exists()
    with tarfile.open(out["bundle"]) as tf:
        assert all(not m.name.startswith("/") for m in tf.getmembers())
    # tampering is detected
    bad = tmp_path / "bundles" / "visionkb-test.0.tar"
    data = bytearray(bad.read_bytes())
    i = data.find(b"RHO")
    data[i : i + 3] = b"RHX"
    bad.write_bytes(bytes(data))
    assert verify(bad)["status"] == "fail"


def test_check_relative(tmp_path):
    from visionkb_kg.release import check_relative

    (tmp_path / "neo4j-import").mkdir()
    (tmp_path / "manifest.json").write_text(
        json.dumps({"sources": [{"path_or_url": "https://example.org/x"}, {"path_or_url": "summaries/a.parquet"}]})
    )
    (tmp_path / "neo4j-import" / "import.args").write_text(
        "--array-delimiter=;\n--nodes=nodes_Gene.csv\n--relationships=rels_x.csv\n"
    )
    assert check_relative(tmp_path) == []
    (tmp_path / "manifest.json").write_text(json.dumps({"sources": [{"path_or_url": "/abs/data/x.txt"}]}))
    (tmp_path / "neo4j-import" / "import.args").write_text("--nodes=Gene=/abs/nodes.csv\n--relationships=../rels.csv\n")
    assert len(check_relative(tmp_path)) == 3


def _prop(col: str) -> str:
    """Property a neo4j-admin header column writes: 'id:ID' -> 'id', 'name:string' -> 'name'; ':LABEL' etc. -> itself."""
    return col if col.startswith(":") else col.split(":", 1)[0]


def test_neo4j_headers_have_unique_property_names(tmp_path, mini_graph):
    rel = tmp_path / "rel" / "visionkb-test.0"
    (rel / "kgx").mkdir(parents=True)
    n, e = _write(rel / "kgx", *mini_graph)
    write_import(n, e, rel / "neo4j-import", toolkit())
    files = sorted((rel / "neo4j-import").glob("*.csv"))
    assert any(f.name.startswith("nodes_") for f in files) and any(f.name.startswith("rels_") for f in files)
    for f in files:
        cols = f.read_text().splitlines()[0].split(",")
        props = [_prop(c) for c in cols]
        dup = sorted({x for x in props if props.count(x) > 1})
        assert not dup, f"{f.name}: duplicate header names {dup} in {cols}"
        if f.name.startswith("nodes_"):
            assert cols[:2] == ["id:ID", ":LABEL"], cols


def test_prop_helper_rejects_id_vs_id_ID():
    props = [_prop(c) for c in ["id:ID", ":LABEL", "id", "name:string"]]
    assert props.count("id") == 2
