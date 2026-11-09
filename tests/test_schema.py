import json

import jsonschema
import pytest
import yaml

from visionkb_kg.biolink import BIOLINK_SHA256, SCHEMA_DIR, VISIONKB_YAML, biolink_sha256, toolkit

SCHEMA = json.loads((SCHEMA_DIR / "generated" / "visionkb.schema.json").read_text())


def _validator(cls):
    return jsonschema.Draft7Validator({"$ref": f"#/$defs/{cls}", "$defs": SCHEMA["$defs"]})


def test_biolink_pinned():
    assert biolink_sha256() == BIOLINK_SHA256
    assert toolkit().get_model_version() == "4.4.6"


def test_schema_imports_biolink_and_defines_classes():
    s = yaml.safe_load(VISIONKB_YAML.read_text())
    assert "biolink/biolink-model_v4.4.6" in s["imports"]
    for c in (
        "CellState",
        "VisionKBDataset",
        "ExpressionStudyResult",
        "LiteratureStudyResult",
        "GeneExpressedInCellAssociation",
        "ImagingModality",
        "ImageDerivedFeature",
        "RetinalLocation",
    ):
        assert c in s["classes"], c
    assert s["classes"]["CellState"]["is_a"] == "cell"


def test_json_schema_instances():
    _validator("CellState").validate(
        {
            "id": "visionkb:cs-syn-scrna-rod-1",
            "category": ["biolink:Cell", "visionkb:CellState"],
            "in_dataset": "visionkb:ds-syn-scrna",
            "parent_cell_type": "CL:0000604",
            "summary_level": "celltype",
        }
    )
    _validator("ExpressionStudyResult").validate(
        {
            "id": "visionkb:sr-1",
            "category": ["visionkb:ExpressionStudyResult"],
            "edge_id": "visionkb:e-abc",
            "dataset_id": "visionkb:ds-syn-scrna",
            "pct_expressing": 0.5,
            "n_cells": 30,
            "summary_level": "majorclass",
        }
    )
    with pytest.raises(jsonschema.ValidationError):
        _validator("ExpressionStudyResult").validate(
            {"id": "visionkb:sr-1", "edge_id": "x", "dataset_id": "visionkb:ds-x", "pct_expressing": 1.5}
        )
    lit = {
        "id": "visionkb:sr-2",
        "category": ["visionkb:LiteratureStudyResult"],
        "edge_id": "visionkb:e-abc",
        "pmid": "PMID:123",
        "grade": "A",
        "supporting_text": ["x" * 300],
    }
    _validator("LiteratureStudyResult").validate(lit)
    with pytest.raises(jsonschema.ValidationError):
        _validator("LiteratureStudyResult").validate({**lit, "supporting_text": ["x" * 301]})
    with pytest.raises(jsonschema.ValidationError):
        _validator("LiteratureStudyResult").validate({**lit, "pmid": "123"})
    _validator("VisionKBDataset").validate(
        {
            "id": "visionkb:ds-syn-scrna",
            "category": ["visionkb:VisionKBDataset"],
            "species": "NCBITaxon:9606",
            "tissue": ["UBERON:0000966"],
            "assay": "EFO:0008913",
            "access_level": "lab",
            "bind": {"bind_path": "/path/to/x.h5ad", "bind_name": "syn", "bind_role": "reference"},
        }
    )
