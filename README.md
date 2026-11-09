# visionkb-kg

Build of the VisionKB knowledge graph: a LinkML extension of the [Biolink Model](https://biolink.github.io/biolink-model/) 4.4.6,
one ingest per source writing [KGX](https://github.com/biolink/kgx) JSONL, validation, and a deterministic release folder for
`neo4j-admin database import full`, packed as a relocatable, checksummed bundle.

## Pieces

| Part | Where | What |
|---|---|---|
| Schema | `src/visionkb_kg/schema/visionkb.yaml` | imports the vendored Biolink 4.4.6 (`schema/biolink/`, byte-identical to the v4.4.6 tag, sha256 checked in tests); adds `CellState` (dataset-specific, `subclass_of` a CL class), `VisionKBDataset` registry slots, `ExpressionStudyResult` / `LiteratureStudyResult`, `GeneExpressedInCellAssociation`, `CuratedGeneToConditionAssociation`, imaging stubs `ImagingModality`, `ImageDerivedFeature`, `RetinalLocation` |
| JSON Schema | `src/visionkb_kg/schema/generated/visionkb.schema.json` | `python scripts/gen_json_schema.py` (definitions reachable from the visionkb classes) |
| Crosswalk | `src/visionkb_kg/mappings/heca-celltypes.sssom.tsv` | SSSOM: atlas labels → Cell Ontology, every CL id checked in OLS (`scripts/build_crosswalk.py`) |
| Ingests | `visionkb-kg ingest <source>` | `hgnc` (approved protein-coding genes + any HGNC id used), `heca` (expressed_in edges + StudyResults + CellStates), `retigene` (curated gene → MONDO/HP edges with PMIDs), `datasets` (Dataset nodes), `ontologies` (CL, UBERON, MONDO, HP slims from OBO releases; EFO, NCBITaxon from OLS) |
| Build | `visionkb-kg build <release>` | merge, validation report, `neo4j-import/`, `summaries/`, `registry/`, `manifest.json`, `RELEASE` |
| Bundle | `visionkb-kg pack <release>`, `visionkb-kg verify <bundle>` | `<release>.tar` + `.sha256`, relative paths only, every file checksummed in `manifest.json` |

Edges carry Biolink provenance natively (`primary_knowledge_source`, `knowledge_level`, `agent_type`); see `docs/contracts.md` §5–§7.

| Source | primary_knowledge_source | knowledge_level | agent_type |
|---|---|---|---|
| atlas expression | `infores:visionkb-heca` | observation | data_analysis_pipeline |
| CellState → CL (crosswalk / majorclass rule) | `infores:visionkb-heca` | knowledge_assertion | manual_agent / automated_agent |
| curated haploinsufficiency review | `infores:visionkb-retigene` | knowledge_assertion | manual_validation_of_automated_agent |
| ontology hierarchy | `infores:cl`, `infores:uberon`, `infores:mondo`, `infores:hpo`, ... | knowledge_assertion | manual_agent |

## Validation

`build` writes `validation/report.json` and `report.md`: every CURIE checked against its source term set (HGNC approved ids and the
symbol of each id; non-obsolete CL/UBERON/MONDO/HP terms of the recorded releases; EFO/NCBITaxon verified in OLS), Biolink categories
and predicates (bmt on the vendored 4.4.6), domain/range of every category-predicate triple, dangling edges, study-result ↔ edge
links, and the KGX validator (streaming). A self-test checks that the design proposal's wrong id for APLN (`HGNC:16693`; APLN is
`HGNC:16665`) is flagged.

`visionkb-kg post-build <release> [--neo4j]` (before `pack`) rewrites `validation/report.md` (small-group tables, Known issues)
from `report.json` / `known_issues.json`, optionally regenerates `neo4j-import/` from the release KGX, and refreshes the manifest
(counts, file checksums, `code.visionkb-kg.git_sha`, a `post_build` record); KGX, summaries and registry are not touched.

Known issues (report only; nothing is flagged or removed) are computed from the release summaries and written to
`validation/known_issues.json`, a "Known issues" section of `report.md` and the manifest `known_issues` list. KI-1
(`src/visionkb_kg/known_issues.py`): per atlas, the top-50 markers of each majorclass are looked up in every cell-type group of
another majorclass; a group is listed when the median ratio of its marker pseudobulk to the majorclass level is >= 0.10 and >= 5x
that of its sibling cell types, with the number of D18 expressed_in study results it contributes for those markers.

## Install and run

```bash
pip install git+https://github.com/lijinbio/visionkb-kg.git@v0.1.0
cp config/sources.example.yaml $VISIONKB_DATA_ROOT/config/sources.yaml   # set real paths
for s in heca retigene datasets hgnc ontologies; do visionkb-kg ingest $s --release visionkb-YYYY.MM.N; done
visionkb-kg build visionkb-YYYY.MM.N
visionkb-kg pack visionkb-YYYY.MM.N && visionkb-kg verify $VISIONKB_DATA_ROOT/kg/bundles/visionkb-YYYY.MM.N.tar
```

Neo4j: run `neo4j-admin database import full neo4j @import.args` from inside `neo4j-import/` (paths there are relative), then
apply `post_import.cypher` (uniqueness constraint on `NamedThing.id`, `StudyResult.edge_id` index, full-text index `entity_search`).

Environment: `VISIONKB_DATA_ROOT` (releases, bundles, config), `VISIONKB_WORK_DIR` or `VISIONKB_KG_WORK` (build working files).

## Tests

`pytest` on synthetic data only (no lab data in this repository): schema and JSON Schema instances, ids, OBO parsing, ingests,
validation (including the APLN check), Neo4j CSVs, pack/verify with tamper detection.

## License

MIT, see [LICENSE](LICENSE). The vendored Biolink Model files are CC0 (Biolink Model authors).
