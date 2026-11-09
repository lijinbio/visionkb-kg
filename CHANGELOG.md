# Changelog

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow [Semantic Versioning](https://semver.org/).

## [0.1.0] - 2026-11-09

### Added

- LinkML schema `visionkb` 0.1.0 importing a vendored, checksummed Biolink Model 4.4.6; pruned JSON Schema.
- SSSOM crosswalk of HECA atlas labels to the Cell Ontology (OLS-checked).
- Ingests: HGNC, ontology slims (CL, UBERON eye branch, MONDO, HP; EFO and NCBITaxon via OLS), atlas expression
  (expressed_in + StudyResult + CellState), curated haploinsufficiency review (gene → MONDO/HP), datasets.
- Validation (term sets, HGNC symbol/id consistency, bmt categories/predicates/domain-range, dangling edges, KGX validator, APLN self-test).
- Release build: merged KGX, neo4j-admin import CSVs + `import.args` + `post_import.cypher`, summaries, registry, manifest.
- Portable bundle `pack` / `verify` (relative paths, per-file sha256, tamper detection); directory names via `dirname_of`.
