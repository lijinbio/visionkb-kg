# visionkb-kg

LinkML schema on Biolink 4.4.6, per-source ingests to KGX, validation and release bundles for the VisionKB KG.

- `docs/contracts.md` (sanitized project contracts) defines the release layout, Neo4j layout, ids and provenance; change it first.
- Biolink names only (`biolink:expressed_in`, ...); no custom relationship names. bmt must load the vendored 4.4.6 (`visionkb_kg.biolink.toolkit`).
- Validate every id against its source; never guess a mapping (OLS exact label/synonym only; unmapped items are reported).
- Every filter or exclusion is counted in the ingest/validation reports.
- No lab data, personal paths, accounts or tokens in the repo; real inputs come from `$VISIONKB_DATA_ROOT/config/sources.yaml`.
- File and directory names from CURIEs only via `sources.dirname_of`; release files use relative paths only.
- After schema edits: `python scripts/gen_json_schema.py`. Before a commit: `ruff check . && ruff format . && pytest`.
