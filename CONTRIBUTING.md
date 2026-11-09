# Contributing

```bash
git clone https://github.com/lijinbio/visionkb-kg.git && cd visionkb-kg
uv sync && uv run pytest          # or: pip install -e . pytest ruff
pre-commit install                # optional
```

| Topic | Rule |
| --- | --- |
| Tests | unit tests on synthetic data only |
| Lint | `ruff check .` and `ruff format --check .` (settings in `pyproject.toml`) |
| Schema | edit `src/visionkb_kg/schema/visionkb.yaml`, then `python scripts/gen_json_schema.py` |
| Ids | every id validated against its source; mappings only by exact match, unmapped items reported |
| Hygiene | no personal paths, emails, cluster accounts, tokens or lab data; use `/path/to/...`, `<partition>`, `<account>` |
| Versions | bump `pyproject.toml` and `__init__.py`, add a CHANGELOG line, tag `vX.Y.Z` |
