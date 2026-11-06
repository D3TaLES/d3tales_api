# D3TaLES API pytest infrastructure

## Goals

This suite converts the legacy executable scripts into isolated pytest tests. It separates fast unit tests from parser integration, database, network, and slow scientific regression tests.

## Install and run

```bash
pip install -e .
pip install pytest pytest-cov pytest-mock pytest-timeout
pytest
pytest -m "not database and not network and not slow"
pytest --cov=d3tales_api --cov-report=term-missing
```

## Migration rules

1. Replace prints with assertions.
2. Replace current-working-directory paths with `data_dir` fixtures.
3. Use `pytest.approx` for floating-point scientific outputs.
4. Do not contact the production REST API or MongoDB in the default suite.
5. Mock database classes at the import location used by the code under test.
6. Add one valid case, boundary cases, malformed input, missing metadata, and unit-conversion variants for every calculator/parser.
7. Mark real database/network tests and require explicit opt-in.

## Recommended test matrix

- Calculators: nominal equations, unit conversion, connector traversal, malformed units, missing keys, zero/negative inputs, non-mutation, numerical tolerances.
- Processors: each supported file flavor, metadata precedence, path handling, malformed/truncated files, schema output, deterministic IDs/timestamps through mocks.
- Database: validation, insert/update behavior, duplicate IDs, query construction, auth failures, no production connections.
- Back-to-front: computed property mapping, conditions provenance, backend ID linkage, idempotency, insert=False behavior.
- Plotters: returned trace structure and labels; avoid pixel-perfect image tests.

## Legacy script mapping

- `concentration_test.py` -> `tests/unit/calculators/test_concentration.py`
- `cv_test.py`, `diffusion_test.py`, `static_diffusion_test.py` -> `tests/unit/calculators/test_cv.py`
- `dft_test.py` -> `tests/unit/calculators/test_dft.py`
- parser scripts -> `tests/integration/processors/` once their raw fixtures are copied
- database insertion scripts -> opt-in database tests using an isolated database or mocks
