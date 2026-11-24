# Contributing to MissensExplorer

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,gui]"
pre-commit install          # runs ruff and mypy on every commit
```

## Checks

Every change must pass the same four checks as CI:

```bash
ruff check .                # lint
ruff format --check .       # formatting (ruff format . to fix)
mypy                        # types (configured in pyproject.toml)
pytest                      # offline tests; pytest -m network for the live API tests
```

The architecture (modules, data flow, external services) is described in
[docs/architecture.md](docs/architecture.md), and the model in [docs/model.md](docs/model.md).

## Code style

**Formatting.** `ruff format`, 120 columns. Do not align comments or values by hand.

**Types.** Every function in `src/` is fully annotated, and `mypy` must pass.
- Use built-in generics and unions: `list[str]`, `dict[str, float]`, `X | None`.
- Import `Callable`, `Iterable`, `Iterator`, `Mapping` and `Sequence` from `collections.abc`.
- Use the shared aliases in `missensexplorer._types`:
  - `Row`: one output-table row;
  - `JSON`: raw web-API JSON;
  - `Progress`: the progress callback;
  - `KeyValueCache`: the cache interface.
- Use a `dataclass` (or `TypedDict`) when a dict has a fixed set of keys that callers rely on.
- Use `Any` only where the value really is untyped, such as JSON from a web API.

**Docstrings.** Google style (`Args:`, `Returns:`, `Raises:`). Required on:
- every module;
- every public class, function and method.

Writing them:
- Start with a one-line summary in the imperative mood ("Return …", "Fetch …").
- Add `Args:`/`Returns:` when the names and types do not already say everything: units, allowed values, `None` meaning "missing", side effects (network, cache).
- Private helpers (`_name`) get a one-line docstring when their purpose is not obvious.

**Comments.**
- Explain *why*, not *what*.
- Cite the source (paper, guideline, release) of every threshold or constant, as the existing code does, e.g. `# Pejaver et al., AJHG 2022`.

**Constants.** `UPPER_CASE` at the top of the module, with their source.

**Errors and logging.**
- Network failures raise `ApiError` or `OfflineError` (`api.py`) and are recorded in `SourceStatus`, so one failing service never stops a run.
- Use `log = logging.getLogger(__name__)`.
- Never use a bare `except:`.

**Generated files.** Never edit these by hand; regenerate them with their benchmark script:
- `src/missensexplorer/model_params.py`
- `src/missensexplorer/calibration_params.py`
- `src/missensexplorer/frequency_params.py`

**Language.** Code, comments and docstrings are in English. User-facing text in the desktop app goes through `gui/i18n.py` (English and Spanish).

## Tests

- Offline by default. Network calls are mocked (see `fake_services` in `tests/test_robustness.py`) or replayed from fixtures in `tests/data/`.
- Live tests carry `@pytest.mark.network` and run nightly in CI.
- A change in behaviour needs a test. A refactor must keep the results table identical.

## Benchmarks

`benchmarks/` holds the scripts that train, calibrate and validate the model. Each one:
- states in its module docstring what it reads and writes;
- writes its results to `benchmarks/results/`.

Retraining changes the shipped numbers, so a retrain is its own commit, together with its results.

## Commits

- One logical change per commit.
- The message has a short summary line, then the reason for the change.
