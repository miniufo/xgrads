# Contributing to xgrads

## Development setup

Create a Python 3.11 or newer environment and install the package with its
test dependencies:

```bash
python -m pip install -e ".[test]"
```

Run the tests and lint checks before opening a pull request:

```bash
python -m pytest --cov --cov-report=term-missing
ruff check .
```

Documentation can be built locally with:

```bash
python -m pip install -r docs/requirements.txt
python -m sphinx -W --keep-going -b html docs/source docs/_build/html
```

## Pull requests

Please include tests for parser, dataset, projection, or interpolation changes.
Keep changes focused and explain compatibility implications for supported GrADS
descriptor syntax.

## Release checklist

1. Update `xgrads/__init__.py` and `CHANGELOG.md`.
2. Run the full test and documentation builds.
3. Build and check the wheel and source distribution.
4. Create a Git tag and GitHub Release.
5. Verify installation from PyPI and update the conda-forge recipe if needed.
