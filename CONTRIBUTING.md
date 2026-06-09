# Contributing to weatherisk

Thank you for your interest in `weatherisk`. Contributions, bug reports, and suggestions are welcome.

## Reporting issues

Please use [GitHub Issues](https://github.com/wkchikach99-byte/weatherisk/issues) to report bugs or request features. Include:

- Your operating system and Python version (`python --version`)
- The full error traceback
- A minimal reproducible example if possible

## Development setup

```bash
git clone https://github.com/wkchikach99-byte/weatherisk.git
cd weatherisk
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
```

## Running the tests

```bash
pytest                             # all tests
pytest -m "not slow"               # skip long-running tests
pytest tests/test_clustering.py    # single module
```

## Code style

The project uses [Ruff](https://docs.astral.sh/ruff/) for linting and [mypy](https://mypy.readthedocs.io) for type checking:

```bash
ruff check weatherisk/
mypy weatherisk/
```

Please ensure both pass before submitting a pull request.

## Submitting a pull request

1. Fork the repository and create a feature branch from `main`.
2. Make your changes with a clear, focused commit message.
3. Add or update tests for any changed behaviour.
4. Open a pull request describing what was changed and why.

## Citation

If you use `weatherisk` in your research, please cite the associated SoftwareX paper and the Zenodo archive — see [CITATION.cff](CITATION.cff) for the full reference.
