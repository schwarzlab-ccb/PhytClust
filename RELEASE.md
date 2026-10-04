# Preparing PhytClust 1.0.0

The package version is defined in `pyproject.toml`. Keep `CITATION.cff` and
`CHANGELOG.md` in sync with it before creating a release tag.

Run these checks from the repository root:

```sh
MPLBACKEND=Agg pytest -q
ruff check src/ tests/
mkdocs build --strict
python -m build
python -m twine check dist/phytclust-1.0.0*
```

Install the wheel in a separate environment and check `phytclust --version`,
clustering, and the GUI before uploading. The release files are
`dist/phytclust-1.0.0-py3-none-any.whl` and `dist/phytclust-1.0.0.tar.gz`.

Local Clonetrac inputs, notebooks, generated results, and review figures are
ignored by Git. `MANIFEST.in` also excludes analysis folders from source
archives. The analysis scripts remain in the Git repository for reproducibility.

After reviewing the committed release and passing CI on `master`, create the
`v1.0.0` tag and upload these two files to PyPI. Publishing and tagging are
separate from preparing the release.
