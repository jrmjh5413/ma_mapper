# Releasing

What a release involves, in order, and which steps need an account only the
maintainer has. Everything above the line is already automated or written; the
numbered steps are what is left to do by hand.

## Before the first release: decisions

| decision | why it blocks the release |
|---|---|
| **The name** | *Bioinformatics* requires the tool's name in the title, and the PyPI and Bioconda package names are set at first upload. Renaming after release means a second package and a redirect. |
| **The author list and its order** | goes in `CITATION.cff`, `pyproject.toml`, `.zenodo.json` and the paper. |
| **Which repository is canonical** | the paper cites one URL that has to stay valid for at least two years. Today the original is `pkncsk/ma_mapper` and the work is in `jrmjh5413/ma_mapper`. |

## Release steps

1. **Bump the version** in `ma_mapper/__init__.py` (`__version__`). Everything
   else reads it from there: `pyproject.toml` via `dynamic = ["version"]`, the
   `--version` flag, `ma-mapper info`, and every run record.

2. **Update `CHANGELOG.md`**, and `docs/CHANGES_FROM_UPSTREAM.md` if the fork's
   relationship to the original has changed.

3. **Check the tree is clean and the tests pass** on a supported Python:

   ```bash
   pytest -q                     # 111 tests
   ruff check ma_mapper tests
   python examples/toy_end_to_end.py /tmp/release_check
   ```

4. **Build and validate both distributions.** Write to project space, not
   `/home`, which has a 50 GB quota:

   ```bash
   S=/rds/project/rds-XrHDlpCeVDg/users/jeremi/99-scratch/mamapper_release
   TMPDIR=$S/tmp python -m build --outdir $S/dist .
   python -m twine check $S/dist/*
   ```

   The sdist must contain `docs/`, `examples/` and `tests/` (see
   `MANIFEST.in`), because that is what the journal means by test data and
   documentation being available.

5. **Install the wheel into an empty environment and run it**, which is the
   only honest check of "installable without significant investment of time":

   ```bash
   python -m venv $S/venv
   $S/venv/bin/pip install $S/dist/ma_mapper-<version>-py3-none-any.whl
   $S/venv/bin/ma-mapper info
   $S/venv/bin/python examples/toy_end_to_end.py $S/clean_check
   ```

6. **Tag and push:**

   ```bash
   git tag -a v<version> -m "ma_mapper <version>"
   git push origin v<version>
   ```

7. **Upload to PyPI** (needs a PyPI API token):

   ```bash
   python -m twine upload $S/dist/*
   ```

   Then confirm `pip install ma-mapper` works from a clean environment.

8. **Zenodo** (needs the GitHub–Zenodo link enabled once, under Settings →
   Webhooks for the repository): enabling it means every GitHub *release* is
   archived automatically and gets a DOI. Create the release from the tag, then
   put the resulting DOI in `CITATION.cff` (`doi:`), the README and the paper's
   availability statement. `.zenodo.json` supplies the metadata.

9. **Bioconda** (needs a GitHub fork): follow `recipe/README.md` — the recipe
   is written, and only the version and the sdist checksum have to be filled in.
   Bioconda builds from PyPI, so step 7 comes first.

10. **Check the availability claims the paper will make**: the repository is
    public, the licence file is present, the DOI resolves, `pip install` works,
    and `conda install -c bioconda ma-mapper` works.

## What the journal asks for, and where it is

| requirement (checked against the author guidelines, 2026-10-11) | where it stands |
|---|---|
| Source at a stable URL at submission | public GitHub repository |
| Available two years after publication | GitHub plus the Zenodo archive |
| Open-source licence (encouraged) | MIT, `LICENSE` |
| Archived at Zenodo/Figshare/Software Heritage | step 8 |
| Installable without significant investment | PyPI and Bioconda, steps 7 and 9; verified by step 5 |
| Test data and documentation, including every optional parameter | `tests/`, `examples/`, `docs/package_docs.md`, `--help` on every subcommand |
| How to reproduce the paper's figures | one script per figure, under `paper/` — **not written yet, needs the figure** |
| Data availability statement | the case study must use public data only |
| AI-use declaration (ISCB policy) | the supplement, the cover letter and the Methods or Acknowledgements. Note the policy's limit: *"It is not acceptable to use LLMs or related technologies to draft papers"* — the manuscript text is written by the authors |

## Version numbering

`0.9.0` is the first release with tests, a command line and a build. `1.0.0`
should be the version the paper cites, cut once the name is settled, since the
package name cannot change afterwards without a second PyPI project.
