# Bioconda recipe

Ready to submit; two placeholders have to be filled at release time.

1. Build the sdist and get its checksum:

   ```bash
   python -m build
   openssl sha256 dist/ma_mapper-<version>.tar.gz
   ```

2. Publish to PyPI first. Bioconda builds from the PyPI sdist, not from git.

3. Put the checksum in `meta.yaml` (`sha256:`), set `version:`, and leave
   `build: number: 0` for a new version.

4. Fork `bioconda/bioconda-recipes`, copy this file to
   `recipes/ma-mapper/meta.yaml`, and open a pull request. Their CI builds and
   tests it; the `test:` block above is what it runs.

`noarch: python` is correct here, because the package itself is pure Python. The
compiled dependencies (pysam, cyvcf2, pyBigWig) carry their own builds, and
`bedtools` is the binary `pybedtools` shells out to, which is the main reason
to offer a conda package at all.
