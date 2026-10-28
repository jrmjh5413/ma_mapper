#!/usr/bin/env bash
# Run the checks and record that this exact commit passed them.
#
# The same two the GitHub Actions workflow ran: ruff over the package and the
# tests, then pytest. The record is what `githooks/pre-push` looks for, so a
# commit cannot be pushed until its checks have been run here. It is written
# only when the working tree is clean, because otherwise "this commit passed"
# is not what was checked.
#
#   bash scripts/run_tests.sh              # run them now
#   sbatch scripts/run_tests.slurm         # the same, on a CSD3 compute node
#
# Records live in .git/test-records/, which is per clone and never committed.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

if [ -n "$(git status --porcelain -uno)" ]; then
    echo "The working tree has uncommitted changes."
    echo "Commit them first: a record says that a commit passed, not that a draft did."
    git status --short -uno
    exit 1
fi

sha=$(git rev-parse HEAD)
records=.git/test-records
mkdir -p "$records"

echo "commit   $sha"
echo "python   $(command -v python)"
echo "bedtools $(command -v bedtools || echo 'NOT ON PATH - the tests need it')"
echo

# ruff is a dev extra and is not in every environment. It is run when it is
# there, skipped loudly when it is not, and either way the record says which
# checks actually ran, so a record is never read as more than it is.
if command -v ruff > /dev/null; then
    ruff check ma_mapper tests
    checks="ruff,pytest"
else
    echo "ruff is not on PATH, so the lint was NOT run (pip install -e '.[dev]' adds it)."
    checks="pytest only - ruff not installed"
fi

python -m pytest -q

{
    echo "commit=$sha"
    echo "when=$(date -Is)"
    echo "host=$(hostname)"
    echo "python=$(python -V 2>&1)"
    echo "checks=$checks"
} > "$records/$sha"

echo
echo "PASSED - recorded in $records/$sha"
echo "This commit can now be pushed."
