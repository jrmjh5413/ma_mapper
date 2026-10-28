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

ruff check ma_mapper tests
python -m pytest -q

{
    echo "commit=$sha"
    echo "when=$(date -Is)"
    echo "host=$(hostname)"
    echo "python=$(python -V 2>&1)"
} > "$records/$sha"

echo
echo "PASSED - recorded in $records/$sha"
echo "This commit can now be pushed."
