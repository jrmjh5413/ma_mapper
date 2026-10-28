#!/bin/bash
#SBATCH -J ma_mapper_tests -A IMBEAULT-SL3-CPU -p icelake
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --time=00:30:00 --mem=12G
#SBATCH -o /rds/project/rds-XrHDlpCeVDg/users/jeremi/99-scratch/ma_mapper_tests/tests_%j.out
#SBATCH -e /rds/project/rds-XrHDlpCeVDg/users/jeremi/99-scratch/ma_mapper_tests/tests_%j.err
# The checks on a compute node, never on a login node.
#
#   sbatch scripts/run_tests.slurm
#
# On success it writes the record that githooks/pre-push requires. Logs go to
# /rds/project/rds-XrHDlpCeVDg/users/jeremi/99-scratch/ma_mapper_tests/.
# Environment: /home/jj515/miniconda3/envs/ma_mapper_tests  (migration debt: not yet under users/jeremi/envs/)
set -euo pipefail
# The repository root is found from the directory sbatch was run in, never
# hardcoded: the tree is reorganised from time to time and a hardcoded path
# would quietly leave the gate testing somewhere else. Submit from anywhere
# inside the checkout.
REPO=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && git rev-parse --show-toplevel 2>/dev/null) || {
    echo "Submit this from inside the checkout: cd to the repository, then sbatch scripts/run_tests.slurm" >&2
    exit 1
}
mkdir -p /rds/project/rds-XrHDlpCeVDg/users/jeremi/99-scratch/ma_mapper_tests
export TMPDIR=/local/$SLURM_JOB_ID; mkdir -p "$TMPDIR"
export MPLCONFIGDIR=$TMPDIR NUMBA_CACHE_DIR=$TMPDIR PYTHONDONTWRITEBYTECODE=1
export PATH=/home/jj515/miniconda3/envs/ma_mapper_tests/bin:$PATH
echo "env: $(which python)"
cd "$REPO"
/usr/bin/time -v bash scripts/run_tests.sh
rm -rf "$TMPDIR"
