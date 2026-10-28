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
REPO=/rds/project/rds-XrHDlpCeVDg/users/jeremi/02-ABCDom/ma_mapper
mkdir -p /rds/project/rds-XrHDlpCeVDg/users/jeremi/99-scratch/ma_mapper_tests
export TMPDIR=/local/$SLURM_JOB_ID; mkdir -p "$TMPDIR"
export MPLCONFIGDIR=$TMPDIR NUMBA_CACHE_DIR=$TMPDIR PYTHONDONTWRITEBYTECODE=1
export PATH=/home/jj515/miniconda3/envs/ma_mapper_tests/bin:$PATH
echo "env: $(which python)"
cd "$REPO"
/usr/bin/time -v bash scripts/run_tests.sh
rm -rf "$TMPDIR"
