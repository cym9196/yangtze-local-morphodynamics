#!/bin/bash
#SBATCH --job-name=yangtze-iron-63s
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm-%j.out

# Generic Slurm template. Adjust the partition/account/module lines to the
# target cluster; no cloud or paid resource is assumed by this repository.
set -euo pipefail

CASE_DIR="${SLURM_SUBMIT_DIR:?Submit this script from the case directory}"
cd "$CASE_DIR"

# Sites that provide a module may replace this with: module load openfoam/2412
OPENFOAM_BASHRC="${OPENFOAM_BASHRC:-/usr/lib/openfoam/openfoam2412/etc/bashrc}"
source "$OPENFOAM_BASHRC"

# Compile once on the cluster architecture before launching the case.
(
    cd solver/immersedSedExnerFoam
    ./Allwmake
)

# Allrun63s recognizes Slurm and uses srun. It also updates the decomposition
# count to match the allocation and reconstructs only the final time.
NPROCS="$SLURM_NTASKS" ./Allrun63s
