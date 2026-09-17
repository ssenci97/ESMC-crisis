#!/bin/bash
#SBATCH --job-name=mafft_global
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/merge.log

source /orfeo/cephfs/scratch/area/ssenci/venvs/ml_transformers/bin/activate

python3 src/npz_merge.py --npzs data/reps/aegypti_reps.npz,data/reps/albopictus_reps.npz --out_npz data/reps/aedes_reps.npz

python3 src/npz_merge.py --npzs data/reps/viridula_reps.npz,data/reps/halys_reps.npz --out_npz data/reps/stinkbug_reps.npz

python3 src/npz_merge.py --npzs data/reps/melanogaster_reps.npz,data/reps/suzukii_reps.npz --out_npz data/reps/drosophila_reps.npz
