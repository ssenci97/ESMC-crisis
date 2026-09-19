#!/bin/bash
#SBATCH --job-name=scanpy
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/scanpy.log
#SBATCH --open-mode=append

source /orfeo/scratch/area/ssenci/venvs/ml_transformers/bin/activate
python3 src/scanpy_clusters.py --npzs data/reps/suzukii_reps.npz,data/reps/melanogaster_reps.npz
python3 src/scanpy_clusters.py --npzs data/reps/halys_reps.npz,data/reps/viridula_reps.npz
python3 src/scanpy_clusters.py --npzs data/reps/aegypti_reps.npz,data/reps/albopictus_reps.npz