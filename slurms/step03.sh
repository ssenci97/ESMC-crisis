#!/bin/bash
#SBATCH --partition=DGX
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --gres=gpu:A100:1
#SBATCH --time=04:00:00
#SBATCH --job-name=ESMC300m
#SBATCH --output=logs/esm300m.log
#SBATCH --array=0-15
#SBATCH --open-mode=append

FASTA_DIR="data/UniProtKB"

echo "[SLURM-INFO] Job started on $(date)"

mkdir -p logs data/reps

source ~/scratch/miniconda/etc/profile.d/conda.sh
module load cuda
conda activate esm_env

echo "[SLURM-INFO] Active Python: $(which python)"

shopt -s nullglob
files=("$FASTA_DIR"/*_dedup.fasta)
shopt -u nullglob

if [ "${#files[@]}" -eq 0 ]; then
    echo "[SLURM-INFO] Error: No matching fasta files found in $FASTA_DIR/"
    exit 1
elif [ "$SLURM_ARRAY_TASK_ID" -ge "${#files[@]}" ]; then
    echo "[SLURM-INFO] Task ID $SLURM_ARRAY_TASK_ID exceeds available files (${#files[@]}). Exiting cleanly."
    exit 0
fi

input_file="${files[$SLURM_ARRAY_TASK_ID]}"
filename=$(basename "$input_file")
base_name="${filename%_dedup.fasta}"
base_name="${base_name%_seqs}"
output_file="data/reps/${base_name}_reps.npz"

echo "[SLURM-INFO] Target input: $input_file"
echo "[SLURM-INFO] Target output: $output_file"

if [ -f "$output_file" ]; then
    echo "[SLURM-INFO] SUCCESS: $output_file already exists. Skipping computation."
    exit 0
fi

time python src/extract_reps.py --fasta "$input_file" --out "$output_file"

echo "[SLURM-INFO] Done $(date)"