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

FASTA_DIR="data/mmseqs_clusters/dedup"
CONFIG_FILE="config/query_configs.tsv"

echo "[SLURM-INFO] Job started on $(date)"

mkdir -p logs data/reps

source ~/scratch/miniconda/etc/profile.d/conda.sh
module load cuda
conda activate esm_env

echo "[SLURM-INFO] Active Python: $(which python)"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "[SLURM-INFO] Error: Config file $CONFIG_FILE not found."
    exit 1
fi

# Extract 5th column (label), skipping header (NR>1) and empty values
files=()
while IFS= read -r label; do
    target_file="${FASTA_DIR}/${label}_dedup.fasta"
    if [ -f "$target_file" ]; then
        files+=("$target_file")
    else
        echo "[SLURM-INFO] Warning: $target_file (label: $label) not found."
    fi
done < <(awk -F'\t' 'NR > 1 && $5 != "" {print $5}' "$CONFIG_FILE")

if [ "${#files[@]}" -eq 0 ]; then
    echo "[SLURM-INFO] Error: No matching fasta files found for labels in $CONFIG_FILE."
    exit 1
elif [ "$SLURM_ARRAY_TASK_ID" -ge "${#files[@]}" ]; then
    echo "[SLURM-INFO] Task ID $SLURM_ARRAY_TASK_ID exceeds available matched files (${#files[@]}). Exiting cleanly."
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