#!/bin/bash
#SBATCH --job-name=dim32GB
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --array=0-10
#SBATCH --output=logs/dimred.log
#SBATCH --open-mode=append

mkdir -p slurms/logs data/dim_red
set -e

echo "[SLURM-INFO] .................................................."
echo "[SLURM-INFO] Job Name:    ${SLURM_JOB_NAME}"
echo "[SLURM-INFO] Master ID:   ${SLURM_ARRAY_JOB_ID}"
echo "[SLURM-INFO] Task ID:     ${SLURM_ARRAY_TASK_ID}"
echo "[SLURM-INFO] Node:        $(hostname)"
echo "[SLURM-INFO] Started at:  $(date)"
echo "[SLURM-INFO] .................................................."

source ~/scratch/miniconda/etc/profile.d/conda.sh
conda activate esm_env

echo "[SLURM-INFO] Active Python: $(which python3)"

# Dynamically find matching files and filter for size <= 80MB (81920 KiB)
files=()
for f in data/reps/*_reps.npz; do
    [ -e "$f" ] || continue
    # 80 MB in binary (MiB) is 80 * 1024 = 81920 KiB
    size_kb=$(du -k "$f" | cut -f1)
    if [ "$size_kb" -le 81920 ]; then
        files+=("$f")
        echo "[SLURM-INFO] Included (<= 80MB): $f (${size_kb} KiB)"
    else
        echo " ... \n"
    fi
done

total_files=${#files[@]}

if [ "$total_files" -eq 0 ]; then
    echo "[SLURM-INFO] ERROR: No files found matching criteria (> 80MB)! Exiting."
    exit 1
fi

if [ "$SLURM_ARRAY_TASK_ID" -ge "$total_files" ]; then
    echo "[SLURM-INFO] SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) exceeds total matching files ($total_files). Exiting with success."
    exit 0
fi

input_file="${files[$SLURM_ARRAY_TASK_ID]}"

# Extract the base name (e.g., 'aegypti' from 'data/reps/aegypti_reps.npz')
filename=$(basename "$input_file")
current_base="${filename%_reps.npz}"

tsne_out="data/dim_red/${current_base}_tsne.npz"
umap_out="data/dim_red/${current_base}_umap.npz"

echo "[SLURM-INFO] Processing input: $input_file"
echo "[SLURM-INFO] Target base: $current_base"

# Skip if targets already exist
if [ -f "$tsne_out" ] && [ -f "$umap_out" ]; then
    echo "[SLURM-INFO] SUCCESS: Both $tsne_out and $umap_out already exist. Skipping run."
    exit 0
fi

if [ -f "$tsne_out" ]; then
    echo "[SLURM-INFO] SUCCESS: $tsne_out already exists. Skipping t-SNE."
else
    echo "[SLURM-INFO] Starting t-SNE for $current_base..."
    python3 src/npz_tsne.py --npz "$input_file" --out "$tsne_out"
fi

if [ -f "$umap_out" ]; then
    echo "[SLURM-INFO] SUCCESS: $umap_out already exists. Skipping UMAP."
else
    echo "[SLURM-INFO] Starting UMAP for $current_base..."
    python3 src/npz_umap.py --npz "$input_file" --out "$umap_out"
fi

echo "[SLURM-INFO] .................................................."
echo "[SLURM-INFO] Finished dimensionality reduction for $current_base on $(date)"
echo "[SLURM-INFO] .................................................."