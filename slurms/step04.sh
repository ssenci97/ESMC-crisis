#!/bin/bash
#SBATCH --job-name=npzmerge
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/merge_npz.log

###################################### PARAMETERS
CONFIG_PATH="config/query_configs.tsv"
REPS_DIR="data/reps"
MERGE_SCRIPT="src/npz_merge.py"

###################################### ENVIRONMENT
start_time=$(date +%s)
source /orfeo/cephfs/scratch/area/ssenci/venvs/ml_transformers/bin/activate

###################################### MAIN
awk -F'\t' -v dir="$REPS_DIR" 'NR>1 && $8!="" && $6!="" {groups[$8] = groups[$8] ? groups[$8] "," dir "/" $6 "_reps.npz" : dir "/" $6 "_reps.npz"} END {for (g in groups) print g, groups[g]}' "$CONFIG_PATH" | while read -r group npzs; do
    out_npz="${REPS_DIR}/${group}_reps.npz"
    if [ -f "$out_npz" ]; then
        echo "[SLURM-INFO] Skipping group '${group}': target '${out_npz}' already exists"
    else
        echo "[SLURM-INFO] Merging group '${group}' from [${npzs}] into '${out_npz}'"
        python3 "$MERGE_SCRIPT" --npzs "$npzs" --out_npz "$out_npz"
    fi
done

###################################### LOGGING
end_time=$(date +%s)
executiontime="$((end_time - start_time))s"
maxmemory="$(ps -o rss= -p $$ | tr -d ' ')KB"
time_date=$(date "+%Y-%m-%d %H:%M:%S")
echo "[SLURM-INFO] ${executiontime} ${maxmemory} ${time_date}"