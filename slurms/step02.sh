#!/bin/bash
#SBATCH --job-name=mmseqs
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=9
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/mmseqs2_clust.log
#SBATCH --array=1-15
#SBATCH --open-mode=append
set -euo pipefail
MMSEQS='/orfeo/cephfs/scratch/area/ssenci/pkgs/mmseqs_new/mmseqs/bin/mmseqs'
THREADS=${SLURM_CPUS_PER_TASK:-16}
IN_DIR="data/UniProtKB"
OUT_BASE_DIR="data/mmseqs_clusters"
declare -A THRESHOLDS_LABELS
THRESHOLDS=("0.90" "0.70" "0.50" "0.40" "0.35" "0.30" "0.25" "0.20" "0.15")
THRESHOLDS_LABELS=([0.90]="strains" [0.70]="subfamily" [0.50]="family" [0.40]="superfamily" [0.35]="broad" [0.30]="verybroad" [0.25]="superbroad" [0.20]="ultrabroad" [0.15]="extreme")
###################################### <SETUP>
mkdir -p logs "${IN_DIR}" "${OUT_BASE_DIR}"
mapfile -t FASTA_FILES < <(find "${IN_DIR}" -maxdepth 1 -type f -name "*.fasta" ! -name "*_dedup.fasta" | sort)
ARRAY_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
fp="${FASTA_FILES[$ARRAY_INDEX]:-}"
if [ -z "${fp}" ]; then
  echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Task ${SLURM_ARRAY_TASK_ID} exceeds number of FASTA files found. Exiting."
  exit 0
fi
base_name=$(basename "$fp" .fasta)
base_name=${base_name%_seqs}
base_name=${base_name%_dedup}
dedup_fasta="${IN_DIR}/${base_name}_dedup.fasta"
out_dir="${OUT_BASE_DIR}/${base_name}"
tmp_dir="${out_dir}/tmp"
mkdir -p "$out_dir" "$tmp_dir"
echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Starting Task ${SLURM_ARRAY_TASK_ID} on host $(hostname) for ${base_name}"
###################################### <DEDUPLICATION>
if [ ! -f "$dedup_fasta" ]; then
    "$MMSEQS" easy-linclust "$fp" "${out_dir}/${base_name}_dedup" "$tmp_dir" --min-seq-id 1.0 -c 0.8 --cov-mode 5 --threads "$THREADS" || exit 1
    mv "${out_dir}/${base_name}_dedup_rep_seq.fasta" "$dedup_fasta"
fi
if [ ! -f "${out_dir}/dedupDB.index" ]; then
    "$MMSEQS" createdb "$dedup_fasta" "${out_dir}/dedupDB" || exit 1
fi
total_seqs=$(grep -c "^>" "$dedup_fasta" 2>/dev/null || echo "0")
original_seqs=$(grep -c "^>" "$fp" 2>/dev/null || echo "0")
duplicates_removed=$((original_seqs - total_seqs))
reduction_pct=$(echo "scale=2; (($original_seqs - $total_seqs) * 100) / $original_seqs" | bc -l 2>/dev/null || echo "0")
metrics_file="${out_dir}/${base_name}_cluster_metrics.tsv"
if [ ! -f "$metrics_file" ]; then
    echo -e "threshold\tclusters\tsingleton_clusters\tmax_cluster_size\tavg_cluster_size\tcoverage" > "$metrics_file"
fi
###################################### <CLUSTERING>
for seq_id in "${THRESHOLDS[@]}"; do
    threshold_clean="${seq_id//./_}"
    label="${THRESHOLDS_LABELS[${seq_id}]}"
    cluster_db="${out_dir}/clusterDB_id${threshold_clean}_${label}"
    tsv_out="${out_dir}/${base_name}_clusters_id${threshold_clean}_${label}.tsv"
    fasta_out="${out_dir}/${base_name}_representatives_id${threshold_clean}_${label}.fasta"
    if [ -f "$tsv_out" ] && [ -f "$fasta_out" ]; then
        continue
    fi
    cov=0.5
    "$MMSEQS" cluster "${out_dir}/dedupDB" "$cluster_db" "$tmp_dir" --min-seq-id "$seq_id" -c "$cov" --cov-mode 0 --threads "$THREADS" || continue
    "$MMSEQS" createtsv "${out_dir}/dedupDB" "${out_dir}/dedupDB" "$cluster_db" "$tsv_out" || continue
    "$MMSEQS" result2repseq "${out_dir}/dedupDB" "$cluster_db" "${cluster_db}_rep" || continue
    "$MMSEQS" result2flat "${out_dir}/dedupDB" "${out_dir}/dedupDB" "${cluster_db}_rep" "$fasta_out" --use-fasta-header 1 || continue
    num_clusters=$(awk -F'\t' '{print $1}' "$tsv_out" | sort -u | wc -l)
    singleton_clusters=$(awk -F'\t' '{count[$1]++} END {singletons=0; for (c in count) if (count[c]==1) singletons++; print singletons+0}' "$tsv_out")
    max_size=$(awk -F'\t' '{count[$1]++} END {max=0; for (c in count) if (count[c]>max) max=count[c]; print max}' "$tsv_out")
    avg_size=$(echo "scale=2; $total_seqs / $num_clusters" | bc -l)
    if ! grep -q "^${seq_id}\t" "$metrics_file"; then
        echo -e "${seq_id}\t${num_clusters}\t${singleton_clusters}\t${max_size}\t${avg_size}\t${cov}" >> "$metrics_file"
    fi
    rm -f "${cluster_db}"*
done
###################################### <CLEANUP>
summary_file="${out_dir}/clustering_summary.txt"
{
    echo -e "=== Clustering Summary for ${base_name} ===\nTimestamp: $(date +'%Y-%m-%d %H:%M:%S')\n\n--- Deduplication Stats ---"
    echo -e "Original sequences: ${original_seqs}\nDuplicates & fragments consolidated: ${duplicates_removed}"
    echo -e "Representative sequences retained: ${total_seqs}\nReduction: ${reduction_pct}%\n\n--- Clustering Metrics ---"
    [ -f "$metrics_file" ] && cat "$metrics_file" || echo "No metrics file found."
} > "$summary_file"
rm -rf "$tmp_dir"
rm -f "${out_dir}/dedupDB"*
executiontime="${SECONDS}s"
maxmemory="$(grep VmPeak /proc/$$/status 2>/dev/null | awk '{print $2$3}' || echo 'N/A')"
echo "[PYTHON-INFO] ${executiontime} ${maxmemory} $(date +'%Y-%m-%d %H:%M:%S')"