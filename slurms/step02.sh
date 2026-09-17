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

mkdir -p logs "${IN_DIR}" "${OUT_BASE_DIR}"

# Download initial data if not present
TARGET_FASTA="${IN_DIR}/halys_seqs.fasta"

if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
    if [ ! -f "${TARGET_FASTA}" ]; then
        echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Task 1 downloading and processing sequence data..."
        curl -sSLf --retry 3 --retry-delay 5 \
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/696/795/GCF_000696795.3_Hhal_1.1/GCF_000696795.3_Hhal_1.1_protein.faa.gz" \
            -o "${IN_DIR}/halys_seqs.fasta.gz"
        gzip -d -c "${IN_DIR}/halys_seqs.fasta.gz" | \
        awk '/^>/ { if (NR > 1) printf("\n"); print $0; next } { printf("%s", $0) } END { printf("\n") }' \
        > "${TARGET_FASTA}.tmp"
        # Atomic rename ensures waiting jobs don't read a partially written file
        mv "${TARGET_FASTA}.tmp" "${TARGET_FASTA}"
        rm -f "${IN_DIR}/halys_seqs.fasta.gz"
    fi
else
    # Non-task-1 jobs wait until Task 1 finishes writing the target FASTA
    #echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Task ${SLURM_ARRAY_TASK_ID} waiting for Task 1 to complete download..."
    while [ ! -f "${TARGET_FASTA}" ]; do
        sleep 3
        echo "..."
    done
    #echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] ${TARGET_FASTA} detected. Proceeding with Task ${SLURM_ARRAY_TASK_ID}."
fi

# Build list of target FASTA files (excluding deduplicated files)
mapfile -t FASTA_FILES < <(find "${IN_DIR}" -maxdepth 1 -type f -name "*.fasta" ! -name "*_dedup.fasta" | sort)

# Select file corresponding to current array task ID (1-based index)
ARRAY_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
fp="${FASTA_FILES[$ARRAY_INDEX]:-}"

# Exit gracefully if task ID exceeds total number of files
if [ -z "${fp}" ]; then
  echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Task ${SLURM_ARRAY_TASK_ID} exceeds number of FASTA files found (${#FASTA_FILES[@]}). Exiting."
  exit 0
fi

declare -A THRESHOLDS_LABELS
THRESHOLDS=(
"0.90"
"0.70"
"0.50"
"0.40"
"0.35"
"0.30"
"0.25"
"0.20"
"0.15"
)

THRESHOLDS_LABELS[0.90]="strains"
THRESHOLDS_LABELS[0.70]="subfamily"
THRESHOLDS_LABELS[0.50]="family"
THRESHOLDS_LABELS[0.40]="superfamily"
THRESHOLDS_LABELS[0.35]="broad"
THRESHOLDS_LABELS[0.30]="verybroad"
THRESHOLDS_LABELS[0.25]="superbroad"
THRESHOLDS_LABELS[0.20]="ultrabroad"
THRESHOLDS_LABELS[0.15]="extreme"

echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Starting Task ${SLURM_ARRAY_TASK_ID} on host $(hostname)"
echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] MMseqs2 executable: ${MMSEQS}"
echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] CPU threads allocated: ${THREADS}"

base_name=$(basename "$fp" .fasta)
base_name=${base_name%_seqs}
base_name=${base_name%_dedup}

echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] =========================================="
echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Processing file ${SLURM_ARRAY_TASK_ID}: ${base_name}"

dedup_fasta="${IN_DIR}/${base_name}_dedup.fasta"
out_dir="${OUT_BASE_DIR}/${base_name}"
tmp_dir="${out_dir}/tmp"

mkdir -p "$out_dir" "$tmp_dir"

if [ -f "$dedup_fasta" ]; then
    echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Deduplicated file already exists: ${dedup_fasta}. Skipping deduplication step."
else
    echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Running duplicate + fragment consolidation..."
    
    "$MMSEQS" easy-linclust "$fp" "${out_dir}/${base_name}_dedup" "$tmp_dir" \
    --min-seq-id 1.0 \
    -c 0.8 \
    --cov-mode 5 \
    --threads "$THREADS" || {
        echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] ERROR: Deduplication failed for ${base_name}"
        exit 1
    }

    echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Moving deduplicated sequences to ${dedup_fasta}"
    mv "${out_dir}/${base_name}_dedup_rep_seq.fasta" "$dedup_fasta"
fi

if [ -f "${out_dir}/dedupDB.index" ]; then
    echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Database already exists, skipping createdb."
else
    echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Building MMseqs2 database..."
    "$MMSEQS" createdb "$dedup_fasta" "${out_dir}/dedupDB" || {
        echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] ERROR: Database creation failed for ${base_name}"
        exit 1
    }
fi

total_seqs=$(grep -c "^>" "$dedup_fasta" 2>/dev/null || echo "0")
original_seqs=$(grep -c "^>" "$fp" 2>/dev/null || echo "0")
duplicates_removed=$((original_seqs - total_seqs))
reduction_pct=$(echo "scale=2; (($original_seqs - $total_seqs) * 100) / $original_seqs" | bc -l 2>/dev/null || echo "0")

echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Original sequences: ${original_seqs}"
echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] After dedup: ${total_seqs} (removed: ${duplicates_removed}, reduction: ${reduction_pct}%)"

metrics_file="${out_dir}/${base_name}_cluster_metrics.tsv"
if [ ! -f "$metrics_file" ]; then
    echo -e "threshold\tclusters\tsingleton_clusters\tmax_cluster_size\tavg_cluster_size\tcoverage" > "$metrics_file"
fi

for seq_id in "${THRESHOLDS[@]}"; do
    threshold_clean="${seq_id//./_}"
    label="${THRESHOLDS_LABELS[${seq_id}]}"
    cluster_db="${out_dir}/clusterDB_id${threshold_clean}_${label}"
    tsv_out="${out_dir}/${base_name}_clusters_id${threshold_clean}_${label}.tsv"
    fasta_out="${out_dir}/${base_name}_representatives_id${threshold_clean}_${label}.fasta"

    if [ -f "$tsv_out" ] && [ -f "$fasta_out" ]; then
        echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Clustering at ${seq_id} (${label}) already completed. Extracting metrics..."
        num_clusters=$(awk -F'\t' '{print $1}' "$tsv_out" | sort -u | wc -l)
        singleton_clusters=$(awk -F'\t' '{count[$1]++} END {for (c in count) if (count[c]==1) singletons++; print singletons+0}' "$tsv_out")
        max_size=$(awk -F'\t' '{count[$1]++} END {max=0; for (c in count) if (count[c]>max) max=count[c]; print max}' "$tsv_out")
        avg_size=$(echo "scale=2; $total_seqs / $num_clusters" | bc -l)
        cov="0.5"

        if ! grep -q "^${seq_id}\t" "$metrics_file"; then
            echo -e "${seq_id}\t${num_clusters}\t${singleton_clusters}\t${max_size}\t${avg_size}\t${cov}" >> "$metrics_file"
        fi

        echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] → ${label} (${seq_id}): ${num_clusters} clusters | singletons: ${singleton_clusters} | max size: ${max_size} | avg: ${avg_size}"
        continue
    fi

    echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Clustering at ${seq_id} identity (${label})..."
    cov=0.5
    "$MMSEQS" cluster "${out_dir}/dedupDB" "$cluster_db" "$tmp_dir" \
    --min-seq-id "$seq_id" \
    -c "$cov" \
    --cov-mode 0 \
    --threads "$THREADS" || {
        echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] ERROR: Clustering at ${seq_id} failed"
        continue
    }

    echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Exporting TSV..."
    "$MMSEQS" createtsv "${out_dir}/dedupDB" "${out_dir}/dedupDB" "$cluster_db" "$tsv_out" || {
        echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] ERROR: TSV export at ${seq_id} failed"
        continue
    }

    echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Exporting representative FASTA..."
    "$MMSEQS" result2repseq "${out_dir}/dedupDB" "$cluster_db" "${cluster_db}_rep" || {
        echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] ERROR: result2repseq at ${seq_id} failed"
        continue
    }
    "$MMSEQS" result2flat "${out_dir}/dedupDB" "${out_dir}/dedupDB" "${cluster_db}_rep" "$fasta_out" --use-fasta-header 1 || {
        echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] ERROR: FASTA export at ${seq_id} failed"
        continue
    }

    num_clusters=$(awk -F'\t' '{print $1}' "$tsv_out" | sort -u | wc -l)
    singleton_clusters=$(awk -F'\t' '{count[$1]++} END {singletons=0; for (c in count) if (count[c]==1) singletons++; print singletons}' "$tsv_out")
    max_size=$(awk -F'\t' '{count[$1]++} END {max=0; for (c in count) if (count[c]>max) max=count[c]; print max}' "$tsv_out")
    avg_size=$(echo "scale=2; $total_seqs / $num_clusters" | bc -l)

    if ! grep -q "^${seq_id}\t" "$metrics_file"; then
        echo -e "${seq_id}\t${num_clusters}\t${singleton_clusters}\t${max_size}\t${avg_size}\t${cov}" >> "$metrics_file"
    fi

    echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] → ${label} (${seq_id}): ${num_clusters} clusters | singletons: ${singleton_clusters} | max size: ${max_size} | avg: ${avg_size}"

    rm -f "${cluster_db}"*
done

summary_file="${out_dir}/clustering_summary.txt"
{
    echo "=== Clustering Summary for ${base_name} ==="
    echo "Timestamp: $(date +'%Y-%m-%d %H:%M:%S')"
    echo ""
    echo "--- Deduplication Stats ---"
    echo "Original sequences: ${original_seqs}"
    echo "Duplicates & fragments consolidated: ${duplicates_removed}"
    echo "Representative sequences retained: ${total_seqs}"
    echo "Reduction: ${reduction_pct}%"
    echo ""
    echo "--- Clustering Metrics ---"
    if [ -f "$metrics_file" ]; then
        cat "$metrics_file"
    else
        echo "No metrics file found."
    fi
} > "$summary_file"

echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Summary saved to: ${summary_file}"
echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Cleaning up temporary directory and databases..."
rm -rf "$tmp_dir"
rm -f "${out_dir}/dedupDB"*
echo "[SLURM-INFO] [$(date +'%Y-%m-%d %H:%M:%S')] Task ${SLURM_ARRAY_TASK_ID} successfully completed processing: ${base_name}"