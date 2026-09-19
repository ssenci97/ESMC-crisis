#!/bin/bash
#SBATCH --job-name=mmseqs_aggregate
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --output=logs/dedup_merge.log
#SBATCH --open-mode=append
set -euo pipefail

DEDUP_DIR="data/mmseqs_clusters/dedup/"
CONFIG_PATH="config/query_configs.tsv"

# 1. Merge Group FASTA files
declare -A GROUP_FILES
if [ -f "$CONFIG_PATH" ]; then
  while IFS=$'\t' read -r p_id ord gen sp label tax group alt; do
    if [[ "$label" != "label" && -n "$group" && -n "$label" ]]; then
      label_file="${DEDUP_DIR}/${label}_dedup.fasta"
      if [ -f "$label_file" ]; then
        GROUP_FILES["$group"]+="$label_file "
      fi
    fi
  done < "$CONFIG_PATH"

  for grp in "${!GROUP_FILES[@]}"; do
    grp_fasta="${DEDUP_DIR}/${grp}_dedup.fasta"
    cat ${GROUP_FILES[$grp]} > "$grp_fasta"
  done
fi

# 2. Merge Union FASTA file
union_fasta="${DEDUP_DIR}/union_dedup.fasta"
mapfile -t INDIV_DEDUP < <(find "${DEDUP_DIR}" -maxdepth 1 -type f -name "*_dedup.fasta" ! -name "union_dedup.fasta" | sort)

SPECIES_FILES=()
for file in "${INDIV_DEDUP[@]}"; do
  base=$(basename "$file" _dedup.fasta)
  if [[ ! " ${!GROUP_FILES[*]} " =~ " ${base} " ]]; then
    SPECIES_FILES+=("$file")
  fi
done

if [ ${#SPECIES_FILES[@]} -gt 0 ]; then
  cat "${SPECIES_FILES[@]}" > "$union_fasta"
fi

executiontime="${SECONDS}s"
maxmemory="$(grep VmPeak /proc/$$/status 2>/dev/null | awk '{print $2$3}' || echo 'N/A')"
echo -e "\n[SLURM-INFO] Task:AGGREGATE GroupsMerged:${!GROUP_FILES[*]} UnionTarget:${union_fasta} Time:${executiontime} Mem:${maxmemory} Date:$(date +'%Y-%m-%d %H:%M:%S')\n"