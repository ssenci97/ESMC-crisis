#!/bin/bash
#SBATCH --job-name=UPquery
#SBATCH --partition=THIN
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/uniprot.log
#SBATCH --array=1-15
#SBATCH --open-mode=append

OUTDIR="data/UniProtKB"
CONFIG="config/query_configs.tsv"

mkdir -p logs "${OUTDIR}"

ROW=$(tail -n +2 "${CONFIG}" | sed -n "${SLURM_ARRAY_TASK_ID}p")

PROTEOME_ID=$(echo "${ROW}" | cut -f1)
EXACT_SPECIES=$(echo "${ROW}" | cut -f4)
LABEL=$(echo "${ROW}" | cut -f5)

if [ -z "${PROTEOME_ID}" ]; then
  echo "[SLURM-INFO] Task ${SLURM_ARRAY_TASK_ID} exceeds number of entries in TSV. Exiting."
  exit 0
fi

echo "Job started on $(date)"
echo "[SLURM-INFO] Task ${SLURM_ARRAY_TASK_ID} targeting Proteome ID: ${PROTEOME_ID} (${EXACT_SPECIES})"

python src/download_seqs.py \
  --query "proteome:${PROTEOME_ID}" \
  --label "${LABEL}" \
  --fp_seqs "${OUTDIR}/${LABEL}_seqs.tsv" \
  --fp_fasta "${OUTDIR}/${LABEL}_seqs.fasta"

python src/download_features.py \
  --query "proteome:${PROTEOME_ID}" \
  --label "${LABEL}" \
  --fp_feats "${OUTDIR}/${LABEL}_feats.tsv"