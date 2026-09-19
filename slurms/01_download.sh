#!/bin/bash
# Slurm array job script to parse configuration rows, check for existing target 
# files to exit success on skip, download custom alt_dw links using label naming, 
# uncompress, unwrap FASTA sequences, extract IDs, and run standard pipelines.

#SBATCH --job-name=UPquery
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=03:00:00
#SBATCH --output=logs/uniprot.log
#SBATCH --array=1-16
#SBATCH --open-mode=append

###################################### initialization
OUTDIR="data/raw_seqs/"
CONFIG="config/query_configs.tsv"
START_TIME=$(date +%s)
mkdir -p logs "${OUTDIR}"
echo "Job started on $(date)"

###################################### parse-config
eval "$(python3 -c "
import pandas as pd
df = pd.read_csv('${CONFIG}', sep='\t', dtype=str)
idx = ${SLURM_ARRAY_TASK_ID} - 1
print('EXCEEDED=true' if idx >= len(df) else 'EXCEEDED=false')
row = df.iloc[idx] if idx < len(df) else None
if row is not None:
    for k in ['proteome_id', 'exact_species', 'label', 'alt_dw']:
        val = str(row.get(k, '')).strip()
        if val in ['nan', '\"\"', '']: val = ''
        print(f'{k.upper()}={val}')
")"

if [ "${EXCEEDED}" = "true" ]; then
    echo "[SLURM-INFO] Task ${SLURM_ARRAY_TASK_ID} exceeds number of entries in TSV. Exiting."
    exit 0
fi

if [ -z "${PROTEOME_ID}" ] && [ -z "${ALT_DW}" ]; then
    echo "[SLURM-INFO] Task ${SLURM_ARRAY_TASK_ID}: both proteome_id and alt_dw are missing for label '${LABEL}'. Skipping."
    exit 0
fi

###################################### download-alt-dw
if [ -n "${ALT_DW}" ]; then
    if [ -n "${LABEL}" ]; then
        FASTA_OUT="${OUTDIR}/${LABEL}_seqs.fasta"
        IDS_OUT="${OUTDIR}/${LABEL}_ids.txt"
    else
        FILENAME=$(basename "${ALT_DW}")
        FASTA_OUT="${OUTDIR}/${FILENAME%.gz}"
        IDS_OUT="${OUTDIR}/${FASTA_OUT%.*}_ids.txt"
    fi
    
    if [ -f "${FASTA_OUT}" ] && [ -f "${IDS_OUT}" ]; then
        echo "[SLURM-INFO] Task ${SLURM_ARRAY_TASK_ID} target files for label '${LABEL:-custom}' already exist. Skipping and exiting success."
        exit 0
    fi
    
    TMP_DOWNLOAD="${OUTDIR}/temp_$(basename "${ALT_DW}")"
    wget -q "${ALT_DW}" -O "${TMP_DOWNLOAD}"
    
    if [[ "${TMP_DOWNLOAD}" == *.gz ]]; then
        gunzip -f "${TMP_DOWNLOAD}"
        TMP_DOWNLOAD="${TMP_DOWNLOAD%.gz}"
    fi
    
    awk '/^>/ {if (n++) print ""; print; next} {printf "%s", $0} END {print ""}' "${TMP_DOWNLOAD}" > "${FASTA_OUT}"
    rm -f "${TMP_DOWNLOAD}"
    
    grep '^>' "${FASTA_OUT}" | sed 's/^>//' | awk '{print $1}' > "${IDS_OUT}"
    
    echo "[SLURM-INFO] Task ${SLURM_ARRAY_TASK_ID} successfully processed alt_dw link and extracted IDs for: ${LABEL:-custom}"
fi

###################################### run-pipelines
if [ -n "${PROTEOME_ID}" ]; then
    if [ -f "${OUTDIR}/${LABEL}_seqs.tsv" ] && [ -f "${OUTDIR}/${LABEL}_seqs.fasta" ] && [ -f "${OUTDIR}/${LABEL}_ids.txt" ] && [ -f "${OUTDIR}/${LABEL}_feats.tsv" ]; then
        echo "[SLURM-INFO] Task ${SLURM_ARRAY_TASK_ID} target files for label '${LABEL}' already exist. Skipping and exiting success."
        exit 0
    fi
    
    echo "[SLURM-INFO] Task ${SLURM_ARRAY_TASK_ID} targeting Label: ${LABEL} (Species: ${EXACT_SPECIES}, Proteome: ${PROTEOME_ID})"
    
    python src/download_seqs.py \
        --label "${LABEL}" \
        --fp_seqs "${OUTDIR}/${LABEL}_seqs.tsv" \
        --fp_fasta "${OUTDIR}/${LABEL}_seqs.fasta" \
        --fp_ids "${OUTDIR}/${LABEL}_ids.txt"
        
    python src/download_features.py \
        --label "${LABEL}" \
        --fp_feats "${OUTDIR}/${LABEL}_feats.tsv"
fi

###################################### finalize
END_TIME=$(date +%s)
EXEC_TIME=$((END_TIME - START_TIME))
MAX_MEM=$(free -m | awk '/Mem:/ {print $3"MB"}')
echo "[SLURM-INFO] Task ${SLURM_ARRAY_TASK_ID} finished in ${EXEC_TIME}s, Memory used: ${MAX_MEM} at $(date)"