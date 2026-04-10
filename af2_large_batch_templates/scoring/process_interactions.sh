#!/usr/bin/env bash
set -euo pipefail

INPUT_BUCKET=$1
CUSTOM_THRESH=$2

GCP_INPUT_PATH="${INPUT_BUCKET}/predictions"
GCP_OUTPUT_PATH="${INPUT_BUCKET}/scored"
WORKDIR="/home/miles"
INPUT="${WORKDIR}/predictions"
OUTPUT="${WORKDIR}/scored"
## Permitting adjustable thresholds
if(${CUSTOM_THRESH} == 1){
    MINA=$3
    MAXA=$4
} else {
    MINA="1.5"
    MAXA="5.0"
}

## Reporting to log
echo "====== RB FF50-Multimer Predictions ======"
echo "Parameters:"
echo " - Date:        `date`"
echo " - Input Path:  ${GCP_INPUT_PATH}"
echo " - Output Path: ${GCP_OUTPUT_PATH}"
echo "==== Starting Run: `date` ===="

## generating file list and removing log files
gcloud storage ls ${GCP_INPUT_PATH} > temp_file_list
grep "/batch" temp_file_list > file_list
for temp_batch in $(cat file_list)
do
    batch_num=$(basename ${temp_batch})
    echo "== Processing ${batch_num}: `date`"

    ## Move files to VM
    gcloud storage cp $temp_batch"*" ${INPUT}/
    
    ## Scoring step 1 - pDockQ, iPAE, contact nums
    python ${WORKDIR}/scoring_interactions.py ${INPUT} \
        ${OUTPUT}/${batch_num}_rawinteractions.csv \
        ${OUTPUT}/${batch_num}_scores.csv \
        --min_distance ${MINA} \
        --max_distance ${MAXA}

    ## Scoring step 2 - pTM, ipTM
    grep -o "ptm.*" ./predictions/*.json | sed "s/_scores_/,/" | sed "s/_alphafold2.*json.ptm.. /,/" | sed "s/.*\\///" | sed "s/ .iptm...//" | sed "s/}//" >> ${OUTPUT}/${batch_num}_tmscores.csv

    ## Secondary structure
    python ${WORKDIR}/secondary_structure.py ${INPUT} \
        ${OUTPUT}/${batch_num}_ssfull.csv

    python ${WORKDIR}/merge_results.py ${OUTPUT}/${batch_num}_scores.csv \
        ${OUTPUT}/${batch_num}_tmscores.csv \
        ${OUTPUT}/${batch_num}_ssfull.csv \
        ${OUTPUT}/${batch_num}_final_summary.csv

    rm ${INPUT}/*
    gcloud storage cp ${OUTPUT}/ ${GCP_OUTPUT_PATH}/
    #rm ${OUTPUT}/*
done
