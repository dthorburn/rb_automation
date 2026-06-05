#!/usr/bin/env bash
set -euo pipefail

INPUT_BUCKET=$1
## pDockQ thresholds
MINA=$2
MAXA=$3
## ipSAE thresholds
PAEC=$4
DISTC=$5

GCP_INPUT_PATH="${INPUT_BUCKET}/predictions"
GCP_OUTPUT_PATH="${INPUT_BUCKET}/scored"
WORKDIR="/home/miles"
INPUT="${WORKDIR}/predictions"
OUTPUT="${WORKDIR}/scored"

## exporting PATH variable
export PATH="/home/miles/miniconda3/envs/af2_scoring/bin:${PATH}"

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

    ## Scoring step 3 - ipSAE
    ## python ipsae.py <pae_json> <pdb> <pae_cutoff> <dist_cutoff>   
    ## Does this have to be a loop? - could use xargs, but this is fine for now.  
    for pdb_file in ${INPUT}/*.pdb
    do
        base_name=$(basename ${pdb_file})
        json_file=`echo ${pdb_file} | sed -e "s/_unrelaxed.*/_predicted_aligned_error_v1.json/"`
        python ${WORKDIR}/ipsae.py ${json_file} ${pdb_file} ${PAEC} ${DISTC}
    done
    #echo "" > ${OUTPUT}/${batch_num}_ipsaefull.tsv
    grep -h "max" ${INPUT}/*txt >> ${OUTPUT}/${batch_num}_ipsaefull.csv
    sed -Ei 's/[[:space:]]+/,/g' ${OUTPUT}/${batch_num}_ipsaefull.csv

    ## STRETCH GOAL: Eventually add some kind of TMscore to estimate if the complex is stable or not. Not for version 1.

    ## Scoring step 4 - Secondary structure
    python ${WORKDIR}/secondary_structure.py ${INPUT} \
        ${OUTPUT}/${batch_num}_ssfull.csv

    python ${WORKDIR}/merge_results.py ${OUTPUT}/${batch_num}_scores.csv \
        ${OUTPUT}/${batch_num}_tmscores.csv \
        ${OUTPUT}/${batch_num}_ssfull.csv \
        ${OUTPUT}/${batch_num}_ipsaefull.csv \
        ${OUTPUT}/${batch_num}_final_summary.csv

    rm ${INPUT}/*
    gcloud storage cp ${OUTPUT}/${batch_num}_final_summary.csv ${GCP_OUTPUT_PATH}/
    mv ${OUTPUT}/${batch_num}_final_summary.csv ${WORKDIR}/completed/
    gcloud storage cp ${OUTPUT}/${batch_num}_* ${GCP_OUTPUT_PATH}/individual_scores/
    mv ${OUTPUT}/${batch_num}_* ${WORKDIR}/completed/
    #rm ${OUTPUT}/*
done
