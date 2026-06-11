#!/usr/bin/env bash
set -euo pipefail

## To run use: ./process_interactions.sh "$ROOT_PATH" "$MINA" "$MAXA" "$PAEC" "$DISTC"
## Use the af2_ss or af2_scoring conda environment, they should both have the same dependencies.

INPUT_BUCKET=$1
## pDockQ thresholds
MINA=$2
MAXA=$3
## ipSAE thresholds
PAEC=$4
DISTC=$5

GCP_INPUT_PATH="${INPUT_BUCKET}/predictions"
GCP_OUTPUT_PATH="${INPUT_BUCKET}/scored"
WORKDIR="/home/resurrect/resurrectbio/01_projects/04_rb_automation/rb_automation/af2_large_batch_templates/scoring"
INPUT="${WORKDIR}/predictions"
OUTPUT="${WORKDIR}/scored"

## exporting PATH variable
#export PATH="/home/miles/miniconda3/envs/af2_scoring/bin:${PATH}"

## Reporting to log
echo "====== RB FF50-Multimer Predictions ======"
echo "Parameters:"
echo " - Date:        `date`"
echo " - Input Path:  ${GCP_INPUT_PATH}"
echo " - Output Path: ${GCP_OUTPUT_PATH}"
echo "==== Starting Run: `date` ===="

## Scoring step 1 - pDockQ, iPAE, contact nums
python ${WORKDIR}/scoring_interactions.py ${INPUT} \
    ${OUTPUT}/run_rawinteractions.csv \
    ${OUTPUT}/run_scores.csv \
    --min_distance ${MINA} \
    --max_distance ${MAXA} &&
    echo "Scoring step 1 complete"

## Scoring step 2 - pTM, ipTM
grep -o "ptm.*" ./predictions/*.json | sed "s/_scores_/,/" | sed "s/_alphafold2.*json.ptm.. /,/" | sed "s/.*\\///" | sed "s/ .iptm...//" | sed "s/}//" >> ${OUTPUT}/run_tmscores.csv
echo "Scoring step 2 complete"

## Scoring step 3 - ipSAE
## python ipsae.py <pae_json> <pdb> <pae_cutoff> <dist_cutoff>   
## Does this have to be a loop? - could use xargs, but this is fine for now.  
for pdb_file in ${INPUT}/*.pdb
do
    base_name=$(basename ${pdb_file})
    json_file=`echo ${pdb_file} | sed -e "s/_unrelaxed.*/_predicted_aligned_error_v1.json/"`
    python ${WORKDIR}/ipsae.py ${json_file} ${pdb_file} ${PAEC} ${DISTC}
done
grep -h "max" ${INPUT}/*txt >> ${OUTPUT}/run_ipsaefull.csv
sed -Ei 's/[[:space:]]+/,/g' ${OUTPUT}/run_ipsaefull.csv &&
echo "Scoring step 3 complete"

## STRETCH GOAL: Eventually add some kind of TMscore to estimate if the complex is stable or not. Not for version 1.

## Scoring step 4 - Secondary structure
python ${WORKDIR}/secondary_structure.py ${INPUT} \
    ${OUTPUT}/run_ssfull.csv &&
echo "Scoring step 4 complete"

python ${WORKDIR}/merge_results.py ${OUTPUT}/run_scores.csv \
    ${OUTPUT}/run_tmscores.csv \
    ${OUTPUT}/run_ssfull.csv \
    ${OUTPUT}/run_ipsaefull.csv \
    ${OUTPUT}/run_final_summary.csv &&
echo "All scoring steps complete"