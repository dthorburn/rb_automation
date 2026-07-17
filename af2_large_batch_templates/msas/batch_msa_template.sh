#!/usr/bin/env bash
set -euo pipefail

## Taking the env variables
BATCH_ID=$1
INPUT_BUCKET=$2

GCP_INPUT_PATH="${INPUT_BUCKET}/fastas"
GCP_OUTPUT_PATH="${INPUT_BUCKET}/msas/batch${BATCH_ID}"
WORKDIR="/home/miles"
DB_PATH=${WORKDIR}/dbs/alphafold
FF_PATH=${WORKDIR}/dbs/florafold
INPUT="${WORKDIR}/input"
OUTPUT="${WORKDIR}/msas"
THREADS=26

## Setting environment
export PATH="/home/miles/localcolabfold/.pixi/envs/default/bin/:${PATH}"

## Reporting to log
echo "====== RB FF50-Multimer MSA Generation ======"
echo "Parameters:"
echo " - Date:        `date`"
echo " - Batch:       ${BATCH_ID}"
echo " - Input Path:  ${GCP_INPUT_PATH}"
echo " - Output Path: ${GCP_OUTPUT_PATH}"
echo "==== Starting Run: `date` ===="

## Importing data
echo "== Importing data: `date` =="
gsutil -m cp ${GCP_INPUT_PATH}/*0${BATCH_ID}.fasta ${INPUT}/

echo "== Splitting input: `date` =="
seqkit split2 -s 1 -O ${INPUT}/split ${INPUT}/*fasta
for f in ${INPUT}/split/*.fasta; do
    id=$(grep -m1 "^>" "$f" | sed 's/^>//; s/[ \t].*$//; s/[^A-Za-z0-9_]/_/g')
    mv -f "$f" "${INPUT}/split/${id}.fasta"
done

echo "== Starting Predictions: `date` =="
colabfold_search \
  --db1 ${DB_PATH}/uniref30_2302_db \
  --db3 ${FF_PATH}/florafolddb_reduced50 \
  --threads ${THREADS} \
  --db-load-mode 2 \
  --use-env 1 \
  ${INPUT}/split ${DB_PATH} ${OUTPUT} &&
    echo "== Finished predicitons: `date` =="

echo "== Exporting data to GCP: `date` =="
gsutil -m cp ${OUTPUT}/* ${GCP_OUTPUT_PATH}/
echo "==== Run Complete: `date` ===="
