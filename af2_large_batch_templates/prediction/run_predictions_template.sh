#!/usr/bin/env bash
set -euo pipefail

## Taking the env variables
BATCH_ID=$1
MODELS=$2
RECYCLES=$3
INPUT_BUCKET=$4

GCP_INPUT_PATH="${INPUT_BUCKET}/msas/batch${BATCH_ID}"
GCP_OUTPUT_PATH="${INPUT_BUCKET}/predictions/batch${BATCH_ID}"
WORKDIR="/home/miles"
TEMP_INPUT="${WORKDIR}/temp_msas"
INPUT="${WORKDIR}/msas"
OUTPUT="${WORKDIR}/predictions"

## Setting environment
export PATH="/home/miles/localcolabfold/.pixi/envs/default/bin/:${PATH}"
export LD_LIBRARY_PATH="/home/miles/localcolabfold/.pixi/envs/default/lib:/usr/local/cuda/lib64:/usr/local/nccl2/lib:/usr/local/cuda/extras/CUPTI/lib64"

## Reporting to log
echo "====== RB FF50-Multimer Predictions ======"
echo "Parameters:"
echo " - Date:       `date`"
echo " - Batch:      ${BATCH_ID}"
echo " - Models:     ${MODELS}"
echo " - Recycles:   ${RECYCLES}"
echo " - Input Path: ${GCP_INPUT_PATH}"
echo " - Output Path: ${GCP_OUTPUT_PATH}"
echo "==== Starting Run: `date` ===="

## Importing data
echo "== Importing data: `date` =="
mkdir -p ${TEMP_INPUT}
gsutil -m cp ${GCP_INPUT_PATH}/*.a3m ${TEMP_INPUT}/

## Cancelling script if no files were uploaded for any reason. 
if ! ls "${TEMP_INPUT}"/*.a3m 1> /dev/null 2>&1; then
    echo "==== No input files detected in ${TEMP_INPUT} -- exiting cleanly ===="
    exit 0
fi

## Chunking concept
cd ${WORKDIR}
ls -1 "${TEMP_INPUT}"/*.a3m | split -l 100 -d - batch_
for batch_num in ${WORKDIR}/batch_*
do
    echo $batch_num
    echo "== Starting ${batch_num} Predictions: `date` =="
    xargs -a ${batch_num} -I {} mv "{}" ${INPUT}/
    colabfold_batch \
        --disable-unified-memory \
        --num-recycle ${RECYCLES} \
        --num-models ${MODELS} \
        ${INPUT} ${OUTPUT} &&
    echo "== Finished ${batch_num} predicitons: `date` =="
    free -h
    sync && echo 3 | sudo -n tee /proc/sys/vm/drop_caches || true
    free -h
    cat ${WORKDIR}/predictions/log.txt >> ${WORKDIR}/RB_FloraFoldPrediction_batch${BATCH_ID}.log
    ## Removing any failed MSAs so they're not duplicated.
    rm ${INPUT}/*a3m || true
    sleep 15
done 

## Handling when vRAM wall reached errors - should be fixed with the watcher script now
if ls ${WORKDIR}/*.log 1> /dev/null 2>&1 && grep -q "Killed" ${WORKDIR}/*.log; then
    num_complete=`ls -1 ${OUTPUT}/done/*.done.txt | wc -l`
    num_expected=`ls -1 ${INPUT}/*.a3m | wc -l`
    echo "== ColabFold crashed. Completed ${num_expected}/${num_complete}. `date`  =="
fi

## Should be handled by watcher now instead. Will add sleep command to ensure it's caught and processed. 
sleep 30
## Just in case the watcher silently stops.
echo "== Exporting data to GCP: `date` =="
gsutil -m cp ${OUTPUT}/*.pdb ${GCP_OUTPUT_PATH} || true
gsutil -m cp ${OUTPUT}/*0.json ${GCP_OUTPUT_PATH} || true
gsutil -m cp ${OUTPUT}/*v1.json ${GCP_OUTPUT_PATH} || true

echo "==== Run Complete: `date` ===="