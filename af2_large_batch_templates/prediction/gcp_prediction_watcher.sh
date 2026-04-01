#!/usr/bin/env bash
set -euo pipefail

## Taking the env variables
BATCH_ID=$1
INPUT_BUCKET=$2
#OUTPUT_BUCKET=$3

WORKDIR="/home/miles"
INPUT="${WORKDIR}/msas"
OUTPUT="${WORKDIR}/predictions"
GCP_INPUT_PATH="${INPUT_BUCKET}/msas/batch${BATCH_ID}"
GCP_OUTPUT_PATH="${INPUT_BUCKET}/predictions/batch${BATCH_ID}"
GCP_FAIL_PATH="${INPUT_BUCKET}/preempted_jobs"

## Preparing a cleaup on SIGTERM
cleanup() {
    echo "==== Caught SIGTERM at $(date) ===="
    touch "${WORKDIR}/vm_interupt_batch${BATCH_ID}"
    gsutil -m cp "${WORKDIR}/vm_interupt_batch${BATCH_ID}" "${GCP_FAIL_PATH}" || true
    gsutil -m cp ${WORKDIR}/*.log "${GCP_OUTPUT_PATH}" || true
    echo "==== Cleanup done, exiting ===="
    exit 0
}
## Failing SIGTERM/SIGINT, checking preempted status of the VM.
check_preemption() {
    PREEMPTED=$(curl -s -f \
      -H "Metadata-Flavor: Google" \
      http://metadata.google.internal/computeMetadata/v1/instance/preempted || echo "FALSE")

    if [[ "$PREEMPTED" == "TRUE" ]]; then
        echo "==== Preemption detected via metadata ===="
        cleanup
    fi
}
## Run the cleaup command on catching a SIGTERM/SIGINT
trap cleanup SIGTERM SIGINT

## Setting up a log file and starting. 
LOGFILE="${WORKDIR}/RB_FloraFoldWatcher_batch${BATCH_ID}.log"
exec > >(tee -a "$LOGFILE") 2>&1

echo "====== Starting watcher at $(date) ======"
#mkdir -p "${INPUT}/done"
#mkdir -p "${OUTPUT}/done"

while true; do
    check_preemption
    ## Find done markers
    for done_file in "${OUTPUT}"/*.done.txt
    do
        [ -e "$done_file" ] || continue  ## skip if none found
        base_name=$(basename "$done_file" .done.txt)
        echo "== Processing completed model: ${base_name} at $(date) =="

        ## Move output files
        gsutil -m cp ${OUTPUT}/${base_name}*pdb "${GCP_OUTPUT_PATH}/"
        gsutil -m cp ${OUTPUT}/${base_name}*json "${GCP_OUTPUT_PATH}/"
        gsutil -m cp "$done_file" "${GCP_OUTPUT_PATH}/"
        mv ${OUTPUT}/${base_name}* ${OUTPUT}/done || continue

        ## Move input to done folder
        input_file="${INPUT}/${base_name}.a3m"
        if [ -f "$input_file" ]; then
            mv "$input_file" "${INPUT}/done/"
            gsutil -m mv "${GCP_INPUT_PATH}/${base_name}.a3m" "${GCP_INPUT_PATH}/done/"
        fi
    done
    ## Check every 10 seconds - needs to be considerably more often than the ~30 shutdown signal. 
    sleep 10 &
    wait $!
done
