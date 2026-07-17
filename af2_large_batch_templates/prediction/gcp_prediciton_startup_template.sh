#!/bin/bash
set -e

## Getting env variables from GCP metadata
WORKDIR="/home/miles"
BATCH_ID=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/attributes/batch_id -H "Metadata-Flavor: Google")
MODELS=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/attributes/models -H "Metadata-Flavor: Google")
RECYCLES=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/attributes/recycles -H "Metadata-Flavor: Google")
ROOT_BUCKET=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/attributes/root_bucket -H "Metadata-Flavor: Google")
INSTANCE_NAME=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/name -H "Metadata-Flavor: Google")
ZONE=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/zone -H "Metadata-Flavor: Google" | awk -F/ '{print $NF}')
LOGFILE="${WORKDIR}/RB_FloraFoldStartup_batch${BATCH_ID}.log"

## Starting logging
exec > >(tee -a "$LOGFILE") 2>&1
echo "====== Startup $(date) ======"

## These are only for debugging, this is explicitely printed in the script below.
#echo "Run: $RUN_NAME"
#echo "Batch: $BATCH_ID"

# Install drivers if needed
if ! nvidia-smi &> /dev/null; then
    echo "==== Installing NVIDIA drivers ===="
    /opt/deeplearning/install-driver.sh
fi

## Making file executable
#chmod +x /home/miles/run_predictions.sh
#chmod +x /home/miles/testing_watcher.sh

## Set up the watcher
nohup sudo -u miles bash ${WORKDIR}/watcher.sh ${BATCH_ID} ${ROOT_BUCKET} &
## This took a moment to get up and running, so just trying to sync things up in the logs
sleep 5 

## Run the job
sudo -u miles bash ${WORKDIR}/run_predictions.sh "$BATCH_ID" "$MODELS" "$RECYCLES" "$ROOT_BUCKET" || true
cp ${WORKDIR}/predictions/log.txt ${WORKDIR}/RB_FloraFoldPrediction_batch${BATCH_ID}.log
#pkill -9 -f watcher.sh

## Shutting down the VM
echo "==== Runner complete - Shutting down VM ===="
gsutil -m cp ${WORKDIR}/*.log ${ROOT_BUCKET}/predictions || true
gcloud compute instances delete "$INSTANCE_NAME" --zone="$ZONE" --quiet
