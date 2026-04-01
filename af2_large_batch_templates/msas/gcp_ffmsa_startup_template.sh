#!/bin/bash
set -e

## Getting env variables from GCP metadata
WORKDIR="/home/miles"
BATCH_ID=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/attributes/batch_id -H "Metadata-Flavor: Google")
ROOT_BUCKET=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/attributes/root_bucket -H "Metadata-Flavor: Google")
INSTANCE_NAME=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/name -H "Metadata-Flavor: Google")
ZONE=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/zone -H "Metadata-Flavor: Google" | awk -F/ '{print $NF}')
LOGFILE="${WORKDIR}/RB_FloraFoldMSAStartup_batch${BATCH_ID}.log"

## Starting logging
exec > >(tee -a "$LOGFILE") 2>&1
echo "====== Startup $(date) ======"

## Run the job
sudo -u miles bash ${WORKDIR}/batch_msa.sh "$BATCH_ID" "$ROOT_BUCKET" || true

## Shutting down the VM
echo "==== Runner complete - Shutting down VM $(date) ===="
gsutil -m cp ${WORKDIR}/*.log ${ROOT_BUCKET}/msas/ || true
gcloud compute instances delete "$INSTANCE_NAME" --zone="$ZONE" --quiet
