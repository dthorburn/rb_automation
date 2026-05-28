#!/bin/bash
set -e

## Getting env variables from GCP metadata
## MINA, MAXA are distance cutoffs for the pDockQ scoring
## PAEC and DISTC are the PAE and distance cutoffs for the ipSAE scoring.
WORKDIR="/home/miles"
ROOT_BUCKET=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/attributes/root_bucket -H "Metadata-Flavor: Google")
MINA=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/attributes/mina -H "Metadata-Flavor: Google")
MAXA=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/attributes/maxa -H "Metadata-Flavor: Google")
PAEC=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/attributes/paec -H "Metadata-Flavor: Google")
DISTC=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/attributes/distc -H "Metadata-Flavor: Google")
DELETE=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/attributes/delete -H "Metadata-Flavor: Google")
INSTANCE_NAME=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/name -H "Metadata-Flavor: Google")
ZONE=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/zone -H "Metadata-Flavor: Google" | awk -F/ '{print $NF}')
LOGFILE="${WORKDIR}/RB_FloraFoldScoring_batch${BATCH_ID}.log"

## Starting logging
exec > >(tee -a "$LOGFILE") 2>&1
echo "====== Startup $(date) ======"

## These are only for debugging, this is explicitely printed in the script below.
#echo "Run: $RUN_NAME"
#echo "Batch: $BATCH_ID"

## Run the job
sudo -u miles bash ${WORKDIR}/process_interactions.sh "$ROOT_BUCKET" "$MINA" "$MAXA" "$PAEC" "$DISTC" || true

## Shutting down the VM
echo "==== Runner complete - Shutting down VM ===="
gsutil -m cp ${WORKDIR}/*.log ${ROOT_BUCKET}/predictions || true
if [ "${DELETE}" = "true" ]; then
    gcloud compute instances delete "$INSTANCE_NAME" --zone="$ZONE" --quiet
else
    gcloud compute instances stop "$INSTANCE_NAME" --zone="$ZONE" --quiet
fi
