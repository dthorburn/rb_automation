#!/bin/bash
set -e

## Parameters - ALL ARE MANDATORY
gcp_image_project_id=""
gcp_project_id=""
gcp_service_account=""
batches=
root_bucket=""
run_name=""

## Zone list is not close to exhaustive for CPU only jobs. Generally, this isn't an issue.
zones=("us-east4-a" "us-west4-a" "us-east1-b" "us-central1-c" "us-central1-a" "us-west1-a")
for batch_id in {1..${batches}}
do
    launched=false
    for zone in "${zones[@]}"
    do
        echo "== Attempting to launch batch${batch_id} in zone ${zone}..."
        if gcloud compute instances create rb-msaflorafold-${run_name}-${batch_id} \
            --source-instance-template=projects/${gcp_image_project_id}/global/instanceTemplates/florafold-msa-template \
            --project=${gcp_project_id} \
            --service-account=${gcp_service_account} \
            --zone=${zone} \
            --metadata=startup-script=/home/miles/startup.sh,batch_id=${batch_id},root_bucket=${root_bucket} 2>&1; then
            echo "== Successfully launched batch${batch_id} in zone ${zone}"
            launched=true
            break
        else
            echo "WARN: Failed to launch in ${zone}, trying next zone..."
        fi
    done

    if [ "$launched" = false ]; then
        echo "ERROR: Failed to launch batch${batch_id} in all zones. Exiting."
        exit 1
    fi
done