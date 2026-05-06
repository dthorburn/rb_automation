#!/bin/zsh
set -e

## Parameters - ALL ARE MANDATORY
gcp_image_project_id=""
gcp_project_id=""
gcp_service_account=""
batches=
root_bucker=""
run_name=""
## General baselines for our current use case
models=3
recycles=2
region_quota=32
start_vm_num=1

## NOTE: If deployment of large jobs are going to be staggered, you'll need to adjust the ${batches}
##       and ${start_vm_num} values to reflect which batches to run in each execution.
##       e.g., start_vm_num=40; batches=80 will launch all batches between 40-80. 
## I tried to make this one shorter, but I thought I would also try to automate deployment and failed zones/region caps. 
############################################
##        Production - Don't change       ##
############################################
## Strips zone identified to get region id.
get_region() {
    echo "${1%-*}"
}
## Zone list is exhaustive for G2-Nvidia-L4 instances.
## Zone list is interleaved to try and first avoid region caps.
zones=(
    "us-central1-a" "us-east1-b" "us-east4-a" "us-west1-a" "us-west4-a" 
    "us-central1-b" "us-east1-c" "us-east1-d" "us-west1-b" "us-west4-c" 
    "us-central1-c" "us-east4-b" "us-west1-c" 
    "us-central1-f" "us-east4-c"
)
zones=("us-central1-a" "us-east1-b" "us-east4-a" "us-west1-a" "us-west4-a" "us-central1-b" "us-east1-c" "us-east1-d" "us-west1-b" "us-west4-c" "us-central1-c" "us-east4-b" "us-west1-c" "us-central1-f" "us-east4-c")
## bash is 0-indexed, and zsh is 1-indexed. Change if needed. 
zone_index=1
num_zones=${#zones[@]}
for batch_id in $(seq ${start_vm_num} ${batches})
do
    launched=false
    attempts=0
    while [[ "${launched}" = false && ${attempts} -lt ${num_zones} ]]
    do
        zone="${zones[${zone_index}]}"
        region=$(get_region "${zone}")
        used_var="region_used_${region//-/_}"
        used=${!used_var:-0}
        if (( used + 1 >= region_quota )); then
            echo "Skipping ${zone} — region ${region} at quota (${used}/${region_quota} GPUs used)"
            zone_index=$(( (zone_index + 1) % num_zones ))
            attempts=$(( attempts + 1 ))
        else
            echo "Attempting to launch batch ${batch_id} in zone ${zone}..."
            if gcloud compute instances create rb-gpucfold-l4-${run_name}-${batch_id} \
                --source-instance-template=projects/${gcp_image_project_id}/global/instanceTemplates/florafold-l4-template \
                --project=${gcp_project_id} \
                --service-account=${gcp_service_account} \
                --zone="${zone}" \
                --metadata=startup-script=/home/miles/startup.sh,batch_id=${batch_id},models=${models},recycles=${recycles},root_bucket=${root_bucket} 2>&1; then
                eval "${used_var}=$(( used + 1 ))"
                echo "Launched batch ${batch_id} in ${zone} — region ${region} now at $(( used + 1 ))/${region_quota} GPUs"
                launched=true
            else
                echo "Failed in ${zone}, trying next..."
                zone_index=$(( (zone_index + 1) % num_zones ))
                attempts=$(( attempts + 1 ))
            fi
        fi
    done
    if [[ "${launched}" = false ]]; then
        echo "ERROR: Failed to launch batch ${batch_id} in all zones. Exiting."
        break
    fi
done
