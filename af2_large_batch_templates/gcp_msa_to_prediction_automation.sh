## Automating the deployment of Nvidia L4 Screens
                            ######################################
                            ##            Parameters            ##
                            ######################################
gcp_project_id=""
## Adjust zones based on general availability. These have served me well
zones=("us-central1-a" "us-west1-b" "us-west4-a" "us-east1-d" "us-east4-a")
## Do not include the "/msas"  or "/predictions" ending of bucket paths. This is automtically added.
target_bucket="gs://florafold/03-screens/asr_working_nlrs/batch2"
msas_per_vm=300
vm_per_zone=16
run_name="asr_wnlrsb4"


## Or manually do this step if the path is incorrect. 
gsutil -m cp ${target_bucket}/msas/"*"log ./
split -l ${msas_per_vm} ./*log

## Splitting MSAs into managable chunks
counter=0
for i in ./xa*
do
    echo $i
    counter=`expr ${counter} + 1`
    cat $i | gsutil -m mv -I ${target_bucket}/msas/batch${counter}/ &&
        rm $i
done

zone_counter=0
num_zones=${#zones[@]}
## Automating the change of buckets
updated_bucket_path=$(echo "${target_bucket}" | sed -e 's/\//\\\//g')

## Known issue: If there is no zone availability, it will cycle through 12 attempts each time iterating up before moving to the next one. 
#for i in {8..8..1}
for i in $(seq 1 ${counter})
do 
    ## Sone information 
    current_zone=${zones[`expr $zone_counter + 1`]} ## This may be a little buggy. 
    echo "~~~~ Launching L4 VM${i} in zone ${current_zone}: `date`"

    # Create the temp_entrypoint.sh script
    echo '#!/bin/bash' > temp_entrypoint.sh
    echo "
    # Modify the run_predictions.sh script
    #sed -i \"28s/1/0/\" /home/miles/run_predictions.sh
    sed -i \"81s/\.\//\/home\/miles\//\" /home/miles/run_predictions.sh
    ## NB, This is the full path in the default file: gs://florafold/03-screens/asr/msas
    sed -i \"33s//gs:\/\/florafold\/03-screens\/asr/${updated_bucket_path}/\" /home/miles/run_predictions.sh
    sed -i \"34s//gs:\/\/florafold\/03-screens\/asr/${updated_bucket_path}/\" /home/miles/run_predictions.sh
    sed -i \"37s/batch1/batch${i}/\" /home/miles/run_predictions.sh
    sudo chmod +x /home/miles/run_predictions.sh

    # Create a systemd service to run the script after boot
    echo \"[Unit]
    Description=Run Post Startup Script
    After=network.target

    [Service]
    Type=simple
    ExecStart=/home/miles/run_predictions.sh
    StandardOutput=append:/home/miles/${run_name}_florafold_batch${i}.log
    StandardError=append:/home/miles/${run_name}_florafold_batch${i}.log
    User=miles
    Group=miles
    Environment=PATH=\"/usr/local/cuda/bin:/opt/conda/bin:/opt/conda/condabin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games\"

    [Install]
    WantedBy=multi-user.target\" > /etc/systemd/system/run_predictions.service

    # Enable and start the systemd service
    sudo systemctl daemon-reload
    systemctl enable run_predictions.service
    systemctl start run_predictions.service
    " >> temp_entrypoint.sh
    sed -i 's/^[ \t]*//g' temp_entrypoint.sh

    gcloud compute instances create rb-gpufold-l4-${run_name}-${i} \
        --project=${gcp_project_id} \
        --zone=${current_zone} \
        --machine-type=g2-custom-8-49152 \
        --network-interface=network-tier=PREMIUM,stack-type=IPV4_ONLY,subnet=default \
        --maintenance-policy=TERMINATE \
        --provisioning-model=STANDARD \
        --service-account=382883280368-compute@developer.gserviceaccount.com \
        --scopes=https://www.googleapis.com/auth/cloud-platform \
        --accelerator=count=1,type=nvidia-l4 \
        --create-disk=auto-delete=yes,boot=yes,device-name=rb-gpufold-l4-${run_name}-${i},image=projects/${gcp_project_id}/global/images/rb-colabfoldbootdisk-florafold1,mode=rw,size=100,type=projects/${gcp_project_id}/zones/us-west1-c/diskTypes/pd-balanced \
        --no-shielded-secure-boot \
        --shielded-vtpm \
        --shielded-integrity-monitoring \
        --labels=goog-ec-src=vm_add-gcloud \
        --reservation-affinity=any \
        --metadata-from-file startup-script=temp_entrypoint.sh
    rm temp_entrypoint.sh ## Just being paranoid

    # Increment the zone counter every ${vm_per_zone} iterations
    if (( i % vm_per_zone == 0 )); then
        zone_counter=$(( (zone_counter + 1) % num_zones ))
    fi
done

