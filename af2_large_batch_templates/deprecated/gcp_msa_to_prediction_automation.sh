## Automating the deployment of Nvidia L4 Screens
                            ######################################
                            ##            Parameters            ##
                            ######################################
## Adjust zones based on general availability. These have served me well
gcp_project_id=""
gcp_service_account=""
zones=("us-west4-a") 
zones=("us-west1-a" "us-east4-a" "us-west4-a" "us-east1-d" "us-central1-c") 
## Do not include the "/msas"  or "/predictions" ending of bucket paths. This is automtically added.
target_bucket="gs://rb-interactions/soleracea_peffusa"
run_name="soleracea_peffusa"
msas_per_vm=200 ## 255 top
vm_per_zone=40

## Options to override the iteration of VMs from 1 to x and custom vms per zone to fill out quota.
custom_vms_per_zone=0 ## boolean use custom vms per zone array
#vm_array=(21 18 11 6 6) #158
vm_array=(21) #158
start_vm_num=1
end_vm_num=148
zone_counter=0
num_zones=${#zones[@]}
updated_bucket_path=$(echo "${target_bucket}" | sed -e 's/\//\\\//g')
updated_run_name=$(echo "${run_name}" | sed -e "s/_/-/g")

## Splitting MSAs into managable chunks
## Or manually do this step if the path is incorrect. 
gsutil -m cp ${target_bucket}/msas/"*"log ./
split -l ${msas_per_vm} file_list

#gsutil -m ls ${target_bucket}/msas/batch"*"/"*"a3m > ${run_name}_files.log
gsutil -m ls ${target_bucket}/"*"a3m > ${run_name}_files.log
split -l ${msas_per_vm} ${run_name}_files.log
end_vm_num=0
for i in ./x*
do
    echo $i
    end_vm_num=`expr ${end_vm_num} + 1`
    cat $i | gsutil -m mv -I ${target_bucket}/msas/batch${end_vm_num}/
done

#for i in 126 127 $(seq 129 134) 137 138
for i in $(seq ${start_vm_num} ${end_vm_num})
do 
    ## Setting zone information
    current_zone=${zones[$zone_counter]} ## This may be a little buggy. 
    if [ $i == $start_vm_num ]
    then
        inter_counter=0
    fi
    inter_counter=`expr $inter_counter + 1`
    
    if [ $custom_vms_per_zone == 1 ]
    then
        vm_per_zone=${vm_array[$zone_counter]}
    fi
    echo "~~~~ Launching L4 VM${i} (${inter_counter} of ${vm_per_zone}) in zone ${current_zone}: `date`"

    # Create the temp_entrypoint.sh script
    echo '#!/bin/bash' > temp_entrypoint.sh
    echo "
    # Modify the run_predictions.sh script
    # sed -i \"28s/1/0/\" /home/miles/run_predictions.sh
    sed -i \"81s/\.\//\/home\/miles\//\" /home/miles/run_predictions.sh
    ## NB, This is the full path in the default file: gs://florafold/03-screens/asr/msas
    sed -i \"33s/gs:\/\/florafold\/03-screens\/asr/${updated_bucket_path}/\" /home/miles/run_predictions.sh
    sed -i \"34s/gs:\/\/florafold\/03-screens\/asr/${updated_bucket_path}/\" /home/miles/run_predictions.sh
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

    gcloud compute instances create rb-gpufold-l4-${updated_run_name}-${i} \
        --project=${gcp_project_id} \
        --zone=${current_zone} \
        --machine-type=g2-custom-8-49152 \
        --network-interface=network-tier=PREMIUM,stack-type=IPV4_ONLY,subnet=default \
        --maintenance-policy=TERMINATE \
        --provisioning-model=STANDARD \
        --service-account=${gcp_service_account} \
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
    if (( inter_counter % vm_per_zone == 0 )); then
        zone_counter=$(( (zone_counter + 1) % num_zones ))
        inter_counter=0
    fi
done

                            ######################################
                            ##     Handling Uncaught Crashes    ##
                            ######################################

in_bucket="gs://florafold/03-screens/scn_working_nlrs/msas/batch107"
out_bucket=`echo ${in_bucket} | sed -e "s/msas/predictions/"`

gsutil -m ls ${out_bucket}/"*"rank_001"*"pdb > temp_file_list
sed -i "s/predictions/msas/g" temp_file_list
sed -i "s/_unrelaxed.*/.a3m/g" temp_file_list
cat temp_file_list | gsutil -m mv -I ${in_bucket}/done/

