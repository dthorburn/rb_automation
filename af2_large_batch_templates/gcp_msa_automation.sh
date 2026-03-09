## Automating the deployment of large GCP MSA VMs
                            ######################################
                            ##            Parameters            ##
                            ######################################
## Adjust zones based on general availability. These have served me well
gcp_image_project_id=""
gcp_project_id=""
gcp_service_account=""
zones=("us-west1-a" "us-east4-a" "us-west4-a" "us-east1-d" "us-central1-c") 
## Do not include the "/msas"  or "/predictions" ending of bucket paths. This is automtically added.
## Do not include the gs:// in the target bucket path.
target_bucket="" 
run_name=""
batches=

zone_counter=0
num_zones=${#zones[@]}
#updated_bucket_path=$(echo "${target_bucket}" | sed -e 's/\//\\\//g')
for i in $(seq 1 ${batches})
do 
    ## Sone information 
    current_zone=${zones[`expr $zone_counter + 1`]} ## This may be a little buggy. 
    echo "~~~~ Launching MSA VM${i} in zone ${current_zone}: `date`"

    # Create the temp_entrypoint.sh script
    echo '#!/bin/bash' > temp_entrypoint.sh
    echo "
    # Modify the run_predictions.sh script
    sed -i \"15s/30/28/\" /home/miles/run_batch_msa.sh
    ## NB, This is the full path in the default file: gs://florafold/03-screens/asr/msas
    sed -i \"21s|florafold/03-screens/asr_working_nlrs|${target_bucket}|\" /home/miles/run_batch_msa.sh
    sed -i \"22s|batch3/msas|batch${i}|\" /home/miles/run_batch_msa.sh
    gsutil -m cp gs://${target_bucket}/*0${i}.fasta /home/miles/input/
    sudo chmod +x /home/miles/run_batch_msa.sh

    # Create a systemd service to run the script after boot
    echo \"[Unit]
    Description=Run Post Startup Script
    After=network.target remote-fs.target

    [Service]
    Type=simple
    ExecStartPre=/bin/sleep 30
    ExecStart=bash -x /home/miles/run_batch_msa.sh
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
    ## Linux GNU sed:
    sed -i 's/^[ \t]*//g' temp_entrypoint.sh
    ## MacOS: 
    #perl -pi -e 's/^[ \t]*//g' temp_entrypoint.sh

    gcloud compute instances create rb-msaflorafold50-${i} \
        --project=${gcp_project_id} \
        --zone=${current_zone} \
        --machine-type=custom-28-184320-ext \
        --network-interface=network-tier=PREMIUM,stack-type=IPV4_ONLY,subnet=default \
        --metadata=^,@^ssh-keys=miles:ecdsa-sha2-nistp256\ AAAAE2VjZHNhLXNoYTItbmlzdHAyNTYAAAAIbmlzdHAyNTYAAABBBP69CEtxDHJuzyNsAVJU88\+nT8\+nvaXamnvMMDBTx8pge5GypN8Vd/lamkXiVc\+EucKwCDLR0LLnnl5fq8QLQkQ=\ google-ssh\ \{\"userName\":\"miles@resurrect.bio\",\"expireOn\":\"2024-08-19T09:01:00\+0000\"\}$'\n'miles:ssh-rsa\ AAAAB3NzaC1yc2EAAAADAQABAAABAG0L\+bB7Fbt7PotObGSQLj1kKDDXySyZ4EPJ5eEwLX0FUNlrFlne71rl/ul3FjQnC4Z0m2v164TSV4ysXd0784Dt/8Tm99gWQr3bFJ4oM8cv87cgN\+OlEFppoJ4wuAhKvKlzoaDhWUtBcAmt/CbD9dGCBLUDGfmCa5\+Don7TlPZ356gcfNxyFDI2aPi5nSfmqiqviiiFRdUYvx6EyHKRiY0uOuz5l4T\+o\+DK2KmFxpm9fAaAoup6pGw5GkxSRaRIuxlYAwaPA8SSzLuPniOUgVSB3OFGLCXfeuROrvK1qIoRIbxBS/DtluKlfwbE0Nei7hKI1kWcVVUvK28S9HsqCO8=\ google-ssh\ \{\"userName\":\"miles@resurrect.bio\",\"expireOn\":\"2024-08-19T09:01:04\+0000\"\} \
        --maintenance-policy=MIGRATE \
        --provisioning-model=STANDARD \
        --service-account=${gcp_service_account} \
        --scopes=https://www.googleapis.com/auth/cloud-platform \
        --tags=http-server,https-server \
        --create-disk=auto-delete=yes,boot=yes,device-name=rb-msacolabfold-${i},image=projects/${gcp_image_project_id}/global/images/rb-msaflorafoldbootdisk-update2,mode=rw,size=2500,type=pd-balanced \
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

