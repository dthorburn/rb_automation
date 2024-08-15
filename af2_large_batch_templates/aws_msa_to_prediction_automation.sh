                            ######################################
                            ##            Parameters            ##
                            ######################################

aws_instance_type=g6.2xlarge
ami_id=ami-06e8a163482ac06a1
target_bucket=s3://rb-florafold/05_screening/toy_data
run_name="asr_wnlrs_b2"
security_key_name=rb-gpufold-l4-1
msas_per_vm=300
vm_per_zone=20
counter=0
start_vm_num=1

                            ######################################
                            ##        Splitting MSAs            ##
                            ######################################

aws s3 cp ${target_bucket}/msas/ . --exclude "*" --include "*log" --recursive
split -l ${msas_per_vm} ./*log
counter=0
for i in ./xa*
do
    echo $i
    counter=`expr ${counter} + 1`
    aws s3 mv $i ${target_bucket}/msas/batch${counter}/
done

zone_counter=0
aws_zones=("us-east-1d" "us-east-2c" "us-west-2a")
num_zones=${#aws_zones[@]}
updated_bucket_path=$(echo "${target_bucket}" | sed -e 's/\//\\\//g')

                            ######################################
                            ##       Launching Instances        ##
                            ######################################

for i in $(seq ${start_vm_num} ${counter})
do 
    current_zone=${aws_zones[`expr $zone_counter + 1`]} 
    aws_region=`echo $current_zone | sed -e "s/[a-z]$//"`
    echo "~~~~ Launching EC2 instance VM${i} in zone ${current_zone}: `date`"

    # Create the temp_entrypoint.sh script
    echo '#!/bin/bash' > temp_entrypoint.sh
    echo "
    sed -i \"81s/\.\//\/home\/ec2-user\//\" /root/run_predictions.sh
    #sed -i \"24s/rb-florafold\/05_screening\/toy_data/${updated_bucket_path}/\" /root/run_predictions.sh
    sed -i \"25s/aws_in.*/aws_in_bucket=${updated_bucket_path}\/01_msas/\" /root/run_predictions.sh
    sed -i \"26s/aws_out.*/aws_in_bucket=${updated_bucket_path}\/02_predictions/\" /root/run_predictions.sh
    sed -i \"34s/batchall/batch${i}/\" /root/run_predictions.sh

    sed -i \"89s/\*log/\/home\/ec2-user\/\*log/\" /root/run_predictions.sh
    sed -i \"66s|./\*.log|/home/ec2-user/\*.log|;94s|./\*.log|/home/ec2-user/\*.log|;106s|./\*.log|/home/ec2-user/\*.log|\" /root/run_predictions.sh    
    sudo chmod +x /root/run_predictions.sh
    sudo chmod 777 /root/
    
    # Create a systemd service to run the script after boot
    echo \"[Unit]
    Description=Run Post Startup Script
    After=network.target

    [Service]
    Type=simple
    ExecStart=/root/run_predictions.sh
    StandardOutput=append:/home/ec2-user/${run_name}_florafold_batch${i}.log
    StandardError=append:/home/ec2-user/${run_name}_florafold_batch${i}.log
    User=ec2-user

    [Install]
    WantedBy=multi-user.target\" > /etc/systemd/system/run_predictions.service

    # Enable and start the systemd service
    sudo systemctl daemon-reload
    sudo systemctl enable run_predictions.service
    sudo systemctl start run_predictions.service
    " >> temp_entrypoint.sh

    aws ec2 run-instances \
        --image-id ${ami_id} \
        --count 1 \
        --instance-type ${aws_instance_type} \
        --key-name ${security_key_name} \
        --network-interfaces '[{"SubnetId":"subnet-0ca24e4a53ef669ec","DeviceIndex":0,"Groups":["sg-027038101fc23a54f"]}]' \
        --launch-template "LaunchTemplateId=lt-0b7e8d39f6e8e3158,Version=3" \
        --region ${aws_region} \
        --placement AvailabilityZone=${current_zone} \
        --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=rb-gpufold-l4-${run_name}-${i}}]" \
        --iam-instance-profile Arn=arn:aws:iam::760901650169:instance-profile/EC2_S3_Access \
        --user-data file://temp_entrypoint.sh 

    rm temp_entrypoint.sh

    # Increment the zone counter every ${vm_per_zone} iterations
    if (( i % vm_per_zone == 0 )); then
        zone_counter=$(( (zone_counter + 1) % num_zones ))
    fi
done

                            ######################################
                            ##     Handling Uncaught Crashes    ##
                            ######################################

in_bucket="s3://florafold/03-screens/asr_working_nlrs/batch2/msas/batch10"
out_bucket=`echo ${in_bucket} | sed -e "s/msas/predictions/"`

aws s3 ls ${out_bucket}/"*"rank_001"*"pdb > temp_file_list
sed -i "s/predictions/msas/g" temp_file_list
sed -i "s/_unrelaxed.*/.a3m/g" temp_file_list
cat temp_file_list | aws s3 mv -I ${in_bucket}/done/
