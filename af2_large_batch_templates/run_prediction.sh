#!/usr/bin/env bash
## Maybe move to .bashrc
export PATH="/home/miles/localcolabfold/colabfold-conda/bin:/home/miles/localcolabfold/conda/bin/:$PATH"

case $1 in
 -[h?] | --help)
        printf "Usage:"
        exit 0;;
esac

## Noticed that sometimes the drivers don't install correctly. So here is an attept to automate fixing that problem
nvs_test=`nvidia-smi`
if  echo $nvs_test | grep -q "has failed"
then
  echo "~~ Reinstalling NVIDIA Drivers: `date`"
  sudo /opt/deeplearning/install-driver.sh
else
  echo "~~ NVIDIA Drivers Fine: `date`"
fi


## Paramaters
models=1
recycles=3

import_from_gcp=1
export_to_gcp=1
terminate_on_complete=1
remove_files_on_complete=1

## Directories
gcp_folder_name="gmax/hgeffectors"
gcp_in_bucket="gs://${gcp_folder_name}/msas"
gcp_out_bucket="gs://${gcp_folder_name}/pdbs"
input="/home/miles/msas"
output="/home/miles/predictions"
batch_name="batch1"

## Commands - don't change
## Options for input source
if [ ${import_from_gcp} == 1 ]
then
        echo "~~ Importing MSAs from GCP: `date`"
        rm ${input}/*
        gsutil -m cp ${gcp_in_bucket}/${batch_name}/*.a3m ${input}/
else
        echo "~~ Using files in ${input}: `ls -1 ${input}/*.a3m | wc -l` MSA files present: `date`"
fi

## Running colabfold AF2 predicitons
rm ${output}/*
echo "~~ Starting preditcions: `date`"
colabfold_batch \
        --num-recycle ${recycles} \
        --num-models ${models} \
        ${input} ${output} &&
        echo "~~ Successfully finished predicitons: `date`"

## Handling when memory wall reached errors
if cat ./*.log | grep -q "Killed"
then
        num_complete=`ls -1 ${output}/*.a3m | wc -l`
        num_expected=`ls -1 ${input}/*.pdb | wc -l`
        echo "~~ ColabFold crashed. Completed ${num_expected}/${num_complete} Moving finished input: `date`"
        for df in ./predictions/*.done.txt
        do
          newname=`basename ${df} | sed "s/.done.txt/.a3m/"`
          if [ -f ./msas/${newname} ]
          then
            echo "Moving ${newname}"
            rm ./msas/${newname}
            gsutil -m mv ${gcp_in_bucket}/${batch_name}/${newname} ${gcp_in_bucket}/${batch_name}/done/${newname}
          fi
        done
fi

## Exporting to GCP options. 
if [ ${export_to_gcp} == 1 ]
then
        echo "~~ Exporting MSAs to GCP: `date`"
        gsutil -m mv ./*.log ${gcp_out_bucket}/
        gsutil -m mv ${output}/*.pdb ${gcp_out_bucket}/${batch_name}/
        gsutil -m mv ${output}/*0.json ${gcp_out_bucket}/${batch_name}/
        gsutil -m mv ${output}/*v1.json ${gcp_out_bucket}/${batch_name}/
fi

if [ ${remove_files_on_complete} == 1 ]
then
        echo "~~ Removing extra files from run: `date`"
        rm ${input}/*
        rm ${output}/*
        rm ./*.log
fi

if [ ${terminate_on_complete} == 1 ]
then
        echo "~~ Terminating VM: `date`"
        sudo shutdown -h now
fi