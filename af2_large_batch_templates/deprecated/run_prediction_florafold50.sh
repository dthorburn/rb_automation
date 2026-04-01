#!/usr/bin/env bash
## Maybe move to .bashrc
export PATH="/home/resurrect/resurrectbio/01_projects/03_colabfold/01_installing/localcolabfold152/colabfold-conda/bin:/home/resurrect/tools/google-cloud-sdk/bin:/home/resurrect/conda/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin"
CUDNN_PATH=$(dirname $(python -c "import nvidia.cudnn;print(nvidia.cudnn.__file__)"))
export LD_LIBRARY_PATH=$CUDNN_PATH/lib:$CONDA_PREFIX/lib/:$LD_LIBRARY_PATH

case $1 in
 -[h?] | --help)
        printf "Usage:"
        exit 0;;
esac

## Paramaters
models=5
recycles=5

import_from_cloud=0
export_to_cloud=0
terminate_on_complete=0
remove_files_on_complete=0

## Directories
cloud_service="aws" ## AWS | GCP - GCP support will end in August 2024.
aws_folder_name="rb-kagome"
aws_in_bucket="s3://${aws_folder_name}/msas"
aws_out_bucket="s3://${aws_folder_name}/predictions"

gcp_folder_name="pathogen-pdb-library"
gcp_in_bucket="gs://${gcp_folder_name}/01_msas"
gcp_out_bucket="gs://${gcp_folder_name}/02_pdbs"

input="/mnt/sdb1/kagome/02_msas"
output="/mnt/sdb1/kagome/03_predicitons"
batch_name="batchall"

## Commands - don't change
## Options for input source
if [ ${import_from_cloud} == 1 ]
then
        rm ${input}/*
	if [ ${cloud_service} == "aws" ]
	then
        	echo "~~ Importing MSAs from AWS: `date`"
		aws s3 cp --recursive ${aws_in_bucket}/${batch_name}/ ${input}/ --exclude "*" --include "*a3m"
	else
	        echo "~~ Importing MSAs from GCP: `date`"
        	gsutil -m cp ${gcp_in_bucket}/${batch_name}/*.a3m ${input}/
	fi
else
        echo "~~ Using files in ${input}: `ls -1 ${input}/*.a3m | wc -l` MSA files present: `date`"
fi

## Running colabfold AF2 predicitons
rm ${output}/*
echo "~~ Starting preditcions: `date`"
colabfold_batch \
	--disable-unified-memory \
        --num-recycle ${recycles} \
        --num-models ${models} \
        ${input} ${output} &&
        echo "~~ Successfully finished predicitons: `date`"

export LD_LIBRARY_PATH=""

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
if [ ${export_to_cloud} == 1 ]
then
	if [ ${cloud_service} == "aws" ]
	then
	        echo "~~ Exporting MSAs to AWS: `date`"
		aws s3 cp --recursive `pwd` ${aws_out_bucket}/${batch_name}/ --exclude "*" --include "*log"
		aws s3 cp --recursive ${output}/ ${aws_out_bucket}/${batch_name}/ --exclude "*" --include "*pdb"
		aws s3 cp --recursive ${output}/ ${aws_out_bucket}/${batch_name}/ --exclude "*" --include "*json"
	else
	        echo "~~ Exporting MSAs to GCP: `date`"
        	gsutil -m mv ./*.log ${gcp_out_bucket}/
        	gsutil -m mv ${output}/*.pdb ${gcp_out_bucket}/${batch_name}/
	        gsutil -m mv ${output}/*0.json ${gcp_out_bucket}/${batch_name}/
        	gsutil -m mv ${output}/*v1.json ${gcp_out_bucket}/${batch_name}/
	fi
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
