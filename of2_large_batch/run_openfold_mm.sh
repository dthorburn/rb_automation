#!/bin/bash -e
export PATH="/home/miles/localopenfold/openfold:/usr/local/cuda/bin:/home/miles/localopenfold/conda/envs/openfold_env/bin:/opt/conda/condabin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games"

case $1 in
 -[h?] | --help)
        printf "Usage:"
        exit 0;;
esac

## NVIDIA driver test
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
batch_name="batchall"
gcp_folder_name="taestivum/alternative-moa"
gcp_in_bucket="gs://${gcp_folder_name}/01_msas"
gcp_out_bucket="gs://${gcp_folder_name}/02_pdbs"
input="/home/miles/msas"
output="/home/miles/predictions"
fasta_dir="/home/miles/fastas"
of_dir="/home/miles/localopenfold/openfold"

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

## Running OpenFold2 predicitons
rm ${output}/*
echo "~~ Starting preditcions: `date`"

cd ${of_dir}
python3 run_pretrained_openfold.py \
	${fasta_dir} \
	\
    --use_precomputed_alignments ${input} \
    --jackhmmer_binary_path /home/miles/localopenfold/conda/envs/openfold_env/bin/jackhmmer \
    --hhblits_binary_path /home/miles/localopenfold/conda/envs/openfold_env/bin/hhblits \
    --hmmsearch_binary_path /home/miles/localopenfold/conda/envs/openfold_env/bin/hmmsearch \
    --hmmbuild_binary_path /home/miles/localopenfold/conda/envs/openfold_env/bin/hmmbuild \
    --kalign_binary_path /home/miles/localopenfold/conda/envs/openfold_env/bin/kalign \
    --config_preset "model_1_multimer_v3" \
    --model_device "cuda:0" \
    --output_dir ${output} 

## Handling when memory wall reached errors
if cat ./*.log | grep -q "Killed"
then
        num_complete=`ls -1 ${output}/*.a3m | wc -l`
        num_expected=`ls -1 ${input}/*.pdb | wc -l`
        echo "~~ OpenFold crashed. Completed ${num_expected}/${num_complete}. Moving finished input: `date`"
        #for df in ./predictions/*.done.txt
        #do
        #  newname=`basename ${df} | sed "s/.done.txt/.a3m/"`
        #  if [ -f ./msas/${newname} ]
        #  then
        #    echo "Moving ${newname}"
        #    rm ./msas/${newname}
        #    gsutil -m mv ${gcp_in_bucket}/${batch_name}/${newname} ${gcp_in_bucket}/${batch_name}/done/${newname}
        #  fi
        #done
fi