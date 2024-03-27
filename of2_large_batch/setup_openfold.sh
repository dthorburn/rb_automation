#!/bin/bash -e

## Reminder: Add the ll alias to .bashrc
## Reminder: Remove conda init from .bashrc
## Reminder: Remove conda installation at /opt

## Paramaters
DOWNLOAD_DBS=0

## Establishing dir paths
CURRENTPATH=`pwd`
OFDIR="${CURRENTPATH}/localopenfold"

## Starting Installation
echo "~~ Installing OpenFold: `date`"
## Checking for wget
type wget 2>/dev/null || { echo "wget is not installed. Please install it using apt or yum." ; exit 1 ; }

## Noticed that sometimes the drivers don't install correctly. So here is an attept to automate fixing that probl>
echo "~~ Checking GCP NVIDIA drivers: `date`"
nvs_test=`nvidia-smi`
if  echo $nvs_test | grep -q "has failed"
then
  echo "~~ Reinstalling NVIDIA Drivers: `date`"
  sudo /opt/deeplearning/install-driver.sh
else
  echo "~~ NVIDIA Drivers Fine: `date`"
fi

## Checking git and installing ninja-build
sudo apt install ninja-build
mkdir -p ${OFDIR}
cd ${OFDIR}
echo "~~ Setting up Mamba: `date`"
sudo apt-get update
type git 2>/dev/null || { echo "Installing git." ; sudo apt-get install git-all ; }
git version

## Setting up paths and installing conda/mamba
echo "~~ Setting up Mamba: `date`"
wget -q -P . https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash ./Mambaforge-Linux-x86_64.sh -b -p "${OFDIR}/conda"
rm Mambaforge-Linux-x86_64.sh

source "${OFDIR}/conda/etc/profile.d/conda.sh"
export PATH="${OFDIR}/conda/condabin:${PATH}"
conda update -n base conda -y
conda init
conda config --set auto_activate_base false

## Setting up OpenFold
echo "~~ Setting up OpenFold: `date`"
git clone https://github.com/aqlaboratory/openfold.git
cd ./openfold
mamba env create -n openfold_env -f environment.yml
conda activate openfold_env
bash scripts/install_third_party_dependencies.sh

## Sanity check - might need to be run when you have a GPU attached. Can you switch GPUs after? 
echo $LD_LIBRARY_PATH
export LIBRARY_PATH=$CONDA_PREFIX/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

## Download weights and templates
echo "~~ Downloading Weights: `date`"
cd ${OFDIR}
bash scripts/download_openfold_params.sh openfold/resources 
bash scripts/download_pdb_mmcif.sh data/
bash scripts/download_pdb_seqres.sh data/


if [ ${DOWNLOAD_DBS} == 0 ]
then
	echo "~~ Downloading DBs: `date`"
	## Why is the path different than the weights above?
	echo "~~~~ AF2 MultimverV3 Weights: `date`"
	bash scripts/download_alphafold_params.sh openfold/resources
	echo "~~~~ Downloading MMSeqs_DBs: `date`"
	bash scripts/download_mmseqs_dbs.sh data/ \
		&& echo "~~~~ Unpacking MMSeqs_DBs: `date`" \
		&& bash scripts/prep_mmseqs_dbs.sh data/ \
		&& echo "~~~~ Unpacked MMSeqs_DBs: `date`"
fi


echo "~~ Finished OpenFold Installation: `date`"
