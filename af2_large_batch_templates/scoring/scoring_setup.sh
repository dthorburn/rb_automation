## 1. Check for install requirements
sudo apt-get update && sudo apt-get install -y \
  git \
  gcc \
  g++ \
  make \
  wget \
  curl \

## 2. Install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda init bash

## 3. Make directories
mkdir -p predictions
mkdir -p scored

## 4. create conda env
git 
