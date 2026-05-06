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
bash

## 3. Make directories
mkdir -p predictions
mkdir -p scored

## 4. create conda env - might need updating if this branch is merged into main
curl -O https://raw.githubusercontent.com/dthorburn/rb_automation/mt-dev/af2_large_batch_templates/scoring/af2_scoring.yml
curl -O https://raw.githubusercontent.com/dthorburn/rb_automation/mt-dev/af2_large_batch_templates/scoring/scoring_interactions.py
curl -O https://raw.githubusercontent.com/dthorburn/rb_automation/mt-dev/af2_large_batch_templates/scoring/secondary_structure.py
curl -O https://raw.githubusercontent.com/dthorburn/rb_automation/mt-dev/af2_large_batch_templates/scoring/merge_results.py
conda create -n af2_scoring --file af2_scoring.yml

