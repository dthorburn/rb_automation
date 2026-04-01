## Steps to install localcolabfold on a VM instance for template.
## Requirements:
## - 28 CPU
## - 1.2TB Storage

## 1. Check for install requirements
sudo apt-get update && sudo apt-get install -y \
  git \
  gcc \
  g++ \
  make \
  aria2 \
  wget \
  curl \
  bzip2 \
  ca-certificates \
  gnupg \
  lsb-release
sudo apt-get install  -y seqkit
git --version
curl --version
wget --version
gcc --version

## 2. Make directories
mkdir -p input/split
mkdir -p msas
mkdir -p dbs/florafold
mkdir -p dbs/alphafold

## 3. Install pixi & localcolabfold
curl -fsSL https://pixi.sh/install.sh | sh
git clone https://github.com/yoshitakamo/localcolabfold.git
cd localcolabfold
pixi install && pixi run setup
cd ~/

## Downloading DBs
nano setup_databases.sh ## gcp_ffmsa_db_downloads.sh
sudo chmod +x setup_databases.sh
touch /home/miles/dbs/alphafold/SKIP_TEMPLATES 
touch /home/miles/dbs/alphafold/PDB_MMCIF_READY 
touch /home/miles/dbs/alphafold/COLABDB_READY 
touch /home/miles/dbs/alphafold/PDB_READY 
touch /home/miles/dbs/alphafold/PDB100_READY 
/home/miles/setup_databases.sh /home/miles/dbs/alphafold

gsutil -m cp gs://florafold/06_prod_db/v1/flora"*" /home/miles/dbs/florafold
cd /home/miles/dbs/florafold
#mmseqs tsv2exprofiledb \
#        florafolddb_reduced50 \
#        florafolddb_reduced50 &&
#        echo "~~~~ Finished expanding profile: `date`"
ln -s florafolddb_reduced50_h florafolddb_reduced50_seq_h
ln -s florafolddb_reduced50_h.dbtype florafolddb_reduced50_seq_h.dbtype
ln -s florafolddb_reduced50_h.index florafolddb_reduced50_seq_h.index
mkdir -p tmp1
mmseqs createindex florafolddb_reduced50 tmp1 --remove-tmp-files 1 &&
        echo "~~~~ Finished indexing: `date`"
rm -r tmp1

## 4. Create prediction runner and watcher scripts.
nano batch_msa.sh
nano startup.sh
sudo chmod +x batch_msa.sh
sudo chmod +x startup.sh

## 5. Get toy data and test
export PATH="/home/miles/localcolabfold/.pixi/envs/default/bin/:${PATH}"
cd ~/
WORKDIR="/home/miles"
gsutil -m cp gs://rb-interactions/gmax-psojae/toy_data/fastas/"*"fasta ${WORKDIR}/input
colabfold_search \
  --db1 ${WORKDIR}/dbs/alphafold/uniref30_2302_db \
  --db3 ${WORKDIR}/dbs/florafold/florafolddb_reduced50 \
  --threads 26 \
  --db-load-mode 2 \
  --use-env 1 \
  ${WORKDIR}/input ${WORKDIR}/dbs/alphafold ${WORKDIR}/msas

## 6. Cleaning up
rm ${WORKDIR}/input/*
rm ${WORKDIR}/msas/*

## 5. Check metadata draw.
INSTANCE_NAME=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/name -H "Metadata-Flavor: Google")
ZONE=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/zone -H "Metadata-Flavor: Google" | awk -F/ '{print $NF}')
gcloud compute instances stop "$INSTANCE_NAME" --zone="$ZONE" --quiet

## 6. make the disk image and instance template
gcloud compute images create rb-msaflorafold-launchtemplate0326 \
    --project=${gcp_project_id} \
    --description=FloraFold\ \
MSA\ deployment\ disk\ with\ colabfold\ 1.6\ updated\ in\ 03/2026 \
    --source-disk=rb-msaflorafold-template \
    --source-disk-zone=us-west1-b \
    --storage-location=us

gcloud compute instance-templates create florafold-msa-template \
  --machine-type=n2-custom-28-196608 \
  --create-disk=auto-delete=yes,boot=yes,image=projects/${gcp_project_id}/global/images/rb-msaflorafold-launchtemplate0326,size=1200,type=pd-balanced \
  --scopes=https://www.googleapis.com/auth/cloud-platform 

## Testing deployment on SPOT
root_bucket='gs://rb-interactions/gmax-psojae/toy_data'
batch_id=1
run_name=testing
current_zone="us-west1-a"
gcloud compute instances create rb-msaflorafold-${run_name}-${batch_id} \
  --source-instance-template=florafold-msa-template \
  --zone=${current_zone} \
  --metadata=startup-script=/home/miles/startup.sh,batch_id=${batch_id},root_bucket=${root_bucket}
