## Steps to install localcolabfold on a VM instance for template.
## 1. Check for install requirements
git --version
nvcc --version
nvidia-smi
curl --version
wget --version
gcc --version
conda config --set auto_activate_base false
bash

## 2. Make directories
mkdir -p msas/done
mkdir -p predictions/done

## 3. Install pixi & localcolabfold
curl -fsSL https://pixi.sh/install.sh | sh
git clone https://github.com/yoshitakamo/localcolabfold.git
cd localcolabfold
pixi install && pixi run setup
pixi add gcc

## 4. Create prediction runner and watcher scripts.
nano run_predictions.sh
nano watcher.sh
nano startup.sh
sudo chmod +x run_predictions.sh
sudo chmod +x watcher.sh
sudo chmod +x startup.sh

## 5. Get toy data and test
export PATH="/home/miles/localcolabfold/.pixi/envs/default/bin/:${PATH}"
export LD_LIBRARY_PATH=$HOME/localcolabfold/.pixi/envs/default/lib:$LD_LIBRARY_PATH
cd ~/
gsutil -m ls gs://rb-interactions/gmax-asr/CcRPP1/toy_dataset/msas/batch1/"*"a3m | head -n 1 | gsutil -m cp -I msas
colabfold_batch \
    --num-recycle 2 \
    --num-models 1 \
    ./msas ./predictions

ROOT_BUCKET="gs://rb-interactions/gmax-asr/CcRPP1/toy_dataset"
BATCH_ID=1
nohup sudo -u miles bash /home/miles/watcher.sh ${BATCH_ID} ${ROOT_BUCKET} &
pkill -9 -f watcher.sh

## 6. Cleaning up
unset LD_LIBRARY_PATH
rm msas/done/*
rm predictions/*
rm predictions/done/*
rm nohop.out
rm RB_F*

## 5. Check metadata draw.
INSTANCE_NAME=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/name -H "Metadata-Flavor: Google")
ZONE=$(curl -s http://metadata.google.internal/computeMetadata/v1/instance/zone -H "Metadata-Flavor: Google" | awk -F/ '{print $NF}')
gcloud compute instances stop "$INSTANCE_NAME" --zone="$ZONE" --quiet

## 6. make the disk image and instance template
gcloud compute images create rb-florafold-launchtemplate0526 \
    --project=${gcp_project_id} \
    --description=FloraFold\ \
deployment\ disk\ with\ colabfold\ 1.6\ updated\ in\ 05/2026 \
    --source-disk=rb-florafoldpredict-template \
    --source-disk-zone=us-east4-c \
    --storage-location=us
gcloud compute instance-templates create florafold-l4-template-v3 \
  --machine-type=g2-custom-8-49152 \
  --accelerator=type=nvidia-l4,count=1 \
  --create-disk=auto-delete=yes,boot=yes,image=projects/${gcp_project_id}/global/images/rb-florafold-launchtemplate0526,size=150,type=pd-ssd \
  --scopes=https://www.googleapis.com/auth/cloud-platform \
  --maintenance-policy=TERMINATE 

## Testing deployment on SPOT
root_bucket='gs://rb-interactions/gmax-asr/CcRPP1/toy_dataset'
batch_id=2
models=1
recycles=1
run_name=testing
current_zone="us-west1-a"
gcloud compute instances create rb-gpucfold-l4-${run_name}-${batch_id} \
  --source-instance-template=florafold-l4-template \
  --zone=${current_zone} \
  --metadata=startup-script=/home/miles/startup.sh,batch_id=${batch_id},models=${models},recycles=${recycles},root_bucket=${root_bucket}

#### Additional steps added after first running
## TO make it useable by other projects run the following
gcp_service_account=""
gcp_project_id=""
gcp_image_project_id=""

gcloud compute images add-iam-policy-binding rb-florafold-launchtemplate0526 \
    --project="${gcp_image_project_id}" \
    --member="serviceAccount:${gcp_service_account}" \
    --role="roles/compute.imageUser"

## The GCP project ID here is for the new project

gcloud compute instance-templates create florafold-l4-template2 \
  --machine-type=g2-custom-8-49152 \
  --project=${gcp_project_id} \
  --accelerator=type=nvidia-l4,count=1 \
  --create-disk=auto-delete=yes,boot=yes,image=projects/${gcp_image_project_id}/global/images/rb-florafold-launchtemplate0526,size=150,type=pd-ssd \
  --scopes=https://www.googleapis.com/auth/cloud-platform \
  --maintenance-policy=TERMINATE 
