## 0. Launch a VM

gcp_service_account=""
gcp_project_id=""
gcloud compute instances create rb-scoring-template-v1 \
    --project=${gcp_project_id} \
    --zone=us-west1-b \
    --machine-type=e2-standard-8 \
    --network-interface=network-tier=PREMIUM,stack-type=IPV4_ONLY,subnet=default \
    --metadata=enable-osconfig=TRUE \
    --maintenance-policy=MIGRATE \
    --provisioning-model=STANDARD \
    --service-account=${gcp_service_account} \
    --scopes=https://www.googleapis.com/auth/cloud-platform \
    --create-disk=auto-delete=yes,boot=yes,device-name=rb-scoring-template-v1,image=projects/debian-cloud/global/images/debian-12-bookworm-v20260513,mode=rw,size=300,type=pd-balanced \
    --no-shielded-secure-boot \
    --shielded-vtpm \
    --shielded-integrity-monitoring \
    --labels=goog-ops-agent-policy=v2-template-1-7-0,goog-ec-src=vm_add-gcloud \
    --reservation-affinity=any
  
## 1. Check for install requirements
sudo apt-get update && sudo apt-get install -y \
  git \
  gcc \
  g++ \
  make \
  wget \
  curl 

## 2. Install conda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda config --set auto_activate_base false
source ~/.bashrc
conda init bash
bash
rm Miniconda3-latest-Linux-x86_64.sh 

## 3. Make directories
mkdir -p predictions
mkdir -p scored
mkdir -p completed

## 4. create conda env - might need updating if this branch is merged into main
curl -O https://raw.githubusercontent.com/dthorburn/rb_automation/mt-dev/af2_large_batch_templates/scoring/af2_scoring.yml
curl -O https://raw.githubusercontent.com/dthorburn/rb_automation/mt-dev/af2_large_batch_templates/scoring/scoring_interactions.py
curl -O https://raw.githubusercontent.com/dthorburn/rb_automation/mt-dev/af2_large_batch_templates/scoring/secondary_structure.py
curl -O https://raw.githubusercontent.com/dthorburn/rb_automation/mt-dev/af2_large_batch_templates/scoring/merge_results.py
curl -O https://raw.githubusercontent.com/dthorburn/rb_automation/mt-dev/af2_large_batch_templates/scoring/ipsae.py
curl -O https://raw.githubusercontent.com/dthorburn/rb_automation/mt-dev/af2_large_batch_templates/scoring/process_interactions.sh
curl -O https://raw.githubusercontent.com/dthorburn/rb_automation/mt-dev/af2_large_batch_templates/scoring/scoring_startup.sh
chmod 766 ./*py
chmod 766 ./*sh

conda create -n af2_scoring --file af2_scoring.yml

## Once everything is installed and ready, you can make an image and template
gcloud compute images create rb-ffscoring-image0526 \
    --project=${gcp_project_id} \
    --description=FloraFold\ \
scoring\ disk\ updated\ in\ 05/2026 \
    --source-disk=rb-scoring-template-v1 \
    --source-disk-zone=us-west1-b \
    --storage-location=us

gcloud compute instance-templates create rb-ffscoring-template-v1 \
  --machine-type=e2-standard-8  \
  --create-disk=auto-delete=yes,boot=yes,image=projects/${gcp_project_id}/global/images/rb-ffscoring-image0526,size=300,type=pd-ssd \
  --scopes=https://www.googleapis.com/auth/cloud-platform 

## Only run this if you want to change the project it can be deployed one
gcloud compute images add-iam-policy-binding rb-ffscoring-image0526 \
    --project="${gcp_image_project_id}" \
    --member="serviceAccount:${gcp_service_account}" \
    --role="roles/compute.imageUser"

gcloud compute instance-templates create rb-ffscoring-template-v1 \
  --machine-type=e2-standard-8  \
  --project=${gcp_project_id} \
  --create-disk=auto-delete=yes,boot=yes,image=projects/${gcp_image_project_id}/global/images/rb-ffscoring-image0526,size=300,type=pd-ssd \
  --scopes=https://www.googleapis.com/auth/cloud-platform 



## For deployment
## Include up the gs:// and all the way up to the directory with /msas and /predictions folders.
root_bucket=""
run_name=""
mina=3
maxa=5
paec=10
distc=10
gcloud compute instances create rb-ffscoring-${run_name} \
                --source-instance-template=projects/${gcp_project_id}/global/instanceTemplates/rb-ffscoring-template-v1 \
                --project=${gcp_project_id} \
                --service-account=${gcp_service_account} \
                --scopes=https://www.googleapis.com/auth/cloud-platform \
                --zone="us-west1-b" \
                --metadata=startup-script=/home/miles/scoring_startup.sh,root_bucket=${root_bucket},mina=${mina},maxa=${maxa},paec=${paec},distc=${distc}


