gcloud compute instances list \
    --filter="status=RUNNING" \
    --format="csv(name,zone,no-heading)" \
    --project "${gcp_project_id}" > vm_list.txt

tail -n +2 vm_list.txt > vm_list2.txt

while IFS=',' read -r INSTANCE ZONE; do
    ZONE="${ZONE%%,*}"
    echo "===== $INSTANCE ($ZONE) ====="
    gcloud compute ssh "$INSTANCE" \
        --zone "$ZONE" \
        --project "${gcp_project_id}" \
        --command "top -b -n1 | grep 'colabf'" < /dev/null
done < vm_list2.txt > vm_activity_5pm.log

#        --command "top -b -n1 | head -n10" < /dev/null


while IFS=',' read -r INSTANCE ZONE BATCH; do
    echo "===== ${INSTANCE} (${ZONE}; ${BATCH}) ====="
    gcloud compute ssh "$INSTANCE" \
        --zone "$ZONE" \
        --project "${gcp_project_id}" \
        --command "gsutil -m cp /home/miles/*.log gs://rb-corn-tar-spot-mbsp/03_interactions/01_zmays_tarspot/predictions" < /dev/null
done < to_stop

while IFS=',' read -r INSTANCE ZONE BATCH; do
    echo "===== ${INSTANCE} (${ZONE}; ${BATCH}) ====="
    gcloud compute ssh "$INSTANCE" \
        --zone "$ZONE" \
        --project "${gcp_project_id}" \
        --command "gsutil -m cp /home/miles/predictions/*pdb gs://rb-corn-tar-spot-mbsp/03_interactions/01_zmays_tarspot/predictions/batch${BATCH}; gsutil -m cp /home/miles/predictions/*v1.json gs://rb-corn-tar-spot-mbsp/03_interactions/01_zmays_tarspot/predictions/batch${BATCH}; gsutil -m cp /home/miles/predictions/*0.json gs://rb-corn-tar-spot-mbsp/03_interactions/01_zmays_tarspot/predictions/batch${BATCH}; sudo poweroff" < /dev/null 
done < to_stop

while IFS=',' read -r INSTANCE ZONE BATCH; do
    echo "===== ${INSTANCE} (${ZONE}; ${BATCH}) ====="
    gcloud compute ssh "$INSTANCE" \
        --zone "$ZONE" \
        --project "${gcp_project_id}" \
        --command "sudo poweroff" < /dev/null 
done < to_stop