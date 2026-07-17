#!/bin/bash
## Realisticlally, this isn't enough time to move everything. Use the persistent file moving one instead.

preempted_checker() {
  while true; do
    if curl -s -H "Metadata-Flavor: Google" \
      http://metadata.google.internal/computeMetadata/v1/instance/preempted | grep TRUE; then
        echo "VM preempted! Uploading results..."
        gsutil -m cp -r /home/miles/predictions/* ${OUTPUT_BUCKET}
        shutdown -h now
    fi
    sleep 60
  done
}

preempted_checker &