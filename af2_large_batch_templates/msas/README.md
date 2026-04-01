# RB Bulk GCP MSA FloraFold50 v1 Running
These scripts are designed only for bulk MSA generation when the workstations would be time prohibitive. They have been recently rebuilt (03/26) to be considerably earier to deploy and should be a little faster than previous solutions. 

## Requirements
1. Access to GCP with the ability to deploy scripts and access to the RB Main GCP project and any specific JDA or other GCP project for budgeting reasons. 
2. A set of fasta files with up to 3,000 sequences is a good start. Please use `seqkit` to split the fasta files into batches.
3. Fasta files need to be regex readable with this pattern `*0${BATCH_ID}.fasta` where `${BATCH_ID}` is a number that goes sequentially up from 1. 
4. 

## Usage


# Debugging
1. If jobs continue to crash or improvements need to be implemented on a new image, remove `startup-script=/home/miles/startup.sh` from the metadata command. This _should_ stop the autodeletion, but will require manual metadata parsing to start the job. 

Author: Miles Thorburn <miles@resurrect.bio>