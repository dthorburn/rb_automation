export PATH="/home/resurrect/resurrectbio/01_projects/03_colabfold/01_installing/localcolabfold152:$PATH"

case $1 in
 -[h?] | --help)
        printf "Usage:\n\t1. Launch GNU Screen using \`screen\`.\n\t2. Update paramaters if necessary. Otherwise, any fasta in /home/miles/input will be used. \n\t3. Submit jobs using \`bash run_batch_msa.sh &> run.log &\`. \n\t4. Press ctrl+a+d to detach from screen. If terminate_vm option is set to 1, the VM will automatically terminate on completion.\n"
        exit 0;;
esac

## Paramaters
export_to_bucket=1
terminate_vm_on_finish=0
threads=28 ## Limit is 32 for workstation1 CPU

## 01_colabfold_uniref30  02_colabfold_uniref30_colabfoldenv  03_colabfold_florafold50  04_colabfold_florafold70  05_colabfold_florafold90
## 06_colabfold_florafold50_uniref30  07_colabfold_florafold70_uniref30  08_colabfold_florafold90_uniref30
input="/mnt/sdb1/florafold/04_screens/asr/working_nlrs/fasta/Cand_working_nlrs.fasta"
outdir="/mnt/sdb1/florafold/04_screens/asr/working_nlrs/msas"
db_path="/mnt/sda1"
floradb_path="/mnt/sdb1/florafold/02_clusters/reduced/04_profiledbs_ps1"

cloud_storage="aws" ## AWS | GCP
aws_bucket="s3://rb-kagome/predictions/"
gcp_bucket="gs://"

## Unsure why, but it requires you to also put the parent directory of the databases, even if I've used the db1 db3 options.
echo "~~Starting MSAs: `date`"
#  --db1 ${db_path}/uniref30_2202_db \
#  --db3 ${db_path}/colabfold_envdb_202108_db \
colabfold_search \
  --db1 ${db_path}/uniref30_2202_db \
  --db3 ${floradb_path}/florafolddb_reduced50 \
  --threads ${threads} \
  --db-load-mode 2 \
  --use-env 1 \
  ${input} ${db_path} ${outdir} &&
  echo "~~Finished with FloraFold50 MSAs: `date`"

## Renaming MSAs -- NOT NEEDED ANYMORE. THEIR UPDATE FIXED THIS.
echo "~~Renaming MSAs: `date`"
for a3m in ${outdir}/*.a3m
do
  echo $a3m
  ## Hope this works for all sequences. If a multimer, the first line is about sequence lengths, not sequence names.
  newname=`grep ">" $a3m | sed 's/\s.*$//' | head -n 1 | sed -e "s/>//"`
  echo $newname
  mv $a3m ${outdir}/${newname}_florafold.a3m
done

if [ ${export_to_bucket} == 1 ]
then
  echo "~~Exporting MSAs: `date`"
  if [ ${cloud_storage} == "aws" ]
  then
    aws s3 cp --recursive ${outdir}/ ${aws_bucket} --include "*a3m"
    aws s3 cp --recursive *log ${aws_bucket}
  else
    gsutil -m cp ${outdir}/* ${gcp_bucket}/
    gsutil -m cp ~/*.log ${gcp_bucket}/
  fi
  #rm ${input}
fi

if [ ${terminate_vm_on_finish} == 1 ]
then
  echo "~~Shutting Down Worstation: `date`"
  sudo shutdown -h now
fi
