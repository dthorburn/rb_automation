#!/usr/bin/env bash
## Probably overkill.
export PATH="/home/miles/localcolabfold/colabfold-conda/bin:/home/miles/localcolabfold/conda/bin/:$PATH"

case $1 in
 -[h?] | --help)
        printf "Usage:\n\t1. Launch GNU Screen using \`screen\`.\n\t2. Update paramaters if necessary. Otherwise, any fasta in /home/miles/input will be used. \n\t3. Submit jobs using \`bash run_batch_msa.sh &> run.log &\`. \n\t4. Press ctrl+a+d to detach from screen. If terminate_vm option is set to 1, the VM will automatically terminate on completion.\n"
        exit 0;;
esac

## Paramaters
export_to_bucket=1
terminate_vm_on_finish=1
threads=16

input="/home/miles/input/*.fasta"
outdir="/home/miles/msas/"
db_path="/home/miles/databases/"
gcp_bucket="gs://gmax/HgEffectors-longread/msas"
batch_name="batch1"

## Unsure why, but it requires you to also put the parent directory of the databases, even if I've used the db1 db3 options.
echo "~~Starting MSAs: `date`"
colabfold_search \
  --db1 /home/miles/databases/uniref30_2302_db \
  --db3 /home/miles/databases/colabfold_envdb_202108_db \
  --threads ${threads} \
  --db-load-mode 2 \
  --use-env 1 \
  ${input} ${db_path} ${outdir} &&
  echo "~~Finished with MSAs: `date`"

## Renaming MSAs
echo "~~Renaming MSAs: `date`"
for a3m in ${outdir}/*.a3m
do
  echo $a3m
  ## Hope this works for all sequences. If a multimer, the first line is about sequence lengths, not sequence names. 
  newname=`grep ">" $a3m | sed 's/\s.*$//' | head -n 1 | sed -e "s/>//"`
  echo $newname
  mv $a3m ${outdir}/${newname}.a3m
done

if [ ${export_to_bucket} == 1 ]
then
  echo "~~Exporting MSAs: `date`"
  gsutil -m mv ${outdir}/* ${gcp_bucket}/${batch_name}/
  gsutil -m mv ~/*.log ${gcp_bucket}/
  rm ${input}
fi

if [ ${terminate_vm_on_finish} == 1 ]
then
  echo "~~Terminating VM: `date`"
  sudo shutdown -h now
fi