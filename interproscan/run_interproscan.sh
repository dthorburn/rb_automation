export PATH=/home/resurrect/conda/envs/iprscan/bin:$PATH

IPR_DIR=/mnt/sda1/interproscan/interproscan-5.71-102.0
DATA_DIR=/mnt/sdb1/candidate_search/crop_nlrs/osativa/nippobare

fasta=`ls -1 ${DATA_DIR}/*.fasta`
sed -i "s/\*//g" ${fasta}
outfile=`basename $fasta | sed -e "s/.fasta/_interproscan.gff/"`

cd ${IPR_DIR}
./interproscan.sh -i ${fasta} \
  -f gff3 -t "p" -o ${DATA_DIR}/${outfile} \
  -cpu 32 -appl Pfam,Gene3D,SUPERFAMILY,PRINTS,SMART,CDD,ProSiteProfiles -dp

