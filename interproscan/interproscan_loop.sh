#!/bin/bash
export PATH=/home/resurrect/conda/envs/iprscan/bin:$PATH

IPR_DIR=/mnt/sda1/interproscan/interproscan-5.71-102.0
DATA_DIR=/mnt/sdb1/candidate_search/crop_nlrs/gmax
CPUS=16

cd ${IPR_DIR}

fastas=`ls -1 ${DATA_DIR}/*/*pep.fasta`
#sed -i "s/\*//g" ${fasta}
for fasta in $fastas
do
	echo "Starting ${fasta}: `date`"
	outfile=`basename $fasta | sed -e "s/.fasta/_interproscan.gff/"`
	echo $outfile

	./interproscan.sh -i ${fasta} \
  		-f gff3 -t "p" -o ${DATA_DIR}/${outfile} \
		-cpu ${CPUS} -appl Pfam,Gene3D,SUPERFAMILY,PRINTS,SMART,CDD,ProSiteProfiles -dp
done
