#!/usr/bin/env bash 

## Make sure you are in the right conda environment - 'conda activate pymolfree'
TOOL=/home/dthorbur/Resurrect_Bio/Scripts
DATA_DIR=$1
echo $DATA_DIR
mkdir ./montages

for file in $DATA_DIR/*.pdb
do
	echo $file
	newname=`basename $file | sed -e "s/_unrelaxed.*//"`
	python3 ${TOOL}/pymol_make_figure_monomer.py ${file} &&
		montage ${newname}*.png -tile 2x2 -geometry 900x900 ./montages/${newname}_all.png &&
		rm ${newname}*.png
done
