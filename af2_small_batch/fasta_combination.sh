#!/usr/bin/env bash

## Usage: fasta_combination.sh /path/nrl.fasta /path/effector.fasta [max length filter, ~2100] /path/output.fasta  

## Remember to use the correct environment. I set up conda af2_prep env for this. 
if echo $CONDA_PREFIX | grep -ql "af2_prep"
then
	echo "~~ Correct conda env loaded"
else
	echo "ERROR: Please load af2_prep env \`conda activate af2_prep\` and try again"
	exit 1
fi

## Order of arguments should be: nlr fasta, effector fasta, AA length cutoff, and output file name
Rscript ./fasta_combination.R 	${1} ${2} ${3} ${4}
