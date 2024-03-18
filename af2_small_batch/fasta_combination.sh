#!/usr/bin/env bash

## Remember to use the correct environment. I set up conda af2_prep env for this. 

## Order of arguments should be: nlr fasta, effector fasta, AA length cutoff, and output file name
Rscript ./fasta_combination.R 	/path/to/nlr.fasta \
								/path/to/effector.fasta \
								2100 \
								/path/to/output_file.fasta

