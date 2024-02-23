#!/usr/bin/env bash

## Order of arguments should be: nlr fasta, effector fasta, AA length cutoff, and output file name
Rscript ./fasta_combination.R 	/path/to/nlr.fasta \
								/path/to/effector.fasta \
								2100 \
								/path/to/output_file.fasta

