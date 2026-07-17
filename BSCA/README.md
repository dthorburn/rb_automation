# BSAI Site Cut Automation (BSCA)
Unlike simple "search and destroy" scripts, this tool performs **context-aware synonymous recoding**:

1. **Scan:** Detects BsaI recognition sites (`GGTCTC` / `GAGACC`).
2. **Translate:** Maps the site to its corresponding Amino Acid codons.
3. **Recode:** Mutates the recognition site using synonymous codons (e.g., changing `GGT` to `GGC` for Glycine) to break the BsaI motif.
4. **Verify:** Re-scans and translates the new sequence to amino acids to ensure that newely created sequence matches the origial sequence.

# Usage 

Input: 
- Your Fasta file with nucleotide sequences

Output: 
- Fasta file with nucleotide sequences with BSAI sites removed (cleaned_dna.fasta)
- Fasta file with protein sequences for additional manual check 
(proteins.fasta)
```
Python3 BSCA.py /Path/to/Directory/DNA_NLRs.fasta
```

# Requirments

- biopython 