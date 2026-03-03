'''This is BSAIremoveScript. It's purpose to take our NLR nucleotide sequences and remove the BSAI cut sites 
without changing the final protein sequnce. 
The script operates via this logic: 
1. Create a dictionary so the the script knows for which amino acids the codons stand for.
2. Create a reading frame which reads 3 nucleobases recognising the amino acids as it moves through the sequence.
3. Recognise which codons does the BSAI cut site falls into, so that change in codon does not affect the amino acid translation. 
5. Remove stop codons at the end. 
6. Create the final files, which are 
    "cleaned_dna.fasta" - nucleotide sequences of your NLRs without BSAI cut sites. 
    "proteins.fasta" - amino acid sequence of your modified nucleotide sequences without BSAI cut sites. 
6. Double check the results by using biopython to create a translated temporary copy of amino acid sequence and compare with the proteins.fasta.

(the script heavily relies on Biopython libraries. Uncomment and run the first line if you don't have it installed)'''

#Uncomment this line to install biopython
#!python3 -m pip install biopython
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#CODING LIBRARY: A dictionary mapping DNA triplets to Amino Acid letters.
#'*' represents a Stop Codon (TAA, TAG, TGA).
CODON_MAP = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def get_synonym(old_codon):
    """Finds a different DNA triplet that codes for the same Amino Acid.
    This allows us to break restriction sites without changing the protein."""
    target_aa = CODON_MAP.get(old_codon)
    #Search for all codons that match the AA but are NOT the current one.
    options = [c for c, aa in CODON_MAP.items() if aa == target_aa and c != old_codon]
    #Return the first alternative found. If no alternative exists (Met or Trp), return original.
    return options[0] if options else old_codon

def delete_terminal_stop(seq_str):
    """Checks only the very last codon of the sequence. If it is a stop codon (*), it is removed to allow for C-terminal tagging."""
    #Extract the last 3 characters
    last_codon = seq_str[-3:]

    #Check if that specific triplet is a stop codon in our map
    if CODON_MAP.get(last_codon) == '*':
        return seq_str[:-3]
    
    #If it's a normal amino acid, return the sequence as is
    return seq_str

def BSA1rem_and_delete_stops(record):
    """First removes BsaI sites, then deletes stop codons."""
    #Standardise the sequence to uppercase and strip white spaces.
    seq_str = str(record.seq).upper().strip()
    
    #TARGET SITES: BsaI recognition sequences (Forward and Reverse Complement).
    sites = ["GGTCTC", "GAGACC"]
    
    for site in sites:
        while site in seq_str:
            idx = seq_str.find(site)
            #CODON ANCHORING: Find the start of the 3-base codon that contains the site's start.
            start_codon_idx = (idx // 3) * 3
            
            applied_fix = False
            #Try mutating the 1st, 2nd, or 3rd codon that overlaps with the 6bp site.
            for offset in [0, 3, 6]:
                current_pos = start_codon_idx + offset
                if current_pos + 3 > len(seq_str): break
                
                old_codon = seq_str[current_pos : current_pos+3]
                new_codon = get_synonym(old_codon)
                
                #Check 1: Is the new codon different? 
                #Check 2: Does it code for the same thing? (Double-check for safety)
                if new_codon != old_codon and CODON_MAP[old_codon] == CODON_MAP[new_codon]:
                    # Splice the new codon into the string.
                    seq_str = seq_str[:current_pos] + new_codon + seq_str[current_pos+3:]
                    applied_fix = True
                    break # Site is broken, exit the 'offset' loop.
            
            if not applied_fix:
                #If we reach here, the site is unfixable (likely purely Met/Trp/Stops).
                break 

    #Now that DNA is cleared of BSA1 sites, remove all stop codons entirely.
    final_seq = delete_terminal_stop(seq_str)
    return Seq(final_seq)

def run_process(filename):
    """Handles file input/output and verifies that the protein sequence wasn't corrupted."""
    #Checks if files exists to begin with
    if not os.path.exists(filename):
        print(f"File {filename} not found.")
        return

    #Use 'fasta-pearson' to deal with header formatting
    records = list(SeqIO.parse(filename, "fasta-pearson"))
    dna_results, aa_results = [], []

    for rec in records:
        #Process the sequence
        clean_dna = BSA1rem_and_delete_stops(rec)
        
        #VERIFICATION stage: Create a reference by deleting stops from original DNA.
        base_dna_deleted = delete_terminal_stop(str(rec.seq).upper().strip())
        
        #Translate both into amino acids (removing the '*' symbol for comparison).
        orig_pro = str(Seq(base_dna_deleted).translate(table=1)).replace('*', '')
        new_pro = str(clean_dna.translate(table=1)).replace('*', '')
        
        if orig_pro == new_pro:
            print(f"{rec.id}: Match. BsaI sites removed and Stop codons deleted.")
        else:
            print(f"{rec.id}: MISMATCH! Protein changed during BSAI site removal.")
            
        #Wrap the new sequence in SeqRecord objects for saving.
        dna_results.append(SeqRecord(clean_dna, id=rec.id, description="BsaI_Cleaned_Stops_Deleted"))
        aa_results.append(SeqRecord(clean_dna.translate(), id=rec.id, description="Protein_Output"))

    #Save outputs to new FASTA files.
    SeqIO.write(dna_results, "clean_dna.fasta", "fasta-2line")
    SeqIO.write(aa_results, "proteins.fasta", "fasta-2line")
    print("\nBSAI sites were removed")

"""Ensure your filename and path is correct. It should be nucleotide sequences with BSAI cut sites."""

run_process("C:/Users/BioData/sequences/YourFileName.fasta")