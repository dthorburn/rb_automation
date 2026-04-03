# library imports
import argparse
import random
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

# Setup the Logger
class Logger(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self): # Required for Python 3 compatibility
        pass

# Initiate logger at the start of execution
sys.stdout = Logger("optimisation_log.txt")

"""CODON MAPPING AND FREQUECIES """

# IV2_Intron sequence 
IV2_INTRON = (
    "GTAAGTTTCTGCTTCTACCTTTGATATATATATAATAATTATCATTAATTAGTAGTAATATAATATTTCAAATATTTTTTTCAAAATAAAAGAATGTAGTATATAGCAATTGCTTTTCTGTAGTTTATAAGTGTGTATATTTTAATTTATAACTTTTCTAATATATGACCAAAATTTGTTGATGTGCAG"
)

# N. benthamiana weighted table 
CODON_DATA = {
    'A': [('GCT', 0.44), ('GCA', 0.31), ('GCC', 0.17), ('GCG', 0.09)],
    'C': [('TGT', 0.56), ('TGC', 0.44)],
    'D': [('GAT', 0.70), ('GAC', 0.30)],
    'E': [('GAA', 0.53), ('GAG', 0.47)],
    'F': [('TTT', 0.57), ('TTC', 0.43)],
    'G': [('GGT', 0.35), ('GGA', 0.32), ('GGC', 0.16), ('GGG', 0.16)],
    'H': [('CAT', 0.61), ('CAC', 0.39)],
    'I': [('ATT', 0.51), ('ATC', 0.26), ('ATA', 0.23)],
    'K': [('AAG', 0.57), ('AAA', 0.43)],
    'L': [('CTT', 0.26), ('TTG', 0.25), ('TTA', 0.13), ('CTC', 0.13), ('CTG', 0.12), ('CTA', 0.10)],
    'M': [('ATG', 1.0)],
    'N': [('AAT', 0.63), ('AAC', 0.37)],
    'P': [('CCT', 0.39), ('CCA', 0.35), ('CCC', 0.13), ('CCG', 0.13)],
    'Q': [('CAA', 0.51), ('CAG', 0.49)],
    'R': [('AGA', 0.30), ('AGG', 0.25), ('CGT', 0.15), ('CGA', 0.11), ('CGG', 0.10), ('CGC', 0.08)],
    'S': [('TCT', 0.28), ('TCA', 0.21), ('AGT', 0.18), ('TCC', 0.13), ('AGC', 0.13), ('TCG', 0.07)],
    'T': [('ACT', 0.36), ('ACA', 0.32), ('ACC', 0.21), ('ACG', 0.11)],
    'V': [('GTT', 0.42), ('GTG', 0.25), ('GTC', 0.17), ('GTA', 0.16)],
    'W': [('TGG', 1.0)],
    'Y': [('TAT', 0.55), ('TAC', 0.45)],
    'X': [('TGA', 0.39), ('TAA', 0.32), ('TAG', 0.29)]
}

# Pre-calculate CAI relative adaptiveness (w) for each codon (This early calculation improves script efficiency)
CAI_WEIGHTS = {}
for aa, codons in CODON_DATA.items():
    max_val = max(f for c, f in codons)
    for codon, freq in codons:
        CAI_WEIGHTS[codon] = freq / max_val

# Simplified map synonym for cut site repair
CODON_MAP = {c: aa for aa, items in CODON_DATA.items() for c, f in items}




# Stochastic sampling is used for the overall sequence to maintain the natural variety and GC balance 
def get_weighted_codon(aa):
    codons, weights = zip(*CODON_DATA[aa])
    return random.choices(codons, weights=weights)[0]

# Repair is used specifically to find BSAI cut sites (or any other cut sites) and remove them after stochastic rerolling of a sequence
def repair_bsai(dna_seq, args):
    sites = args.exclude.split(',')
    seq_list = list(dna_seq)
    
    
    for site in sites:
        while site in "".join(seq_list):
            current_dna = "".join(seq_list)
            idx = current_dna.find(site)
            # Find the start of the codon that overlaps the start of the site
            start_codon_idx = (idx // 3) * 3
            
            applied_fix = False
            # Try mutating 3 consecutive codons to break the 6bp site
            for offset in [0, 3, 6]:
                pos = start_codon_idx + offset
                if pos + 3 > len(current_dna): break
                
                old_codon = current_dna[pos:pos+3]
                target_aa = CODON_MAP.get(old_codon)
                # Find synonyms that aren't the current codon
                synonyms = [c for c, aa in CODON_MAP.items() if aa == target_aa and c != old_codon]
                
                if synonyms:
                    # Pick the first closest synonym to break the site
                    seq_list[pos:pos+3] = list(synonyms[0])
                    applied_fix = True
                    break 
            if not applied_fix: break # In that case the sequence cannot be modified
    return "".join(seq_list)

def check_other_constraints(dna_seq, args):
    """Stochastic Filter: Checks Homopolymers and GC Floor."""
    # Check for homopolymers in the sequence
    for nuc in ['A', 'T', 'G', 'C']:
        if nuc * args.homopolymer in dna_seq:
            return f"Homopolymer ({nuc * args.homopolymer})"
    # GC Window to check the overal count of GC base pairs in the sequence 
    # Use a dummy replacement to ensure window check doesn't count the Intron
    temp_seq = dna_seq.replace(IV2_INTRON, "N" * len(IV2_INTRON))
    for i in range(len(temp_seq) - args.window + 1):
        window = temp_seq[i:i + args.window]
        if "N" in window: continue 
        gc = (window.count('G') + window.count('C')) / len(window) * 100
        if gc < args.gc_floor:
            return f"Low GC ({gc:.1f}%)"
    return "OK"


"""INTRON INSERTION"""


def insert_ime(dna_seq):
    # We look for 'AG' ending a codon and 'G' starting the next one within first 60bp.
    # To explain the logic, Vancanneyt et al. 1990 / ST-LS1 paper describes that AG|G insertion site is the most optimal for the IV2.
    #  Isertion logic: AG | GT(intron)AG | G
    for i in range(3, 60, 3): 
        if dna_seq[i-2:i] == "AG" and dna_seq[i] == "G":
            return dna_seq[:i] + IV2_INTRON + dna_seq[i:]
    # If no optimal site found, insert after Start Codon
    return dna_seq[:3] + IV2_INTRON + dna_seq[3:]

# Combine all of the previous functions and add methionine in front of the sequences if not added. 
def optimise(aa_seq, args):
    if not aa_seq.startswith('M'): aa_seq = 'M' + aa_seq

    fail_report = {}
    last_dna_attempt = None

    for _ in range(2000): 
        # Phase 1: Stochastic Draft
        dna = "".join([get_weighted_codon(aa) for aa in aa_seq])
        
        # Phase 2: Algorithmic Repair (BsaI)
        dna = repair_bsai(dna, args)
        
        # Phase 3: Stochastic Validation (GC/Homopolymers)
        result = check_other_constraints(dna, args)
        
        if result == "OK":
            if args.use_ime:
                dna = insert_ime(dna)
                dna = repair_bsai(dna, args) # Final check for junctions
            return dna, "Success"
        
        fail_report[result] = fail_report.get(result, 0) + 1
        last_dna_attempt = dna

    # Fallback for "Best Effort"
    final_dna = last_dna_attempt
    if args.use_ime and final_dna:
        final_dna = insert_ime(final_dna)
        final_dna = repair_bsai(final_dna, args)
            
    return final_dna, fail_report

#  MAIN COMMANDS
def main():
    parser = argparse.ArgumentParser(description="Benthi Optimiser")
    parser.add_argument("fasta", help="Input protein FASTA file")
    parser.add_argument("--window", type=int, default=50, help="Specify the size of a window, default=50")
    parser.add_argument("--gc_floor", type=float, default=35.0, help="Specify the GC-count floor percentage in the sequence, default=35")
    parser.add_argument("--exclude", type=str, default="GGTCTC,GAGACC", help="Specify the cut site sequence you want to avoid with a comma separator. By default BSAI cut site is always on. Example BSAI = GGTCTC,GAGACC")
    parser.add_argument("--homopolymer", type=int, default=6, help="Specify the maximum size of a homopolymers, default=6")
    parser.add_argument("--use_ime", action="store_true", help="Turn on/Turn off the IME implementation of IV2 in the sequence")
    parser.add_argument("--save_f", action="store_true", help="Save sequences that fail constraints but match translation")
    
    args = parser.parse_args()
    results = []
    failed_proteins = [] 

    for record in SeqIO.parse(args.fasta, "fasta"):
        print(f"Processing {record.id}...")
        aa = str(record.seq).upper().replace('*', '')
        expected_aa = aa if aa.startswith('M') else 'M' + aa
        
        dna, report = optimise(aa, args)

        if dna:
            # Sanity check: Remove intron ONLY for translation verification
            check_dna = dna.replace(IV2_INTRON, "") if args.use_ime else dna
            translated = str(Seq(check_dna).translate()).rstrip('*')
            
            if translated == expected_aa:
                # GC content calculation check 
                gc_content = (dna.count('G') + dna.count('C')) / len(dna) * 100
                
                # CAI on coding portion
                codons = [check_dna[i:i+3] for i in range(0, len(check_dna), 3)]
                log_w_sum = sum(math.log(CAI_WEIGHTS.get(c, 0.01)) for c in codons)
                cai_score = math.exp(log_w_sum / len(codons))
                
                # Verify length
                intron_detected = "YES (+189bp)" if IV2_INTRON in dna else "NO"

                if report == "Success":
                    # Succesful result when DNA matches the original protein sequence when translated, and passes through all of the constraints
                    results.append(SeqRecord(Seq(dna), id=f"{record.id}", description=""))
                    print(f"[OK] Intron: {intron_detected} | GC: {gc_content:.1f}% | CAI: {cai_score:.3f} | Length: {len(dna)}bp")
                else:
                    # The result in case the DNA matches the original protein sequence when translated, but failed the constraints optimisation
                    save_msg = "RESULT SAVED" if args.save_f else "NOT SAVED"
                    print(f"[FAILED] Intron: {intron_detected} | Failures: {report} | CAI: {cai_score:.3f} | Length: {len(dna)}bp | {save_msg}")
                    failed_proteins.append((record.id, report))
                    
                    if args.save_f:
                        results.append(SeqRecord(Seq(dna), id=f"{record.id}", description=""))
            else:
                # Newely optimised DNA does not match the original protein sequence when translated back 
                print(f"[ERROR] Translation mismatch for {record.id}")
                failed_proteins.append((record.id, "Translation Mismatch"))
        else:
            # Hard failure in case the DNA could not be generated at all, basically script failed 
            print(f"[HARD ERROR] {report}")
            failed_proteins.append((record.id, f"Failed: {report}"))

    #Save the results
    output_name = f"{args.fasta.split('.')[0]}_OPT.fasta"
    if results:
        SeqIO.write(results, output_name, "fasta")

    # Output the summary table for failed proteins
    if failed_proteins:
        print("\n" + "="*50)
        print(f"{'Protein ID':<25} | {'Issue/Report'}")
        print("-"*50)
        for pid, issue in failed_proteins:
            print(f"{pid:<25} | {issue}")
        print("="*50)

    print(f"\nOptimisation complete. {len(results)} sequences saved to {output_name}.")  

if __name__ == "__main__":
    main()