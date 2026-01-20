## Use the af2_ss conda env for this. It's running python 3.10 and it clearly outdated, but this script needs rewriting otherwise.  
import os
import glob
import argparse
import pandas as pd
from tqdm import tqdm
import MDAnalysis as mda
from Bio.PDB import PDBParser, DSSP
from MDAnalysis.analysis import distances

## Functions
def ss3_ss8_conv(ss8_seq):
    ss_map = {
        'H': 'H', 'G': 'H', 'I': 'H',
        'E': 'E', 'B': 'E', 
        'T': 'C', 'S': 'C', 'C': 'C', 'L': 'C', 'P': 'C', 
        '-': '-'
    }
    ## For each ss8 in ss8_seq, it applies the ss_map conversion. - Getting more used to this logic.
    ss3_seq = ''.join([ss_map[ss8] for ss8 in ss8_seq])
    return(ss3_seq)

def find_disulphide(input_pdb, min_distance=1.8, max_distance=2.2):
    disulphides = 0
    temp_uni = mda.Universe(input_pdb)
    cyss = temp_uni.select_atoms("resname CYS")
    sulphur_ats = cyss.select_atoms("name SG")
    ## Seems a little slow on larger molecules. Need a more efficient solution? 
    for i in range(len(sulphur_ats)):
        for j in range(i + 1, len(sulphur_ats)):
            distance = distances.distance_array(sulphur_ats[i].position, sulphur_ats[j].position)
            if min_distance < distance < max_distance:
                disulphides = disulphides + 1
    return(disulphides)

def calculate_2d_structure(sname, input_pdb, is_eff, min_distance, max_distance):
    """
    Calculate the secondary structure of PDB files.
    Args:
      sname: Simple file name.
      input_pdb: A PDB file generated from AlphaFold2.
    Returns
      ss_out: A string of secondary struture information. 
    """
    p = PDBParser()
    temp_struc = p.get_structure(sname, input_pdb)
    model = temp_struc[0]
    dssp = DSSP(model, input_pdb)
    temp_ss = [residue[2] for residue in dssp]

    ss_out = {
        "file": sname,
        "ss8": ''.join(temp_ss),
        "ss3": ss3_ss8_conv(''.join(temp_ss))
    }
    
    ## Boolean handling of effector disulphide estimation
    if is_eff:
        ss_out["ds_bonds"] = find_disulphide(input_pdb, min_distance, max_distance)
        ss_out["prop_unstruc"] = round(ss_out["ss3"].count('-') / len(ss_out["ss3"]), 4)
    if ss_out:
        res = pd.DataFrame([ss_out])
        return(res)

parser = argparse.ArgumentParser(description="Define secondary structure of input pdbs")
# Positional arguments for the PDB directory and output file
parser.add_argument("pdb_dir", type=str, help="path to directory with monomer PDB files")
parser.add_argument("output_file", type=str, help="path for output CSV")
# Optional arguments for molecule type and minimum and maximum distance thresholds
parser.add_argument("--effector", action = 'store_true', help="effector flag to estimate disulphide bonds in molecule")
parser.add_argument("--min_distance", type=float, default=1.8, 
                    help="minimum distance threshold in angstroms (default: 1.8 Å)")
parser.add_argument("--max_distance", type=float, default=2.2, 
                    help="maximum distance threshold in angstroms (default: 2.2 Å)")
args = parser.parse_args()

pdb_dir  = os.path.abspath(args.pdb_dir)
output_file = os.path.abspath(args.output_file)
is_eff = args.effector
min_distance = args.min_distance
max_distance = args.max_distance

pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
total_files = len(pdb_files)
all_ss = []

## Testing
#pdb_dir = os.path.abspath("/home/dthorbur/Resurrect_Bio/Scripts/af2_processing/data/test_eff_pdbs/")
#input_pdb = pdb_files[0]

with tqdm(total=total_files, desc="Processing PDBs") as pbar:
    for temp_pdb in pdb_files:
        sname = os.path.basename(temp_pdb)
        pbar.set_description(f"Processing PDB: {sname}")
        temp_ss = calculate_2d_structure(sname, temp_pdb, is_eff, min_distance, max_distance)
        all_ss.append(temp_ss)
        pbar.update(1)

if all_ss:
    combined_df = pd.concat(all_ss, ignore_index=True)
    print(combined_df)
    combined_df.to_csv(output_file, index=False)
else:
    print("No PDB files found.")