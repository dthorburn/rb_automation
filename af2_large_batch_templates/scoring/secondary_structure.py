import os
import re
import glob
import argparse
import numpy as np
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
    ## For each ss8 in ss8_seq, it applies the ss_map conversion.
    ss3_seq = ''.join([ss_map[ss8] for ss8 in ss8_seq])
    return(ss3_seq)

def count_disulphides(sulphur_atoms, min_distance, max_distance):
    if len(sulphur_atoms) < 2:
        return 0
    positions = sulphur_atoms.positions
    dist_matrix = distances.distance_array(positions, positions)
    ## Only looks at upper triangle to avoid counting pairs twice
    upper = np.triu(dist_matrix, k=1)
    return np.sum((upper > min_distance) & (upper < max_distance))

def calculate_2d_structure(sname, rname, input_pdb, min_distance, max_distance):
    """
    Calculate the secondary structure of PDB files.
    Args:
      sname: Simple file name.
      input_pdb: A PDB file with 2 chains generated from AlphaFold2.
    Returns
      ss_out: A string of secondary struture information. 
    """
    p = PDBParser()
    temp_struc = p.get_structure(sname, input_pdb)
    model = temp_struc[0]
    dssp = DSSP(model, input_pdb)
    dsspA = {k: dssp[k] for k in dssp.keys() if k[0] == "A"}
    dsspB = {k: dssp[k] for k in dssp.keys() if k[0] == "B"}
    temp_ss_A = [dsspA[k][2] for k in dsspA]
    temp_ss_B = [dsspB[k][2] for k in dsspB]
    
    temp_uni = mda.Universe(input_pdb)
    disulphides_A = count_disulphides(temp_uni.select_atoms("segid A and not name H* and resname CYS and name SG"), min_distance, max_distance)
    disulphides_B = count_disulphides(temp_uni.select_atoms("segid B and not name H* and resname CYS and name SG"), min_distance, max_distance)
    
    ss_out = {
        "complex": sname,
        "rank": rname,
        "ss8_chainA": ''.join(temp_ss_A),
        "ss3_chainA": ss3_ss8_conv(''.join(temp_ss_A)),
        "ss8_chainB": ''.join(temp_ss_B),
        "ss3_chainB": ss3_ss8_conv(''.join(temp_ss_B)),
        "ds_bonds_chainA": disulphides_A,
        "ds_bonds_chainB": disulphides_B
    }
    ss_out["prop_unstruc_A"] = round(ss_out["ss3_chainA"].count('-') / len(ss_out["ss3_chainA"]), 4)
    ss_out["prop_unstruc_B"] = round(ss_out["ss3_chainB"].count('-') / len(ss_out["ss3_chainB"]), 4)
    
    if ss_out:
        res = pd.DataFrame([ss_out])
        return(res)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Define secondary structure of input pdbs")
    # Positional arguments for the PDB directory and output file
    parser.add_argument("pdb_dir", type=str, help="path to directory with monomer PDB files")
    parser.add_argument("output_file", type=str, help="path for output CSV")
    # Optional arguments for molecule type and minimum and maximum distance thresholds
    parser.add_argument("--min_distance", type=float, default=1.8, 
                        help="minimum distance threshold in angstroms (default: 1.8 Å)")
    parser.add_argument("--max_distance", type=float, default=2.3, 
                        help="maximum distance threshold in angstroms (default: 2.3 Å)")
    args = parser.parse_args()

    pdb_dir  = os.path.abspath(args.pdb_dir)
    output_file = os.path.abspath(args.output_file)
    min_distance = args.min_distance
    max_distance = args.max_distance
    ## Testing
    #pdb_dir="/mnt/c/Users/miles/Documents/Resurrect_Bio/Projects/09_commercialagreements/02_corn/05_interactions/scoring/files/"
    #min_distance=1.8
    #max_distance=2.2
    #input_pdb = pdb_files[0]

    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    total_files = len(pdb_files)
    all_ss = []

    with tqdm(total=total_files, desc="Processing PDBs") as pbar:
        for temp_pdb in pdb_files:
            fname = os.path.basename(temp_pdb)
            sname = re.sub("_unrelaxed.*$", "", fname)
            rname = re.sub(r"(^.*_unrelaxed_)(rank_[0-9]+)(_alphafold2.*$)", r"\2", fname)
            pbar.set_description(f"Processing PDB: {sname}")
            temp_ss = calculate_2d_structure(sname, rname, temp_pdb, min_distance, max_distance)
            all_ss.append(temp_ss)
            pbar.update(1)

    if all_ss:
        combined_df = pd.concat(all_ss, ignore_index=True)
        print(combined_df)
        combined_df.to_csv(output_file, index=False)
    else:
        print("No PDB files found.")
