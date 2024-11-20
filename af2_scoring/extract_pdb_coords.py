import os
import re
import pdb
import glob
import argparse
import statistics
import numpy as np
import pandas as pd
import MDAnalysis as mda
import MDAnalysis.lib.distances
from collections import defaultdict

## Redundant now. Use the NP implementation. 
def find_residue_interactions(pdb_file, min_distance=1.5, max_distance=5.0, complex_name=None):
    ## If complex_name is not provided, generate it from the file name
    if complex_name is None:
        fname = os.path.basename(pdb_file)
        complex_name = re.sub("_unrelaxed.*$", "", fname)
    ## Load the structure using MDAnalysis
    interactions = []
    u = mda.Universe(pdb_file)
    chain_A = u.select_atoms("segid A")
    ## Iterate through each residue in chain A and find nearby residues
    for residue in chain_A.residues:
        ## Select atoms within max_distance from the current residue
        nearby_residues = u.select_atoms(f"around {max_distance} (not name H* and resid {residue.resid})")
        for nearby_residue in nearby_residues.residues:
            ## Only process interactions on different chains
            if residue.segid != nearby_residue.segid:
                ## Calculate the minimum distance between the residues
                distance = np.min(mda.lib.distances.distance_array(residue.atoms.positions, nearby_residue.atoms.positions))                
                ## Only add the interaction if the distance is within the specified range
                if min_distance <= distance <= max_distance:
                    interactions.append({
                        "complex": complex_name,
                        "residue1": f"{residue.resname}{residue.resid}",
                        "chain1": residue.segid,
                        "residue1_num": residue.resid,
                        "residue2": f"{nearby_residue.resname}{nearby_residue.resid}",
                        "chain2": nearby_residue.segid,
                        "residue2_num": nearby_residue.resid,
                        "distance": round(distance, 3)  # Rounded to 3 decimal places
                    })
    if interactions:
        df = pd.DataFrame(interactions)
        # This isn't the most efficient, but my other solutions didn't work.... 
        # df_subset = df[df["chain1"] == "A"]
    else:
        df = pd.DataFrame([{
            "complex": complex_name,
            "residue1": None,
            "chain1": None,
            "residue1_num": None,
            "residue2": None,
            "chain2": None,
            "residue2_num": None,
            "distance": None
        }])
    return df

def find_residue_interactions_np(pdb_file, min_distance=1.5, max_distance=5.0, complex_name=None):
    ## Generate complex name from file name if not provided
    if complex_name is None:
        fname = os.path.basename(pdb_file)
        complex_name = re.sub("_unrelaxed.*$", "", fname)
    ## Load the structure using MDAnalysis
    u = mda.Universe(pdb_file)
    ## Select non-hydrogen atoms
    chain_A = u.select_atoms("segid A and not name H*")
    chain_B = u.select_atoms("segid B and not name H*")
    ## Precompute atom positions for chain A and chain B
    A_positions = chain_A.positions
    B_positions = chain_B.positions
    ## Calculate all pairwise distances between atoms in chain A and chain B
    distances = mda.lib.distances.distance_array(A_positions, B_positions)
    ## Find atom pairs where the distance is within the specified range
    close_contacts = np.where((distances >= min_distance) & (distances <= max_distance))
    too_close = np.where(distances < min_distance)
    ## For pDockQ calculation
    contact_mask = (distances >= min_distance) & (distances <= max_distance)
    num_contacts = np.sum(contact_mask)
    ## Store interactions, but only one per residue pair
    interactions = []
    processed_pairs = set()  
    ## Iterate over close contacts (i.e., atom pairs within the distance threshold)
    ## zip(*tuple) pairs up the elements in the array objects. (i.e., output=(array([0, 1, 3]), array([4, 2, 5])) will become (0,4);(1,2);(3,5))
    for i, j in zip(*close_contacts):
        atom_A = chain_A[i]
        atom_B = chain_B[j]
        residue_A = atom_A.residue
        residue_B = atom_B.residue
        ## Check if this residue pair has been processed already
        if (residue_A.resid, residue_B.resid) not in processed_pairs:
            interactions.append({
                "complex": complex_name,
                "residue1": f"{residue_A.resname}{residue_A.resid}",
                "chain1": residue_A.segid,
                "residue1_num": residue_A.resid,
                "residue2": f"{residue_B.resname}{residue_B.resid}",
                "chain2": residue_B.segid,
                "residue2_num": residue_B.resid,
                "distance": round(distances[i, j], 3)  # Round to 3 decimal places
            })
            ## Mark this residue pair as processed
            processed_pairs.add((residue_A.resid, residue_B.resid))
    if interactions:
        df = pd.DataFrame(interactions)
    else:
        ## If no interactions, return an empty DataFrame with expected columns
        df = pd.DataFrame([{
            "complex": complex_name,
            "residue1": None,
            "chain1": None,
            "residue1_num": None,
            "residue2": None,
            "chain2": None,
            "residue2_num": None,
            "distance": None
        }])
    ## Calculating number of residues that are too close. 
    processed_too_close = set()
    for i, j in zip(*too_close):
        atom_A = chain_A[i]
        atom_B = chain_B[j]
        residue_A = atom_A.residue
        residue_B = atom_B.residue
        if (residue_A.resid, residue_B.resid) not in processed_too_close:
            processed_too_close.add((residue_A.resid, residue_B.resid))
    if processed_too_close:
        int_close = {
            "complex": complex_name,
            "close_atoms" : len(too_close[0]),
            "close_residues" : len(processed_too_close),
            "min_distance_threshold": min_distance
        }
        df2 = pd.DataFrame([int_close])
    else:
        df2 = pd.DataFrame([{
            "complex": complex_name,
            "close_atoms" : None,
            "close_residues" : None,
            "min_distance_threshold": min_distance
        }])
    return df, df2

## These functions to calculate updated pDockQ score are taken from where the paper indicates: https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py?ref_type=heads
## parse_atm_record and read_pdb are kept the same. calc_pdockq is updated.
def parse_atm_record(line):
    '''Get the atm record
    '''
    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['atm_alt'] = line[17]
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])
    return record

def read_pdb(pdbfile):
    chain_coords, chain_plddt = {}, {}
    with open(pdbfile, 'r') as file:
        for line in file:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            #Get CB - CA for GLY
            if record['atm_name']=='CB' or (record['atm_name']=='CA' and record['res_name']=='GLY'):
                if record['chain'] in [*chain_coords.keys()]:
                    chain_coords[record['chain']].append([record['x'],record['y'],record['z']])
                    chain_plddt[record['chain']].append(record['B'])
                else:
                    chain_coords[record['chain']] = [[record['x'],record['y'],record['z']]]
                    chain_plddt[record['chain']] = [record['B']]
    #Convert to arrays
    for chain in chain_coords:
        chain_coords[chain] = np.array(chain_coords[chain])
        chain_plddt[chain] = np.array(chain_plddt[chain])
    return chain_coords, chain_plddt

def calc_pdockq(chain_coords, chain_plddt, min_distance, max_distance):
    '''Calculate the pDockQ scores
    pdockQ = L / (1 + np.exp(-k*(x-x0)))+b
    L= 0.724 x0= 152.611 k= 0.052 and b= 0.018
    '''
    #Get coords and plddt per chain
    ch1, ch2 = [*chain_coords.keys()]
    coords1, coords2 = chain_coords[ch1], chain_coords[ch2]
    plddt1, plddt2 = chain_plddt[ch1], chain_plddt[ch2]
    #Calc 2-norm
    mat = np.append(coords1, coords2,axis=0)
    a_min_b = mat[:,np.newaxis,:] -mat[np.newaxis,:,:]
    dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
    l1 = len(coords1)
    contact_dists = dists[:l1,l1:] #upper triangular --> first dim = chain 1
    contacts = np.argwhere((contact_dists >= min_distance) & (contact_dists <= max_distance))
    if contacts.shape[0]<1:
        pdockq=0
        ppv=0
    else:
        #Get the average interface plDDT
        avg_if_plddt = np.average(np.concatenate([plddt1[np.unique(contacts[:,0])], plddt2[np.unique(contacts[:,1])]]))
        #Get the number of interface contacts
        n_if_contacts = contacts.shape[0]
        x = avg_if_plddt*np.log10(n_if_contacts)
        pdockq = 0.724 / (1 + np.exp(-0.052*(x-152.611)))+0.018
        #PPV
        PPV = np.array([0.98128027, 0.96322524, 0.95333044, 0.9400192 ,
            0.93172991, 0.92420274, 0.91629946, 0.90952562, 0.90043139,
            0.8919553 , 0.88570037, 0.87822061, 0.87116417, 0.86040801,
            0.85453785, 0.84294946, 0.83367787, 0.82238224, 0.81190228,
            0.80223507, 0.78549007, 0.77766077, 0.75941223, 0.74006263,
            0.73044282, 0.71391784, 0.70615739, 0.68635536, 0.66728511,
            0.63555449, 0.55890174])
        pdockq_thresholds = np.array([0.67333079, 0.65666073, 0.63254566, 0.62604391,
            0.60150931, 0.58313803, 0.5647381 , 0.54122438, 0.52314392,
            0.49659878, 0.4774676 , 0.44661346, 0.42628389, 0.39990988,
            0.38479715, 0.3649393 , 0.34526004, 0.3262589 , 0.31475668,
            0.29750023, 0.26673725, 0.24561247, 0.21882689, 0.19651314,
            0.17606258, 0.15398168, 0.13927677, 0.12024131, 0.09996019,
            0.06968505, 0.02946438])
        inds = np.argwhere(pdockq_thresholds>=pdockq)
        if len(inds)>0:
            ppv = PPV[inds[-1]][0]
        else:
            ppv = PPV[0]
    return pdockq, ppv

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find interactions between amino acids in different chains in all PDB files in a directory.")
    # Positional arguments for the PDB directory and output files
    parser.add_argument("pdb_dir", type=str, help="Path to the directory containing PDB files")
    parser.add_argument("interactions_output", type=str, help="Path for output interactions CSV")
    parser.add_argument("proximity_output", type=str, help="Path for output unrealistic atoms CSV")
    # Optional arguments for minimum and maximum distance thresholds
    parser.add_argument("--min_distance", type=float, default=1.5, 
                        help="Minimum distance threshold in angstroms (default: 1.5 Å)")
    parser.add_argument("--max_distance", type=float, default=5.0, 
                        help="Maximum distance threshold in angstroms (default: 5.0 Å)")
    args = parser.parse_args()
    
    pdb_dir = os.path.abspath(args.pdb_dir)
    interactions_output = os.path.abspath(args.interactions_output)
    proximity_output = os.path.abspath(args.proximity_output)
    min_distance = args.min_distance
    max_distance = args.max_distance

    all_interactions = []
    too_close_output = []
    pdockq_output = []

    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    total_files = len(pdb_files) 

    for idx, pdb_file in enumerate(pdb_files,start=1):
        ## Generate the complex name based on the PDB file name
        fname = os.path.basename(pdb_file)
        sname = re.sub("_unrelaxed.*$", "", fname)
        print(f"Processing file {idx}/{total_files}: {sname}")
        interaction_df, too_close_df = find_residue_interactions_np(pdb_file, min_distance=min_distance, max_distance=max_distance, complex_name=sname)
        all_interactions.append(interaction_df)
        too_close_output.append(too_close_df)

        ## pDockQ calculations
        chain_coords, chain_plddt = read_pdb(pdb_file)
        pdockq, ppv = calc_pdockq(chain_coords, chain_plddt, min_distance, max_distance)
        pdockq_output.append({
                "complex": sname,
                "pdockq" : pdockq,
                "pdockq_confidence" : ppv,
                "chain_A_plddt_mean": round(statistics.mean(chain_plddt['A']), 3),
                "chain_A_plddt_sd": round(statistics.stdev(chain_plddt['A']), 3),
                "chain_B_plddt_mean": round(statistics.mean(chain_plddt['B']), 3),
                "chain_B_plddt_sd": round(statistics.stdev(chain_plddt['B']), 3)
            })

    if all_interactions:
        combined_df = pd.concat(all_interactions, ignore_index=True)
        print(combined_df)
        combined_df.to_csv(interactions_output, index=False)
    else:
        print("No PDB files found or no interactions detected.")
    if 
    if too_close_output and pdockq_output:
        pdockq_df = pd.DataFrame(pdockq_output)
        proximity_df = pd.concat(too_close_output, ignore_index=True)

        merged_df = proximity_df.merge(pdockq_df, on='complex')
        print(merged_df)
        merged_df.to_csv(scoring_output, index=False)