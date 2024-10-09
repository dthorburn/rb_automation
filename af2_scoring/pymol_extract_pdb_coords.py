#!/usr/bin/env python3

## Use the pymolfree conda environment

## Load libraries
import pymol
from pymol import cmd
import pandas as pd
import os
import re
import sys
import argparse
import glob

## Initialize pymol - '-c' for no GUI, '-q' for quiet mode
pymol.finish_launching(['pymol', '-cq'])

## FIrst attempt. Way too long to process. 
#    ## Iterate through each atom in the model and compare distances
#    for i, atom1 in enumerate(all_residues.atom):
#        for j, atom2 in enumerate(all_residues.atom):
#            if i < j:  # Avoid redundant pair checks
#                distance = cmd.get_distance(f"index {atom1.index}", f"index {atom2.index}")
#                
#                # Check if the distance is within the specified range
#                if 0 < distance <= distance_threshold:
#                    interactions.append({
#                        "complex": sname,
#                        "residue1": f"{atom1.resn}{atom1.resi}",
#                        "chain1": atom1.chain,
#                        "atom1": atom1.name,
#                        "residue1_num": atom1.resi,  # Add residue1 number
#                        "residue2": f"{atom2.resn}{atom2.resi}",
#                        "chain2": atom2.chain,
#                        "atom2": atom2.name,
#                        "residue2_num": atom2.resi,  # Add residue2 number
#                        "distance": distance
#                    })

## Define functions
def find_complex_interactions(pdb_file, distance_threshold=5.0):
    """
    Load a PDB file in PyMOL and find all interactions (within 0-x Å) between amino acids.
    
    Args:
        pdb_file (str): Path to the PDB file.
        distance_threshold (float): Distance cutoff for interaction (in angstroms).
    
    Returns:
        pd.DataFrame: DF with interacting residue pairs and their distances.
    """
    cmd.load(pdb_file, "alphafold_complex")
    
    ## Getting file names for output
    fname = pdb_file.split('/')[-1]
    sname = re.sub("_unrelaxed.*$", "", fname)

    ## Select all residues. Populate empty output. 
    interactions = []

    all_atoms = "all_atoms"
    cmd.select(all_atoms, "all")

    ## Iterate through residues and find nearby residues
    resi_list = []  # Define the list at a higher scope
    cmd.iterate(all_atoms, 'resi_list.append((resi, resn, chain, index))', space={'resi_list': resi_list})

    ## Iterate over atoms, select nearby atoms, and calculate distances
    for atom in cmd.get_model(all_atoms).atom:
        ## Select nearby atoms using `byres` spatial selection function
        nearby_selection = f"(byres all_atoms within {distance_threshold} of index {atom.index}) and not index {atom.index}"
        
        ## Get the nearby atoms
        nearby_atoms = cmd.get_model(nearby_selection)
        
        for nearby_atom in nearby_atoms.atom:
            distance = cmd.get_distance(f"index {atom.index}", f"index {nearby_atom.index}")
            if 0 < distance <= distance_threshold:
                interactions.append({
                    "complex": sname,
                    "residue1": f"{atom.resn}{atom.resi}",
                    "chain1": atom.chain,
                    "atom1": atom.name,
                    "residue1_num": atom.resi,
                    "residue2": f"{nearby_atom.resn}{nearby_atom.resi}",
                    "chain2": nearby_atom.chain,
                    "atom2": nearby_atom.name,
                    "residue2_num": nearby_atom.resi,
                    "distance (Å)": distance
                })

    ## Create a pandas DF from the interaction data. If empty, generate an empty DF.
    if interactions:
        df = pd.DataFrame(interactions)
    else:
        df = pd.DataFrame([{
            "complex": sname,
            "residue1": None,
            "chain1": None,
            "atom1": None,
            "residue1_num": None,
            "residue2": None,
            "chain2": None,
            "atom2": None,
            "residue2_num": None,
            "distance": None
        }])
    
    return df

## Ensures code block execution when called from bash
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find interactions between amino acids in all PDB files in a directory.")
    
    ## Positional arguments for the PDB directory and output file
    parser.add_argument("pdb_dir",     type=str, help="Path to the directory containing PDB files")
    parser.add_argument("output_file", type=str, help="Path for output csv")
    
    ## Optional argument for distance threshold
    parser.add_argument("--distance_threshold", type=float, default=5.0, 
                        help="Distance threshold in angstroms (default: 5.0)")
    
    ## Parse arguments
    args = parser.parse_args()
    
    ## Get the full PDB directory path and distance threshold
    pdb_dir = os.path.abspath(args.pdb_dir)
    output_file = os.path.abspath(args.output_file)
    distance_threshold = args.distance_threshold
    
    ## Empty output dataframe
    all_interactions = []

    ## Find all PDBs in input directory
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    
    for pdb_file in pdb_files:
        fname = pdb_file.split('/')[-1]
        sname = re.sub("_unrelaxed.*$", "", fname)
        print(f"Processing: {sname}") 
        interaction_df = find_complex_interactions(pdb_file, distance_threshold)
        
        # Append the current dataframe to the list of all interactions
        all_interactions.append(interaction_df)
    
    # Concatenate all dfs
    if all_interactions:
        combined_df = pd.concat(all_interactions, ignore_index=True)
    
        # Display the combined DataFrame
        print(combined_df)
    
        # Optionally save the combined DataFrame to a CSV
        combined_df.to_csv(output_file, index=False)
    else:
        print("No PDB files found or no interactions detected.")