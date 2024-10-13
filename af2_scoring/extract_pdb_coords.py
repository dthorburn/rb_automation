import MDAnalysis as mda
import pandas as pd
import numpy as np
import os
import re
import argparse
import glob

def find_residue_interactions(pdb_file, min_distance=1.5, max_distance=5.0, complex_name=None):
    ## If complex_name is not provided, generate it from the file name
    if complex_name is None:
        fname = os.path.basename(pdb_file)
        complex_name = re.sub("_unrelaxed.*$", "", fname)

    ## Load the structure using MDAnalysis
    u = mda.Universe(pdb_file)
    interactions = []

    ## Iterate through each residue and find nearby residues from different chains
    for residue in u.residues:
        ## Select atoms within max_distance from the current residue
        #nearby_residues = u.select_atoms(f"around {max_distance} resid {residue.resid}")
        nearby_residues = u.select_atoms(f"around {max_distance} (not name H* and resid {residue.resid})")

        for nearby_residue in nearby_residues.residues:
            ## Only check interactions on different chains
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
        df_subset = df[df["chain1"] == "A"]
    else:
        df_subset = pd.DataFrame([{
            "complex": complex_name,
            "residue1": None,
            "chain1": None,
            "residue1_num": None,
            "residue2": None,
            "chain2": None,
            "residue2_num": None,
            "distance": None
        }])

    return df_subset

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find interactions between amino acids in different chains in all PDB files in a directory.")
    # Positional arguments for the PDB directory and output file
    parser.add_argument("pdb_dir", type=str, help="Path to the directory containing PDB files")
    parser.add_argument("output_file", type=str, help="Path for output CSV")
    # Optional arguments for minimum and maximum distance thresholds
    parser.add_argument("--min_distance", type=float, default=1.5, 
                        help="Minimum distance threshold in angstroms (default: 1.5 Å)")
    parser.add_argument("--max_distance", type=float, default=5.0, 
                        help="Maximum distance threshold in angstroms (default: 5.0 Å)")
    args = parser.parse_args()
    
    pdb_dir = os.path.abspath(args.pdb_dir)
    output_file = os.path.abspath(args.output_file)
    min_distance = args.min_distance
    max_distance = args.max_distance

    all_interactions = []

    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    total_files = len(pdb_files) 

    for idx, pdb_file in enumerate(pdb_files,start=1):
        # Generate the complex name based on the PDB file name
        fname = os.path.basename(pdb_file)
        sname = re.sub("_unrelaxed.*$", "", fname)
        print(f"Processing file {idx}/{total_files}: {sname}")
        interaction_df = find_residue_interactions(pdb_file, min_distance=min_distance, max_distance=max_distance, complex_name=sname)
        all_interactions.append(interaction_df)
    
    if all_interactions:
        combined_df = pd.concat(all_interactions, ignore_index=True)
        print(combined_df)
        combined_df.to_csv(output_file, index=False)
    else:
        print("No PDB files found or no interactions detected.")
