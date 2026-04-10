import os
import argparse
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Define secondary structure of input pdbs")
    # Positional arguments for the PDB directory and output file
    parser.add_argument("raw_scores", type=str, help="path to file_scores.csv")
    parser.add_argument("tmscores", type=str, help="path to file_tmscores.csv")
    parser.add_argument("ssscores", type=str, help="path to file_ssfull.csv")
    parser.add_argument("output_file", type=str, help="path for output csv")
    args = parser.parse_args()

    temp_scores = os.path.abspath(args.raw_scores)
    temp_tmscores = os.path.abspath(args.tmscores)
    temp_ssscores = os.path.abspath(args.ssscores)
    output_file = os.path.abspath(args.output_file)

    ##Testing
    #temp_scores="/mnt/c/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/af2_large_batch_templates/scoring/scored/batch99_scores.csv"
    #temp_tmscores="/mnt/c/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/af2_large_batch_templates/scoring/scored/batch99_tmscores.csv"
    #temp_ssscores="/mnt/c/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/af2_large_batch_templates/scoring/scored/batch99_ssfull.csv"
    
    df1 = pd.read_csv(temp_scores)
    df2 = pd.read_csv(temp_tmscores, header=0, names = ["complex", "rank", "pTM", "ipTM"])
    df3 = pd.read_csv(temp_ssscores)
    df3_slim = df3[['complex','rank','ds_bonds_chainA','ds_bonds_chainB','prop_unstruc_A','prop_unstruc_B']]

    m1 = df1.merge(df2, on = ["complex", "rank"])
    m2 = m1.merge(df3_slim, on = ["complex", "rank"])
    m2.to_csv(output_file, index=False)
