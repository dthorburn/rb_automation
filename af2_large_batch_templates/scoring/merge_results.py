import os
import argparse
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Define secondary structure of input pdbs")
    # Positional arguments for the PDB directory and output file
    parser.add_argument("raw_scores", type=str, help="path to file_scores.csv")
    parser.add_argument("tmscores", type=str, help="path to file_tmscores.csv")
    parser.add_argument("ssscores", type=str, help="path to file_ssfull.csv")
    parser.add_argument("ipscores", type=str, help="path to file_ipsaefull.csv")
    parser.add_argument("output_file", type=str, help="path for output csv")
    args = parser.parse_args()

    temp_scores = os.path.abspath(args.raw_scores)
    temp_tmscores = os.path.abspath(args.tmscores)
    temp_ssscores = os.path.abspath(args.ssscores)
    temp_ipscores = os.path.abspath(args.ipscores)
    output_file = os.path.abspath(args.output_file)

    ##Testing
    #temp_scores="/mnt/c/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/af2_large_batch_templates/scoring/scored/batch99_scores.csv"
    #temp_tmscores="/mnt/c/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/af2_large_batch_templates/scoring/scored/batch99_tmscores.csv"
    #temp_ssscores="/mnt/c/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/af2_large_batch_templates/scoring/scored/batch99_ssfull.csv"
    #temp_ipscores="/mnt/c/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/af2_large_batch_templates/scoring/scored/batch99_ipsaefull.csv"
    
    ## Reading in datasets
    df1 = pd.read_csv(temp_scores)
    df2 = pd.read_csv(temp_tmscores, header=0, names = ["complex", "rank", "pTM", "ipTM"])
    df3 = pd.read_csv(temp_ssscores)
    df3_slim = df3[['complex','rank','ds_bonds_chainA','ds_bonds_chainB','prop_unstruc_A','prop_unstruc_B']]
    df4 = pd.read_csv(temp_ipscores, header=0, names=['Chn1','Chn2','PAE','Dist','Type','ipSAE','ipSAE_d0chn','ipSAE_d0dom','ipTM_af','ipTM_d0chn','pDockQ','pDockQ2','LIS','n0res','n0chn','n0dom','d0res','d0chn','d0dom','nres1','nres2','dist1','dist2','Model'])
    ## Extract complex and rank from Model column, then slim down to relevant columns
    df4['complex'] = df4['Model'].str.split("predictions/").str[1].str.split("_unrelaxed").str[0]
    df4['rank'] = df4['Model'].str.split("_unrelaxed_").str[1].str.split("_alphafold2").str[0]
    df4_slim = df4[['complex', 'rank', 'PAE','Dist', 'ipSAE', 'ipSAE_d0chn', 'ipSAE_d0dom', 'pDockQ','pDockQ2','LIS']]
    ## rename pdockq and pdockq2 to pDockQ_af and pDockQ_d0 for clarity
    df4_slim = df4_slim.rename(columns={'pDockQ': 'pDockQ_ipsae', 'pDockQ2': 'pDockQ_ipsae', 'PAE': 'PAE_cutoff', 'Dist': 'Dist_cutoff'})

    ## Merge output files on complex and rank, then write to output
    m1 = df1.merge(df2, on = ["complex", "rank"])
    m2 = m1.merge(df3_slim, on = ["complex", "rank"])
    m3 = m2.merge(df4_slim, on = ["complex", "rank"])
    m3.to_csv(output_file, index=False)
