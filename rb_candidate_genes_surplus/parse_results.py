import os
import re
import argparse
import pandas as pd
from tqdm import tqdm

parser = argparse.ArgumentParser(description="Parse output of TMBed to summary CSV.")
# Positional arguments for the input and output files.
parser.add_argument("input_file", type=str, help="Path to the TMBed output file (out-format must be 0).")
parser.add_argument("output_file", type=str, help="Path for output CSV.")
args = parser.parse_args()

#in_file = "/home/dthorbur/Resurrect_Bio/Tools/tmbed/examples/Phapa1_MT2006_JGI.aa.tmbed.out"
#out_file = "/home/dthorbur/Resurrect_Bio/Tools/tmbed/examples/Phapa1_MT2006_JGI.aa.tmbed.results.csv"
in_file = os.path.abspath(args.input_file)
out_file = os.path.abspath(args.output_file)

with open(in_file, "r") as tmbd_file:
    ## Removing empty lines and leading/trailing white spaces
    not_empty_lines = [line.strip() for line in tmbd_file if line.strip()]
    ## A loop to read 3 lines at a time since this isn't a conventional fasta
    output = []
    for i in tqdm(range(0, len(not_empty_lines), 3)):
        temp_dat = not_empty_lines[i:i+3]
        #print(temp_dat)
        output.append({
            "seqname": re.sub("^>", "", temp_dat[0]),
            "signalpep": bool(re.match(r'^S', temp_dat[2])),
            "signalpep_length": len(re.match(r'^S*', temp_dat[2]).group()),
            "tmdomains": bool(re.findall(r'[BbHh]', temp_dat[2]))
        })
    if output:
        df = pd.DataFrame(output)
        df.to_csv(out_file, index=False)
        print(f"Output file saves to:", out_file)
        print(df)
    else:
        print("Malformed input or all empty lines. Check input file is correct.")

