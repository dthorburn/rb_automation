![Alt text](od_calculator/www/rb_logo.png?raw=true "logo")

## Resurrect Bio Automation, Small Batch AF2 Notebooks, and Large Batch AF2 Templates
These pipelines and automation protocols were developed to increase throughput of existing methodologies. There is some overlap with the private `af2mm` repo that has nextflow implementations of colabfold, foldseek, and scoring of AF2 multimer complexes.   

## FeliX Automation
CyBio FeliX liquid handing robot protocols are found in subdirectory `felix_protocols`.

Protocols developed and tested:
- Method 1. Resuspension and washing of agro' colonies.
- Method 2. Dilution of stock resuspensions and aliquoting to platereader plates. 
- Method 3. Futher dilution to target OD from 1:5/1:10 dilution plate.
- Method 4. Mixing candidate proteins for infections assay to platereader plates.
- Method 5. All pairwise mixing of agro' colonies carrying antagonist protiens.

Protocols in development:
- Method 6. Aliquoting agar and agro' to seedling trays.

## OD Calculator
Source code for UI and Server for the OD calculator app designed for use with the Thermo Multiskan platereader in Scale Space phase C.

Usage:
1. Select local XLSX file emitted from the platereader.
2. Add a run name, last well in use, and if necessary, change any of the other default calculation parameters.
3. Alter editable target OD table, upload target OD table, or use homogeneous target OD table.
4. Click on Output tab and inspect data.
5. Download summary table which is the input for FeliX.
6. Select new platereader XLSX file following method 3 completion.
7. Click on Platereader table tab to check accuracy of observed ODs.

## AF2 Colabfold Notebook
Notebooks copied from the [colabfold](https://github.com/sokrypton/ColabFold) repo. Coped on 01/02/2024 for version 1.5.5.

The modified notebooks are updated to work with GCP buckets where we store most of our data and found in subdirectory `af2_small_batch`. These notebooks are suitable for batches of up to 100. Any more and it becomes viable to use our existing colabfold local solution deployed on private GCP VMs. There are will be 2 copies of the current small batch notebooks so parallelisation is possible. 

This repo also contains backup script templates for the GCP/AWS colabfold solution in the subdirectory `af2_large_batch_templates`. Instructions to install local colabfold can be found in their own [repo](https://github.com/YoshitakaMo/localcolabfold) with the database setup script found in the original ColabFold [repo](https://github.com/sokrypton/ColabFold).

To use, follow the SOPs instructions on RB's Notion. 

## Resurrect Bio Batch Visualisation of AF2 Monomers and Multimers 
The scripts found in subdirectory `af2_visualisation` are written to batch process a directory of PDB files emitting a PNG file with the monomer or multimer from several different angles. 
To use follow these instructions:
1. Create the conda environment `conda create env --name pymolfree --file pymolfree.yml`
2. Activate the conda environment `conda activate pymolfree`
3. Execute bash script `bash pymol_make_figure.sh /path/to/pdb/dir`
