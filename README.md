![Alt text](od_calculator/www/rb_logo.png?raw=true "logo")

# Resurrect Bio Automation and Small Batch Scripts Repo
These pipelines and automation protocols were developed to increase throughput of existing methodologies.

## FeliX Automation
CyBio FeliX liquid handing robot protocols are found in subdirectory `felix_protocols`.

Protocols developed and tested:
	Method 1. Resuspension and washing of agro' colonies.
	Method 2. Dilution of stock resuspensions and aliquoting to platereader plates. 
	Method 3. Futher dilution to target OD from 1:5/1:10 dilution plate.

Protocols in development:
	Method 4. All pairwise mixing of agro' colonies carrying antagonist protiens.
	Method 5. Aliquoting agar and agro' to seedling trays.

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

The notebooks are updated to work with GCP buckets where we store most of our data. These notebooks are suitable for batches of up to 100. Any more, and it becomes viable to use our existing colabfold local solution deployed on private GCP VMs.
