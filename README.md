# Step Field CCS Automation Package
This python package aims to find reliable features and compute the CCS values for chemical standards. It is assumed that prior to running this software you have: 
* preprocessed the raw data files and each of the seven fields is split into a separate file (https://omics.pnl.gov/software/pnnl-preprocessor)
* obtained the correspondent features from Agilent's MP software (a CEF file from each field)

## How to install python packages for this app.
```bash
pip install -r requirements.txt
```

## How to use the app. (python ccs_auto_app.py –help)
```bash
usage: ccs_auto_app.py [-h] [--target_list_file TARGET_LIST_FILE]
                       [--config_file CONFIG_FILE] [--data_folder DATA_FOLDER]
                       [--output OUTPUT]
 
optional arguments:
  -h, --help            show this help message and exit
  --target_list_file TARGET_LIST_FILE
                        Target list file (Tab-delimited text format)
  --config_file CONFIG_FILE
                        Configuration file
  --data_folder DATA_FOLDER
                        Data folder containing all the cef and meta data files
  --output OUTPUT       Output file to save a output table
```

## Demo
You can test this package with the CEF files here included. This is a subset of the standards from the PNNL CCS DB found at https://panomics.pnnl.gov/metabolites/

## How to run
```bash
python ccs_auto_app.py --target_list_file demo/TargetList.tsv --config_file demo/config.xml --data_folder demo/data
```

## Output files
Please refer to the [demo/results](demo/results) folder.
1. [ccs_table.tsv](demo/results/ccs_table.tsv)

  Running this package provides a summary of all about computing CCS values of compounds listed on a target list file (e.g., [TargetList.tsv](demo/TargetList.tsv)). Each row mean a CCS computation results of a single replicate. And each column means the following items.

| Column Names    | Descriptions  |
| -------------   |:-------------:|
| Compound_id     | Compound ID. It will be the same to `CompID` in the target list file |
| name            | Compound Name. It will be the same to `NeutralName` in the target list file |
| Ionization      | Ionization. It will also be the same to `Ionization` in the target list file |
| adduct          | adduct        |
| mass            | mass for the adduct form      |
| ccs_avg         | an average of CCS values of same adducts for all the replicates |
| ccs_rsd         | a relative stdv of CCS values of same adducts for all the replicates |
| ccs             | a CCS value of this adduct for this replicate |
| intensity_org_# | original intensity of a corresponding feature in #th field      |
| intensity_z_#   | z-score of an intensity of a corresponding feature in #th field      |
| mass_error_#    | mass error of a corresponding feature in #th field (ppm)      |
| num_features    | number of features found in all step fields      |
| intercept       | intercept of a CCS regression line      |
| slope           | slope of a CCS regression line      |
| k0              | K0 computed from a CCS regression line      |
| p_value         | p value of a CCS regression line      |
| r_value         | R value of a CCS regression line    |
| replicate       | original file for this CCS computation (.d file)     |


2. <compound_id>_<ionization>_<adduct>.pdf (e.g., [S0000001_neg_\[M-H\].pdf](demo/results/S0000001_neg_\[M-H\].pdf))
  
![demo/results/S0000001_neg_\[M-H\].pdf](demo/results/S0000001_neg_\[M-H\].jpg?raw=true "CCS regression line for \[M-H\]")
3. <compound_id>_<ionization>_intensity_dist.pdf (e.g., [S0000001_neg_intensity_dist.pdf](demo/results/S0000001_neg_intensity_dist.pdf))

   - x-axis: log-scale of original intensities
   - y-axis: kde estimation (normalized distribution)
 
   - dots: features within 15ppm (according to `mz_tolerance` setting in [config.xml](demo/config.xml)) of (e.g., \[M-H\] in the below)
 
   - Vertical lines
      * Most left line: median of feature intensities
      * Center line: 10 \* median
      * Most right line: 2 \* std

![demo/results/S0000001_neg_\[M-H\].pdf](demo/results/S0000001_neg_intensity_dist.jpg?raw=true "Intensity Distribution")
4. <compound_id>_<ionization>_meta.pdf (e.g., [S0000001_neg_meta.pdf](demo/results/S0000001_neg_meta.pdf))
  
  In the figures, dots represent values in the frames of the selected ranges (`frame_offset` in [config.xml](demo/config.xml)) to compute CCS values.
![demo/results/S0000001_neg_\[M-H\].pdf](demo/results/S0000001_neg_meta.jpg?raw=true "Metadata in step fields")
