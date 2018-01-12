# Step Field CCS Automation Package
This python package aims to find reliable features and compute the CCS values for the standards.

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

## Example
```bash
python ccs_auto_app.py --target_list_file CCS-DB-Automation_workshop/TargetList.tsv --config_file CCS-DB-Automation_workshop/config.xml --data_folder CCS-DB-Automation_workshop/cef-and-metadata
```
