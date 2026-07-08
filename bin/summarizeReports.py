#!/usr/bin/env python
import argparse
import json
import os
import pandas as pd

######################################################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarizes a set of json files for QC parameters")
    parser.add_argument("--g", "-gariOut", type=str, help="path to directory with QC json files for each sample")
    parser.add_argument("--p", "-prefix", type=str, help="prefix for table output")
    args = parser.parse_args()

    outprefix = args.p.split("/")[-1]
    outfile = outprefix + "_QC_report.tsv"
    outfile_xlsx = outprefix + "_QC_report.xlsx"

    json_list = []

    for ofile in os.listdir(args.g):
        if ofile.endswith(".json"):
            print(args.g + '/' + ofile)

            with open(args.g + '/' + ofile, 'rt') as f:
                data = json.load(f)
                
                # format floats to make sure they always show two digits --> changes them to string...
                for val in data["assembly"]:
                    if isinstance(data["assembly"][val], float):
                        data["assembly"][val] = "{:.2f}".format(data["assembly"][val])
                    elif isinstance(data["assembly"][val], int):
                        data["assembly"][val] = "{:.0f}".format(data["assembly"][val])                        
                for val in data["reads"]:
                    if isinstance(data["reads"][val], float) :
                        data["reads"][val] = "{:.2f}".format(data["reads"][val])
                    elif isinstance(data["reads"][val], int):
                        data["reads"][val] = "{:.0f}".format(data["reads"][val])
                for val in data["reference"]:
                    if isinstance(data["reference"][val], float):
                        data["reference"][val] = "{:.2f}".format(data["reference"][val])
                    elif isinstance(data["reference"][val], int):
                        data["reference"][val] = "{:.0f}".format(data["reference"][val])
                data["GARI"]["version"] = "v."+ data["GARI"]["version"]

                json_list.append(pd.json_normalize(data))

    df = pd.concat(json_list)
    # remove the nested structure of the json and rename columns
    data_renamed = {}
    for col in df:
        newCol = col.replace(".", "_")
        data_renamed[col] = newCol
    df_renamed = df.rename(columns=data_renamed) 

    df_renamed.to_csv(outfile, index=False, sep="\t")
    df_renamed.to_excel(outfile_xlsx, sheet_name="GARI_QC")