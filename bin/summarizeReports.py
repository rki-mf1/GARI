#!/usr/bin/env python
import argparse
import json
import os
import datetime
import pandas as pd

######################################################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarizes a set of json files for QC parameters")
    parser.add_argument("--g", "-gariOut", type=str, help="GARI output directory with folder per sample")
    parser.add_argument("--p", "-prefix", type=str, help="prefix for table output")
    args = parser.parse_args()

    #timeStamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    #outfile = "GARI_QC_report_" + timeStamp + ".csv"
    outprefix = args.p.split("/")[-1]
    outfile = outprefix + "_QC_report.csv"

    json_list = []

    for ofile in os.listdir(args.g):
        if ofile.endswith(".json"):
            print(args.g + '/' + ofile)

            with open(args.g + '/' + ofile, 'rt') as f:
                data = json.load(f)
                
                print(data)
                # format floats to make sure they always show two digits --> changes them to string...
                for val in data["assembly"]:
                    if isinstance(data["assembly"][val], float) or isinstance(data["assembly"][val], int):
                        data["assembly"][val] = "{:.2f}".format(data["assembly"][val])
                for val in data["reads"]:
                    if isinstance(data["reads"][val], float)or isinstance(data["reads"][val], int):
                        data["reads"][val] = "{:.2f}".format(data["reads"][val])
                for val in data["inferred"]:
                    if isinstance(data["inferred"][val], float) or isinstance(data["inferred"][val], int):
                        data["inferred"][val] = "{:.2f}".format(data["inferred"][val])
                print(data)

                json_list.append(pd.json_normalize(data))

    df = pd.concat(json_list)
    # remove the nested structure of the json and rename columns
    data_renamed = {}
    for col in df:
        #newCol = col.split(".")[-1]
        newCol = col.replace(".", "_")
        data_renamed[col] = newCol
    df_renamed = df.rename(columns=data_renamed)
    

    df_renamed.to_csv(outfile, index=False)