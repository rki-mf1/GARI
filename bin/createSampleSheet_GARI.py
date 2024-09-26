import argparse
import os
import pandas as pd


def parseMetadata(metaPath):
  meta = pd.read_csv(metaPath, sep=",")
  meta = meta.reset_index() #???
  
  metaHash={}
  
  for i,row in meta.iterrows():
    combName = row["CUSTOMER_SAMPLE_ID"] + "_" + row["SEQUENCED_SAMPLE_ID"]
    metaHash[combName] = row["SPECIES"]
  
  return(metaHash)



parser = argparse.ArgumentParser(description="")
parser.add_argument("--d", "-dir", type=str, required=True, help="input directory containing gzipped fastq files")
parser.add_argument("--m", "-metadata", type=str, required=False, help="metadata sheet as output by seqrepo export tool, if not provided classifies all samples as NA")
parser.add_argument("--o", "-out", type=str, required=True, help="output samplesheet needed for GARI")
parser.add_argument("--s", "-split", type=str, default="_", help="delimiter to split the sample name by (left side will be taken as ID)")
args = parser.parse_args()

oWriter=open(args.o, "wt")
oWriter.write("%s,%s,%s,%s\n"%("sample", "fastq_1", "fastq_2", "species"))

if args.m:
  mHash = parseMetadata(args.m)

for file in os.listdir(args.d):
    #print(file)
    if os.path.isdir(args.d + "/" + file):
        for file2 in os.listdir(args.d + "/" + file):
            if file2.endswith(".fastq.gz") or file2.endswith(".fq.gz"):
                if "R1" in file2:
                    file2_r2 = file2.replace("R1", "R2")
                    sID = file2.split(args.s)[2]+"_"+file2.split(args.s)[1]
                    if args.m:
                      metaSpecies = mHash[sID]
                    else:
                      metaSpecies = "NA"
                    fPath=args.d + "/" + file
                    oWriter.write("%s,%s,%s,%s\n"%(sID, fPath + "/" + file2, fPath + "/" +file2_r2, metaSpecies))

    if file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
        if "R1" in file:
            file_r2 = file.replace("R1", "R2")
            sID = file.split(args.s)[2]+"_"+file.split(args.s)[1]
            if args.m:
              metaSpecies = mHash[sID]
            else:
              metaSpecies = "NA"

            oWriter.write("%s,%s,%s,%s\n"%(sID, args.d + "/" + file, args.d + "/" + file_r2, metaSpecies))
    


oWriter.close()
