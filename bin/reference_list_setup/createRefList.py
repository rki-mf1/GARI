import argparse
import os

parser = argparse.ArgumentParser(description="Script to create the reference list required for GARI based on NCBI genomes")
parser.add_argument("--g", "-genomes", type=str, help="directory with downloaded fna genomes from NCBI")
parser.add_argument("--o", "-output", type=str, help="output reference list required for GARI as tsv file")
args = parser.parse_args()

oWriter = open(args.o , "wt")

for file in os.listdir(args.g):
  if file.endswith(".fna"):
    #print(file)
    with open(args.g + "/" + file, "rt") as myReader:
      firstLine = myReader.readline()
      #print(firstLine.rstrip())
      fPath=os.path.abspath(args.g + "/" + file)
      head_split = firstLine.split(",")[0]
      rName = head_split.split(" ")
      ref_strain = " ".join(rName[1:]).split("chromosome")[0].split("contig")[0].split("scaffold")[0]
      oWriter.write(fPath + "\t" + ref_strain + "\n")
oWriter.close()
