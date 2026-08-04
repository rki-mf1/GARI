#!/usr/bin/env python

import pandas as pd
import json
import argparse as ap 
import ete3

# expected Kraken outputs:   IGS-ID.output.txt

# Kraken target and host reads vs assembly

### INPUT ###
parser = ap.ArgumentParser(description="")
parser.add_argument("-t", "--thresholds", type=str)
parser.add_argument("-s", "--samplesheet", type=str)
parser.add_argument("-d", "--taxdumpDB", type=str, required=False)
args = parser.parse_args()

#### GET KRAKEN OUTPUT
def parseThresholds(inPath, species):
  t_hash ={}
  with open(inPath, 'rt') as f:
    data = json.load(f)

  t_hash = data["default"] # load default params

  if species in data: # if some specific thresholds for a given species exist overwrite the default ones
    for thresh in data[species]:
      t_hash[thresh] = data[species][thresh]
  elif " " in species: # check if the species is really a species name if so check if genus has thresholds
    genus = species.split(" ")[0]
    if genus in data: # catch unknown genus names
      for thresh in data[genus]:
        t_hash[thresh] = data[genus][thresh]   

  return t_hash


def get_kraken(sampleID, kraken_path, tax_target, tax_host):

    with open(sampleID + ".classifiedreads.normalized.txt", 'w') as out:

        if args.taxdumpDB:
          ncbi = ete3.NCBITaxa(dbfile=args.taxdumpDB+"/taxa.sqlite")
        else:
          ncbi = ete3.NCBITaxa(dbfile="taxa.sqlite")

        kraken_output = pd.read_csv(kraken_path, sep='\t', names=['STATUS', 'CONTIG_ID', 'TAX_ID', 'LENGTH', 'INFO'])

        all_length = kraken_output['LENGTH'].sum()
        all_count = len(kraken_output)
        
        for tax_id in [tax_target, tax_host, 0]:
            
            if tax_id != 0:
                tax_ids = ncbi.get_descendant_taxa( tax_id,  intermediate_nodes=True )
                tax_ids.append(tax_id) # else parent id is not included
            else: 
                tax_ids = [0]
            
            kraken_id = kraken_output[kraken_output['TAX_ID'].isin(tax_ids)]
            length, count = kraken_id['LENGTH'].sum(), len(kraken_id)
            
            try:
              perc_norm = round((length / all_length * 100), 2)
            except ZeroDivisionError:
              perc_norm = 0.0

            out.write(f"{perc_norm}\tx\tx\tx\t{tax_id}\tx\n")


samplesheetIn = pd.read_csv(args.samplesheet, sep=',', header=0)
for i, row in samplesheetIn.iterrows():
  print(row)
  dataHash={}
  thresholds = parseThresholds(args.thresholds, row["species"])
  get_kraken(row["sample"], row["kraken2"], thresholds["kraken2_targetID"], thresholds["kraken2_hostID"])