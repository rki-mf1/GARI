#!/usr/bin/env python

import pandas as pd
import json
import argparse as ap 
import ete3

# expected Kraken outputs:   IGS-ID.output.txt

# Kraken target and host reads vs assembly

### INPUT ###
parser = ap.ArgumentParser(description="")
parser.add_argument("-k", "--kraken_output", type=str)
parser.add_argument("-t", "--thresholds", type=str)
parser.add_argument("-s", "--species", type=str)
parser.add_argument("-o", "--outputfile", type=str, default='./KRAKEN_report.txt')
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
    for thresh in data[genus]:
      t_hash[thresh] = data[genus][thresh]   

  return t_hash


def get_kraken(kraken_path, tax_target, tax_host):

    with open(args.outputfile, 'w') as out:

        path=args.kraken_output
        ncbi = ete3.NCBITaxa(dbfile='/scratch/roschekc/ncbi/taxa.sqlite')

        kraken_output = pd.read_csv(path, sep='\t', names=['STATUS', 'CONTIG_ID', 'TAX_ID', 'LENGTH', 'INFO'])

        all_length = kraken_output['LENGTH'].sum()
        all_count = len(kraken_output)
        
        for tax_id in [tax_target, tax_host, 0]:
            
            if tax_id != 0:
                tax_ids = ncbi.get_descendant_taxa( tax_id,  intermediate_nodes=True )
            else: 
                tax_ids = [0]
            
            kraken_id = kraken_output[kraken_output['TAX_ID'].isin(tax_ids)]
            length, count = kraken_id['LENGTH'].sum(), len(kraken_id)
            perc_norm = round((length / all_length * 100), 2)

            out.write(f"{perc_norm}\tx\tx\tx\t{tax_id}\tx\n")


dataHash={}

thresholds = parseThresholds(args.thresholds, args.species)
get_kraken(args.kraken_output, thresholds["kraken2_targetID"], thresholds["kraken2_hostID"])



