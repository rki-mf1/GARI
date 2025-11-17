#!/usr/bin/env python

#####################################################################################

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


def parseFASTP(fastpPath):
  dataHash={} 

  with open(fastpPath, 'rt') as f:
    data = json.load(f)
    dataHash['total'] = data['summary']['before_filtering']['total_reads']
    dataHash['past_QC_filter'] = data['summary']['after_filtering']['total_reads']
    dataHash['past_filter_percentage'] = round(float(data['summary']['after_filtering']['total_reads'])/float(data['summary']['before_filtering']['total_reads'])*100, 2)
    dataHash['past_filter_bps'] = data['summary']['after_filtering']['total_bases']
    dataHash['past_filter_GC_percentage'] = round(data['summary']['after_filtering']['gc_content']*100, 2) 
    dataHash['low_quality'] = data['filtering_result']['low_quality_reads']
    dataHash['too_many_Ns'] = data['filtering_result']['too_many_N_reads']
    dataHash['too_short'] = data['filtering_result']['too_short_reads']
    dataHash['insert_size'] = data['insert_size']['peak']
    dataHash['insert_size_unknown'] = data['insert_size']['unknown']
    dataHash['duplication_rate'] = data['duplication']['rate']

  return(dataHash)


# function to return NA FASTP values for qc_mode
def fakeFASTP():
  dataHash={} 
  dataHash['total'] = 'NA'
  dataHash['past_QC_filter'] = 'NA'
  dataHash['past_filter_percent'] = 'NA'
  dataHash['past_filter_bps'] = 'NA'
  dataHash['past_filter_GC_percent'] = 'NA'
  dataHash['low_quality'] = 'NA'
  dataHash['too_many_Ns'] = 'NA'
  dataHash['too_short'] = 'NA'
  dataHash['insert_size'] = 'NA'
  dataHash['insert_size_unknown'] = 'NA'
  dataHash['duplication_rate'] = 'NA'

  return(dataHash)

def parseCheckM(checkM_Path):
  with open(checkM_Path, 'rt') as f:
    header = f.readline()
    checkM_report = f.readline()
    checkM_split = checkM_report.rstrip().split("\t")
    marker_lin = checkM_split[1]
    completeness = float(checkM_split[11])
    contamin = float(checkM_split[12])
    heterogenity = float(checkM_split[13])

  return([marker_lin, completeness, contamin, heterogenity])

def parseSKANI(skaniPath):

  with open(skaniPath, 'rt') as f:
    header = f.readline() # checking first line in enough since ordered by ident
    data = f.readline() # checking first line in enough since ordered by ident
    if data == "": # case for empty fastANI output
      ident = 'NA'
      ref_file = 'NA'
      ref_name = 'NA'
      aniCov_ref = 'NA'
      aniCov_query = 'NA'
    else:
      ident = float(data.rstrip().split("\t")[2])
      ref_file = data.rstrip().split("\t")[0].split("/")[-1]
      ref_name = data.rstrip().split("\t")[5]
      aniCov_ref = float(data.rstrip().split("\t")[3])
      aniCov_query = float(data.rstrip().split("\t")[4])

  return(ref_file, ref_name, ident, aniCov_ref, aniCov_query)



def parseAssemblyScan(statsPath):
  # values as they are given by assemblyscan output
  values=["total_contig", "total_contig_length", "max_contig_length", "mean_contig_length", "median_contig_length", "min_contig_length", "n50_contig_length", "l50_contig_count", "num_contig_non_acgtn", "contig_percent_a", "contig_percent_c", "contig_percent_g", "contig_percent_t", "contig_percent_n", "contig_non_acgtn", "contigs_greater_1m", "contigs_greater_100k", "contigs_greater_10k", "contigs_greater_1k", "percent_contigs_greater_1m", "percent_contigs_greater_100k", "percent_contigs_greater_10k", "percent_contigs_greater_1k"]
  values_float=["contig_percent_a", "contig_percent_c", "contig_percent_g", "contig_percent_t", "contig_percent_n", "contig_non_acgtn", "percent_contigs_greater_1m", "percent_contigs_greater_100k", "percent_contigs_greater_10k", "percent_contigs_greater_1k"]
  # values as they should be named for GARI output
  values_renamed=["num_contigs", "length", "max_contig_length", "mean_contig_length", "median_contig_length", "min_contig_length", "N50", "L50", "num_contig_non_acgtn", "A_percentage", "C_percentage", "G_percentage", "T_percentage", "N_percentage", "non_ACGTN_percentage", "contigs_greater_1m", "contigs_greater_100k", "contigs_greater_10k", "contigs_greater_1k", "contigs_greater_1m_percentage", "contigs_greater_100k_percentage", "contigs_greater_10k_percentage", "contigs_greater_1k_percentage"]
  
  statHash={}

  with open(statsPath, 'rt') as f:
    data = json.load(f)

    for i,v in enumerate(values):
      if v in values_float:
        statHash[values_renamed[i]] = float(data[v])
      else:
        statHash[values_renamed[i]] = int(data[v])

    statHash["GC_percentage"] = round(statHash["G_percentage"] + statHash["C_percentage"], 2)

  return(statHash)


def parseKRAKEN(krakenPath, t_ID_tar, t_ID_host):
  krakenHash = {}
  krakenHash["kraken2_target"] = 0
  krakenHash["kraken2_host"] = 0
  krakenHash["kraken2_unclassified"] = 0
  with open(krakenPath, 'rt') as f:
    for line in f:
      splitted = line.rstrip().split("\t")
      tarID = int(splitted[4])
      perc = float(splitted[0])

      if tarID == 0:
        krakenHash["kraken2_unclassified"] = perc
      if tarID == t_ID_tar:
        krakenHash["kraken2_target"] = perc
      if tarID == t_ID_host:
        krakenHash["kraken2_host"] = perc

  return(krakenHash)


def parseCoverage(covPath):
  avgCov="NA"
  percMapped="NA"

  with open(covPath, 'rt') as f:
    for line in f:
      if line.startswith("Percent mapped:"):
        percMapped=round(float(line.rstrip().split(":")[1].rstrip(" ").lstrip(" ")), 2)
      elif line.startswith("Average coverage:"):
        avgCov=round(float(line.rstrip().split(":")[1].rstrip(" ").lstrip(" ")), 2)

  return(avgCov, percMapped)


# function to check the QC values
def checkThresholds(QCHash, refSpecies, thresh_Hash):

  t_expected_species = args.sp

  # max flags allowed --> more than x flags = FAILED classification
  #maxFlagsAllowed = 3

  print("thresholds used:")
  print(thresh_Hash)
  print(refSpecies)

  error = False
  flagged = False
  errorList = []

  # 1. More than X fragments/scaffolds --> very fragmented
  if QCHash["assembly"]["num_contigs"] >= int(thresh_Hash["flag_max_total_contigs"]):
    flagged = True
    errorList.append("num_contigs")

  # 2. CheckM completeness
    if QCHash["assembly"]["checkM_completeness"] < float(thresh_Hash["flag_checkM_complete"]):
      errorList.append("checkM_completeness") 
      flagged = True
      if QCHash["assembly"]["checkM_completeness"] < float(thresh_Hash["fail_checkM_complete"]):
        error = True
        
  # 3. CheckM contamination
  if QCHash["assembly"]["checkM_contamination"] >= float(thresh_Hash["flag_checkM_contamination"]):
    flagged = True
    errorList.append("checkM_contamination") 
    if QCHash["assembly"]["checkM_contamination"] >= float(thresh_Hash["fail_checkM_contamination"]):
      error = True

  # 4. Avg. Cov 
  if args.c and QCHash["assembly"]["average_coverage"] < float(thresh_Hash["flag_AvgCov"]):
    flagged = True
    errorList.append("average_coverage")
    if QCHash["assembly"]["average_coverage"] < float(thresh_Hash["fail_AvgCov"]):
      error = True    

  # 5. Perc. Mapped back to ASM
  if args.c and QCHash["reads"]["mapped_asm_percentage"] <= float(thresh_Hash["flag_PercMapped"]):
    flagged = True
    errorList.append("reads_mapped_asm_percentage")

  # 6. No reference found == NA
  if QCHash["reference"]["accession"] == 'NA':
    error = True
    errorList.append("no_reference")
  else: # check reference values only if there is a ref
    # check if reference is the correct one
    if not t_expected_species.upper() == "NA": # dont throw flags/errors for NA expected species
      if t_expected_species.upper() not in refSpecies.upper():
        # if not a full match try to only match genus
        if t_expected_species.split(" ")[0].upper() in refSpecies.upper():
          flagged = True
          errorList.append("wrong_reference")
        else:
          error = True
          errorList.append("wrong_reference")
  
  # 7. ref identity 
    if QCHash["reference"]["identity"] <= float(thresh_Hash["flag_ref_ident"]):
      flagged = True
      errorList.append("reference_identity")

  if args.ck == "true": # only classify kraken2 output if enabled
    # 8. Kraken2 - target
    # check Kraken2 reads
    if args.kr: # only checked if gari not run in QC mode
      if QCHash["reads"]["tax_class_target_percentage"] < float(thresh_Hash["flag_krakenTarget"]):
        flagged=True
        errorList.append("reads_tax_class_target_percentage")
    # check Kraken2 assembly
    if QCHash["assembly"]["tax_class_target_percentage_normalized"] < float(thresh_Hash["flag_krakenTarget"]):
        flagged=True
        errorList.append("assembly_tax_class_target_percentage")

    # 9. Kraken2 - host
    # check Kraken2 reads
    if args.kr: # only checked if gari not run in QC mode
      if QCHash["reads"]["tax_class_host_percentage"] >= float(thresh_Hash["flag_krakenHost"]):
        flagged=True
        errorList.append("reads_tax_class_host_percentage")
    # check Kraken2 assembly
    if QCHash["assembly"]["tax_class_host_percentage_normalized"] >= float(thresh_Hash["flag_krakenHost"]):
        flagged=True
        errorList.append("assembly_tax_class_host_percentage")

  # species specific thresholds (only available if set in QC_thresholds.json for this species)
  # 10. GC range
  if not thresh_Hash["flag_min_GC"] == "NA": # only execute if thresholds not NA
    if not thresh_Hash["flag_min_GC"] <=  QCHash["assembly"]["GC_percentage"] <= thresh_Hash["flag_max_GC"]:
      flagged=True
      errorList.append("GC_percentage")

  # 11. Length range
  if not thresh_Hash["flag_min_length"] == "NA": # only execute if thresholds not NA
    if not thresh_Hash["flag_min_length"] <=  QCHash["assembly"]["length"] <= thresh_Hash["flag_max_length"]:
      flagged=True
      errorList.append("assembly_length")
      if not thresh_Hash["fail_min_length"] <=  QCHash["assembly"]["length"] <= thresh_Hash["fail_max_length"]:
        error = True
  
  # 12. N50
  if not thresh_Hash["flag_min_N50"] == "NA": # only execute if thresholds not NA
    if QCHash["assembly"]["N50"] < thresh_Hash["flag_min_N50"]:
      flagged=True
      errorList.append("N50")
      #if QCHash["assembly"]["N50"] < thresh_Hash["fail_min_N50"]:
      #  error = True


# Total QC warnings allowed before automatical FAIL
#  if len(errorList) > maxFlagsAllowed:
#    error = True
#    errorList.append("too_many_flags")

  classification = "PASSED"
  if error:
    classification = "FAILED"
  elif flagged:
    classification = "FLAGGED"

  return(classification, errorList) 


#####################################################################################

if __name__ == '__main__':

  import argparse
  import json

  parser = argparse.ArgumentParser(description="Create summary report of GR workflow")
  parser.add_argument("--cm", "-checkM", type=str, required=False, help="checkM tsv report")
  parser.add_argument("--sk", "-skani", type=str, required=True, help="skani output")
  parser.add_argument("--s", "-stats", type=str, required=True, help="stats output")
  parser.add_argument("--kr", "-krakenread", type=str, required=False, help="kraken report for reads")
  parser.add_argument("--ka", "-krakenasm", type=str, required=True, help="kraken report for assembly")
  parser.add_argument("--kan", "-krakenasmNorm", type=str, required=True, help="normalized kraken report for assembly")
  parser.add_argument("--p", "-fastp", type=str, required=False, help="fastp output")
  parser.add_argument("--c", "-cov", type=str, required=False, help="coverage log file of bbmap")
  parser.add_argument("--sp", "-species", type=str, required=False, help="expected species to compare to identified reference")
  parser.add_argument("--gv", "-gariVersion", type=str, required=False, help="GARI version")
  parser.add_argument("--ga", "-gariAssembler", type=str, required=False, help="assembler used in GARI")
  parser.add_argument("--gkDB", "-gariKrakenDB", type=str, required=False, help="GARI Kraken DB used")
  parser.add_argument("--skDB", "-skaniDB", type=str, required=False, help="GARI Skani DB used")
  parser.add_argument("--th", "-thresholds", type=str, required=False, help="thresholds file in json format")
  parser.add_argument("--o", "-out", type=str, required=True, help="summary file in json format")
  parser.add_argument("--ck", "-classifyKraken", type=str, required=True, help="if enabled will check if enough reads/contigs have been assigned to the correct phylogentic groups")

  args = parser.parse_args()

  dataHash={}

  thresholds = parseThresholds(args.th, args.sp)

  sampleID = args.sk.split("_skani.")[0]
  dataHash["sample_ID"] = sampleID
  dataHash["expected_species"] = args.sp
  dataHash["GARI"] = {}
  dataHash["GARI"]["version"] = args.gv
  dataHash["GARI"]["assembler"] = args.ga
  dataHash["GARI"]["skani_DB"] = args.skDB.split("/")[-1]
  dataHash["GARI"]["kraken2_DB"] = args.gkDB.split("/")[-1]
  dataHash["GARI"]["kraken2_target_taxid"] = thresholds["kraken2_targetID"]
  dataHash["GARI"]["kraken2_host_taxid"] = thresholds["kraken2_hostID"]

  dataHash["assembly"] = parseAssemblyScan(args.s)

  if args.cm:
    checkM_results = parseCheckM(args.cm)
    dataHash["assembly"]["checkM_marker_lineage"] = checkM_results[0]
    dataHash["assembly"]["checkM_completeness"] = checkM_results[1]
    dataHash["assembly"]["checkM_contamination"] = checkM_results[2]
    dataHash["assembly"]["checkM_strain_heterogeneity"] = checkM_results[3]

  ref_file, ref_name, ref_ident, aniCov_ref, aniCov_query = parseSKANI(args.sk)
  dataHash["reference"] = {}
  dataHash["reference"]["accession"] = ref_file
  dataHash["reference"]["taxon"] = ref_name
  dataHash["reference"]["species"] = ref_name.split(" ")[1] + " " + ref_name.split(" ")[2]
  dataHash["reference"]["identity"] = ref_ident
  dataHash["reference"]["cov_percentage"] = aniCov_ref
  dataHash["reference"]["cov_asm_percentage"] = aniCov_query

  if args.p:
    dataHash["reads"] = parseFASTP(args.p)
  else:
    dataHash["reads"] = fakeFASTP()

  # kraken2 assembly (first to init the kraken2 hash)
  kraken2_asm = parseKRAKEN(args.ka, thresholds["kraken2_targetID"], thresholds["kraken2_hostID"])
  dataHash["assembly"]["tax_class_target_percentage"] = kraken2_asm["kraken2_target"]
  dataHash["assembly"]["tax_class_host_percentage"] = kraken2_asm["kraken2_host"]
  dataHash["assembly"]["tax_class_unclassified_percentage"] = kraken2_asm["kraken2_unclassified"]

   # kraken2 assembly normalized 
  kraken2_asm_norm = parseKRAKEN(args.kan, thresholds["kraken2_targetID"], thresholds["kraken2_hostID"])
  dataHash["assembly"]["tax_class_target_percentage_normalized"] = kraken2_asm_norm["kraken2_target"]
  dataHash["assembly"]["tax_class_host_percentage_normalized"] = kraken2_asm_norm["kraken2_host"]
  dataHash["assembly"]["tax_class_unclassified_percentage_normalized"] = kraken2_asm_norm["kraken2_unclassified"]

  # kraken2 read
  if args.kr:
    kraken2_read = parseKRAKEN(args.kr, thresholds["kraken2_targetID"], thresholds["kraken2_hostID"])
    dataHash["reads"]["tax_class_target_percentage"] = kraken2_read["kraken2_target"]
    dataHash["reads"]["tax_class_host_percentage"] = kraken2_read["kraken2_host"]
    dataHash["reads"]["tax_class_unclassified_percentage"] = kraken2_read["kraken2_unclassified"]
  else: # for qc_mode add NA for read based info
    dataHash["reads"]["tax_class_target_percentage"] = 'NA'
    dataHash["reads"]["tax_class_host_percentage"] = 'NA'
    dataHash["reads"]["tax_class_unclassified_percentage"] = 'NA'

  #####################################################################################
  # 1) Coverage and percentage mapped reads
  if args.c:
    avgCov, percMapped = parseCoverage(args.c)
    dataHash["assembly"]["average_coverage"] = avgCov
    dataHash["reads"]["mapped_asm_percentage"] = percMapped
  else:
    dataHash["assembly"]["average_coverage"] = 'NA'
    dataHash["reads"]["mapped_asm_percentage"] = 'NA'

######################################################################################
  # 2) Check for threshold violations
  qc_class, qc_errorList = checkThresholds(dataHash, ref_name, thresholds)
  
  dataHash["GARI_QC_status"] = qc_class
  dataHash["GARI_QC_warnings"] = "|".join(qc_errorList)
  

######################################################################################
  # X) Write output
  with open(args.o, "wt") as outfile:
    json.dump(dataHash, outfile, indent=4)