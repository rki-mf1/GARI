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

  return t_hash


def parseRefLens(inPath):
  refLens = {}
  with open(inPath, 'rt') as f:
    for line in f:
      splitted = line.rstrip().split("\t")
      refLens[splitted[0]] = [int(splitted[1]), float(splitted[2])]  
  return(refLens)


def parseFASTP(fastpPath):
  dataHash={} 

  with open(fastpPath, 'rt') as f:
    data = json.load(f)
    dataHash['reads_total'] = data['summary']['before_filtering']['total_reads']
    dataHash['reads_past_QC_filter'] = data['summary']['after_filtering']['total_reads']
    dataHash['reads_past_filter_percent'] = round(float(data['summary']['after_filtering']['total_reads'])/float(data['summary']['before_filtering']['total_reads'])*100, 2)
    dataHash['reads_past_filter_bps'] = data['summary']['after_filtering']['total_bases']
    dataHash['reads_past_filter_GC_percent'] = round(data['summary']['after_filtering']['gc_content']*100, 2) 
    dataHash['reads_low_quality'] = data['filtering_result']['low_quality_reads']
    dataHash['reads_too_many_Ns'] = data['filtering_result']['too_many_N_reads']
    dataHash['reads_too_short'] = data['filtering_result']['too_short_reads']
    dataHash['reads_insert_size'] = data['insert_size']['peak']
    dataHash['reads_insert_size_unknown'] = data['insert_size']['unknown']
    dataHash['reads_duplication_rate'] = data['duplication']['rate']

  return(dataHash)


# function to return NA FASTP values for qc_mode
def fakeFASTP():
  dataHash={} 
  dataHash['reads_total'] = 'NA'
  dataHash['reads_past_QC_filter'] = 'NA'
  dataHash['reads_past_filter_percent'] = 'NA'
  dataHash['reads_past_filter_bps'] = 'NA'
  dataHash['reads_past_filter_GC_percent'] = 'NA'
  dataHash['reads_low_quality'] = 'NA'
  dataHash['reads_too_many_Ns'] = 'NA'
  dataHash['reads_too_short'] = 'NA'
  dataHash['reads_insert_size'] = 'NA'
  dataHash['reads_insert_size_unknown'] = 'NA'
  dataHash['reads_duplication_rate'] = 'NA'

  return(dataHash)


def parseBUSCO(buscoPath):

  with open(buscoPath, 'rt') as f:
    data = json.load(f)
    
    BUSCO = data['results']['one_line_summary']
    BUSCO_S = data['results']['Single copy']
    BUSCO_D = data['results']['Multi copy']
    BUSCO_F = data['results']['Fragmented']
    BUSCO_M = data['results']['Missing']

  return([BUSCO, BUSCO_S, BUSCO_D, BUSCO_F, BUSCO_M])


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


def parseFASTANI(fastaniPath):

  with open(fastaniPath, 'rt') as f:
    data = f.readline() # checking first line in enough since ordered by ident
    if data == "": # case for empty fastANI output
      ident = 'NA'
      ref = 'NA'
      aniCov = 'NA'
    else:
      ident = round(float(data.rstrip().split("\t")[2]), 2)
      ref = data.rstrip().split("\t")[1].split("/")[-1]
      aniCov = round(int(data.rstrip().split("\t")[3]) / int(data.rstrip().split("\t")[4]), 2)

  return(ref, ident, aniCov)


def parseAssemblyScan(statsPath):
  # values as they are given by assemblyscan output
  values=["total_contig", "total_contig_length", "max_contig_length", "mean_contig_length", "median_contig_length", "min_contig_length", "n50_contig_length", "l50_contig_count", "num_contig_non_acgtn", "contig_percent_a", "contig_percent_c", "contig_percent_g", "contig_percent_t", "contig_percent_n", "contig_non_acgtn", "contigs_greater_1m", "contigs_greater_100k", "contigs_greater_10k", "contigs_greater_1k", "percent_contigs_greater_1m", "percent_contigs_greater_100k", "percent_contigs_greater_10k", "percent_contigs_greater_1k"]
  values_float=["contig_percent_a", "contig_percent_c", "contig_percent_g", "contig_percent_t", "contig_percent_n", "contig_non_acgtn", "percent_contigs_greater_1m", "percent_contigs_greater_100k", "percent_contigs_greater_10k", "percent_contigs_greater_1k"]
  # values as they should be named for GARI output
  values_renamed=["total_contigs", "assembly_length", "max_contig_length", "mean_contig_length", "median_contig_length", "min_contig_length", "N50", "L50", "num_contig_non_acgtn", "percent_A", "percent_C", "percent_G", "percent_T", "percent_N", "contig_non_acgtn", "contigs_greater_1m", "contigs_greater_100k", "contigs_greater_10k", "contigs_greater_1k", "percent_contigs_greater_1m", "percent_contigs_greater_100k", "percent_contigs_greater_10k", "percent_contigs_greater_1k"]
  
  statHash={}

  with open(statsPath, 'rt') as f:
    data = json.load(f)

    for i,v in enumerate(values):
      if v in values_float:
        statHash[values_renamed[i]] = float(data[v])
      else:
        statHash[values_renamed[i]] = int(data[v])

    statHash["percent_GC"] = round(statHash["percent_G"] + statHash["percent_C"], 2)

  return(statHash)


def parseKRAKEN(krakenPath, t_ID_tar, t_ID_host):
  krakenHash = {}
  krakenHash["Kraken2_target"] = 0
  krakenHash["Kraken2_host"] = 0
  krakenHash["Kraken2_unclassified"] = 0
  with open(krakenPath, 'rt') as f:
    for line in f:
      splitted = line.rstrip().split("\t")
      tarID = int(splitted[4])
      perc = float(splitted[0])

      if tarID == 0:
        krakenHash["Kraken2_unclassified"] = perc
      if tarID == t_ID_tar:
        krakenHash["Kraken2_target"] = perc
      if tarID == t_ID_host:
        krakenHash["Kraken2_host"] = perc

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


def parseRefList(refPath, refID):
  refName = "NA"
  refSp = "NA"

  with open(refPath, "rt") as f:
    for line in f:
      splitted = line.rstrip().split("\t")
      rPath = splitted[0]
      rName = splitted[1]
      
      if rPath.split("/")[-1] == refID:
        refName = rName
        break
  if not refName == "NA":
    refSp = refName.split(" ")[0] + " " + refName.split(" ")[1]

  return(refName, refSp)


# function to check the QC values
def checkThresholds(QCHash, refSpecies, thresh_Hash):

  t_expected_species = args.sp

  # max flags allowed --> more than x flags = FAILED classification
  maxFlagsAllowed = 3

  print("thresholds used:")
  print(thresh_Hash)
  print(refSpecies)

  error = False
  flagged = False
  errorList = []

  # 1. More than X fragments/scaffolds --> very fragmented
  if QCHash["assembly"]["total_contigs"] >= int(thresh_Hash["flag_max_total_contigs"]):
    flagged = True
    errorList.append("Fragmented")

  # 2. BUSCO SingleCopy 
  if QCHash["assembly"]["BUSCO_single_copy"] < float(thresh_Hash["flag_min_BUSCO_S"]):
    errorList.append("BUSCO_single_copy") 
    flagged = True
    if QCHash["assembly"]["BUSCO_single_copy"] < float(thresh_Hash["fail_min_BUSCO_S"]):
      error = True
    

  # 3. BUSCO Duplicates 
  if QCHash["assembly"]["BUSCO_duplicated"] >= float(thresh_Hash["flag_max_BUSCO_D"]):
    flagged = True
    errorList.append("BUSCO_duplicated") 
    if QCHash["assembly"]["BUSCO_duplicated"] >= float(thresh_Hash["fail_max_BUSCO_D"]):
      error = True    
  
  # 4. CheckM completeness
    if QCHash["assembly"]["CheckM_completeness"] < float(thresh_Hash["flag_checkM_complete"]):
      errorList.append("CheckM_completeness") 
      flagged = True
      if QCHash["assembly"]["CheckM_completeness"] < float(thresh_Hash["fail_checkM_complete"]):
        error = True
        
  # 5. CheckM contamination
  if QCHash["assembly"]["CheckM_contamination"] >= float(thresh_Hash["flag_checkM_contamination"]):
    flagged = True
    errorList.append("CheckM_contamination") 
    if QCHash["assembly"]["CheckM_contamination"] >= float(thresh_Hash["fail_checkM_contamination"]):
      error = True

  # 6. Avg. Cov 
  if args.c and QCHash["assembly"]["average_coverage"] < float(thresh_Hash["flag_AvgCov"]):
    flagged = True
    errorList.append("average_coverage")
    if QCHash["assembly"]["average_coverage"] < float(thresh_Hash["fail_AvgCov"]):
      error = True    

  # 7. Perc. Mapped back to ASM
  if args.c and QCHash["assembly"]["percent_reads_mapped"] <= float(thresh_Hash["flag_PercMapped"]):
    flagged = True
    errorList.append("percent_reads_mapped")
    if args.c and QCHash["assembly"]["percent_reads_mapped"] <= float(thresh_Hash["fail_PercMapped"]):
      error = True    

  # 8. No reference found == NA
  if QCHash["assembly"]["reference_file"] == 'NA':
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
  
  # 9. ref identity 
    if QCHash["assembly"]["reference_identity"] <= float(thresh_Hash["flag_ref_ident"]):
      flagged = True
      errorList.append("reference_identity")
      if QCHash["assembly"]["reference_identity"] <= float(thresh_Hash["fail_ref_ident"]):
        error = True

  if args.ck == "true": # only classify kraken2 output if enabled
    # 10. Kraken2 - target
    # check Kraken2 reads
    if args.kr: # only checked if gari not run in QC mode
      if QCHash["Kraken2"]["reads_Kraken2_target"] < float(thresh_Hash["flag_krakenTarget"]):
        flagged=True
        errorList.append("reads_Kraken2_target")
    # check Kraken2 assembly
    if QCHash["Kraken2"]["Kraken2_target"] < float(thresh_Hash["flag_krakenTarget"]):
        flagged=True
        errorList.append("Kraken2_target")

    # 11. Kraken2 - host
    # check Kraken2 reads
    if args.kr: # only checked if gari not run in QC mode
      if QCHash["Kraken2"]["reads_Kraken2_host"] >= float(thresh_Hash["flag_krakenHost"]):
        flagged=True
        errorList.append("reads_Kraken2_host")
    # check Kraken2 assembly
    if QCHash["Kraken2"]["Kraken2_host"] >= float(thresh_Hash["flag_krakenHost"]):
        flagged=True
        errorList.append("Kraken2_host")

  # species specific thresholds (only available if set in QC_thresholds.json for this species)
  # 12. GC range
  if not thresh_Hash["flag_min_GC"] == "NA": # only execute if thresholds not NA
    if not thresh_Hash["flag_min_GC"] <=  QCHash["assembly"]["percent_GC"] <= thresh_Hash["flag_max_GC"]:
      flagged=True
      errorList.append("GC_content")

  # 13. Length range
  if not thresh_Hash["flag_min_length"] == "NA": # only execute if thresholds not NA
    if not thresh_Hash["flag_min_length"] <=  QCHash["assembly"]["assembly_length"] <= thresh_Hash["flag_max_length"]:
      flagged=True
      errorList.append("GC_content")
      if not thresh_Hash["fail_min_length"] <=  QCHash["assembly"]["assembly_length"] <= thresh_Hash["fail_max_length"]:
        error = True
        errorList.append("Assembly_length")

# Total QC warnings allowed before automatical FAIL
  if len(errorList) > maxFlagsAllowed:
    error = True
    errorList.append("too_many_flags")

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
  parser.add_argument("--b", "-busco", type=str, required=True, help="BUSCO report in json format")
  parser.add_argument("--cm", "-checkM", type=str, required=False, help="CheckM tsv report")
  parser.add_argument("--f", "-fastani", type=str, required=True, help="fastANI output")
  parser.add_argument("--s", "-stats", type=str, required=True, help="stats output")
  parser.add_argument("--kr", "-krakenread", type=str, required=False, help="kraken report for reads")
  parser.add_argument("--ka", "-krakenasm", type=str, required=True, help="kraken report for assembly")
  parser.add_argument("--r", "-refLen", type=str, required=True, help="refLen file created by GARI")
  parser.add_argument("--p", "-fastp", type=str, required=False, help="fastp output")
  parser.add_argument("--c", "-cov", type=str, required=False, help="coverage log file of bbmap")
  parser.add_argument("--sp", "-species", type=str, required=False, help="expected species to compare to identified reference")
  parser.add_argument("--gv", "-gariVersion", type=str, required=False, help="GARI version")
  parser.add_argument("--ga", "-gariAssembler", type=str, required=False, help="assembler used in GARI")
  parser.add_argument("--gfDB", "-gariFastaniDB", type=str, required=False, help="GARI fastANI reflist used")
  parser.add_argument("--gbDB", "-gariBuscoLin", type=str, required=False, help="GARI BUSCO lineage used")
  parser.add_argument("--gkDB", "-gariKrakenDB", type=str, required=False, help="GARI Kraken DB used")
  parser.add_argument("--th", "-thresholds", type=str, required=False, help="thresholds file in json format")
  parser.add_argument("--o", "-out", type=str, required=True, help="summary file in json format")
  parser.add_argument("--ck", "-classifyKraken", type=str, required=True, help="if enabled will check if enough reads/contigs have been assigned to the correct phylogentic groups")

  args = parser.parse_args()

  dataHash={}

  thresholds = parseThresholds(args.th, args.sp)

  refHash = parseRefLens(args.r)

  sampleID = args.f.split("_T")[0]
  dataHash["sample_ID"] = sampleID
  dataHash["expected_species"] = args.sp
  dataHash["pipeline"] = {}
  dataHash["pipeline"]["GARI_version"] = args.gv
  dataHash["pipeline"]["assembler"] = args.ga
  dataHash["pipeline"]["FastANI_DB"] = args.gfDB.split("/")[-1]
  dataHash["pipeline"]["BUSCO_lineage"] = args.gbDB
  dataHash["pipeline"]["Kraken2_DB"] = args.gkDB.split("/")[-1]
  dataHash["pipeline"]["Kraken2_target_NCBI_tax_ID"] = thresholds["kraken2_targetID"]
  dataHash["pipeline"]["Kraken2_host_NCBI_tax_ID"] = thresholds["kraken2_hostID"]

  dataHash["assembly"] = parseAssemblyScan(args.s)
  BUSCOs = parseBUSCO(args.b)
  dataHash["assembly"]["BUSCO"] = BUSCOs[0]
  dataHash["assembly"]["BUSCO_single_copy"] = BUSCOs[1]
  dataHash["assembly"]["BUSCO_duplicated"] = BUSCOs[2]
  dataHash["assembly"]["BUSCO_fragmented"] = BUSCOs[3]
  dataHash["assembly"]["BUSCO_missing"] = BUSCOs[4]

  if args.cm:
    checkM_results = parseCheckM(args.cm)
    dataHash["assembly"]["CheckM_marker_lineage"] = checkM_results[0]
    dataHash["assembly"]["CheckM_completeness"] = checkM_results[1]
    dataHash["assembly"]["CheckM_contamination"] = checkM_results[2]
    dataHash["assembly"]["CheckM_strain_heterogeneity"] = checkM_results[3]

  ref, refIdent, refCov = parseFASTANI(args.f)
  dataHash["assembly"]["reference_file"] = ref
  dataHash["assembly"]["reference_identity"] = refIdent
  dataHash["assembly"]["reference_coverage"] = refCov

  if args.p:
    dataHash["reads"] = parseFASTP(args.p)
  else:
    dataHash["reads"] = fakeFASTP()

  # kraken2 assembly (first to init the kraken2 hash)
  dataHash['Kraken2'] = parseKRAKEN(args.ka, thresholds["kraken2_targetID"], thresholds["kraken2_hostID"])
  # kraken2 read
  if args.kr:
    kraken2_read = parseKRAKEN(args.kr, thresholds["kraken2_targetID"], thresholds["kraken2_hostID"])
    dataHash["Kraken2"]["reads_Kraken2_target"] = kraken2_read["Kraken2_target"]
    dataHash["Kraken2"]["reads_Kraken2_host"] = kraken2_read["Kraken2_host"]
    dataHash["Kraken2"]["reads_Kraken2_unclassified"] = kraken2_read["Kraken2_unclassified"]
  else: # for qc_mode add NA for read based info
    dataHash["Kraken2"]["reads_Kraken2_target"] = 'NA'
    dataHash["Kraken2"]["reads_Kraken2_host"] = 'NA'
    dataHash["Kraken2"]["reads_Kraken2_unclassified"] = 'NA'


  #####################################################################################
  # 1) ASM LEN
  asmLen = int(dataHash["assembly"]["assembly_length"])
  if not dataHash["assembly"]["reference_file"] == 'NA':
    refLen = refHash[dataHash["assembly"]["reference_file"]][0]
    lenRatio = round(asmLen / refLen * 100, 2)
  else:
    lenRatio='NA'
  dataHash["assembly"]["reference_length_ratio"] = lenRatio

  #####################################################################################
  # 2) ASM GC
  asmGC = dataHash["assembly"]["percent_GC"]
  if not dataHash["assembly"]["reference_file"] == 'NA':
    refGC = refHash[dataHash["assembly"]["reference_file"]][1]*100 # multiply by 100 to make it comparable
    GCDiff = round(refGC - asmGC, 2)
  else:
    GCDiff = 'NA'

  dataHash["assembly"]["reference_GC_diff"] = GCDiff

  #####################################################################################
  # 3) Coverage and percentage mapped reads
  if args.c:
    avgCov, percMapped = parseCoverage(args.c)
    dataHash["assembly"]["average_coverage"] = avgCov
    dataHash["assembly"]["percent_reads_mapped"] = percMapped
  else:
    dataHash["assembly"]["average_coverage"] = 'NA'
    dataHash["assembly"]["percent_reads_mapped"] = 'NA'

######################################################################################
  # x) derive species name and taxID for the reference found with fastANI
  refName, refSp = parseRefList(args.gfDB, dataHash["assembly"]["reference_file"])
  dataHash["assembly"]["reference_taxon"] = refName
  dataHash["assembly"]["reference_species"] = refSp

######################################################################################
  # 4) Check for threshold violations
  qc_class, qc_errorList = checkThresholds(dataHash, refName, thresholds)
  
  dataHash["QC_status"] = qc_class
  dataHash["QC_warnings"] = "|".join(qc_errorList)
  

######################################################################################
  # X) Write output
  with open(args.o, "wt") as outfile:
    json.dump(dataHash, outfile, indent=4)