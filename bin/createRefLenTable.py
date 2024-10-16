#!/usr/bin/env python

import argparse, gzip

def getLenFA(faPath):
    totalLen = 0
    Gs = 0
    Cs = 0
    if faPath.endswith(".gz"):
        faParser=gzip.open(faPath, "rt")
    else:
        faParser=open(faPath, "rt")

    for line in faParser:
        if not line.startswith(">"):
            totalLen += len(line.rstrip())
            Gs += line.rstrip().upper().count("G")
            Cs += line.rstrip().upper().count("C")
    GC= (Gs + Cs) / totalLen

    faParser.close()

    return(totalLen, GC)
#####################################################################################

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Creates a table with lengths per reference in list")
    parser.add_argument("--r", "-refList", type=str, required=True, help="List of references")
    parser.add_argument("--o", "-output", type=str, required=True, help="output for tsv table of lengths")
    args = parser.parse_args()

    oWriter = open(args.o, "wt")

    with open(args.r, 'rt') as f:
        for line in f:         
            refID = line.rstrip().split("/")[-1]
            if not refID == "":
                refLen, refGC = getLenFA(line.rstrip())
                oWriter.write("%s\t%i\t%f\n"%(refID, refLen, refGC))     

    oWriter.close()