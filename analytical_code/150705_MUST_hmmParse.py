#!/usr/bin/env python

import os
import sys
import argparse
import math

parser = argparse.ArgumentParser(description='Select significant annotations from HMM-output.')
parser.add_argument('inputFile', help='one of the output files from consolidate_hmmscan_results.pl')
parser.add_argument('annotationName', help='what ahould be the name of the annotation? eg. KO,pfamID,tigrID,swissProtEC')
parser.add_argument('-g','--numberOfGenes', type=int,help='number of genes used as input to hmmer, score cut-off is calculated as log2 of this')
parser.add_argument('-s','--scoreCutoff', type=float,default=20.0,help='lower cut-off for score, defaults to 20.0, however use of -g is recommended and -s will be ignored, if -g is used')
parser.add_argument('-p','--percentageCutoff', type=float,default=100.0,help='lower cut-off for keeping multiple annotations per gene; percentage of highest scoring genes, defaults to 80.0; 100.0 would keep only best hit per gene')
parser.add_argument('-v','--validIDs', action='store_true',help='set -v, if a list of all annotations that were assigned should be written to a file')


args = parser.parse_args()
hmmFile = args.inputFile
annN = args.annotationName
if args.numberOfGenes:
    sigVal = math.log(args.numberOfGenes,2)
else:
    sigVal = args.scoreCutoff
multiHitPerc = args.percentageCutoff/100

outFile = hmmFile[:-3] + "besthits.tsv"
if args.validIDs:
    outFile2 = hmmFile[:-3] + "valid" + annN + ".txt"

gene_dict = {}
hmm_file = open(hmmFile, "r")
header = 1
while 1:
    linek = hmm_file.readline()
    if linek == "":
        break
    if header == 1:
        header = 0
    else:
        linek = linek.rstrip()
        tabi = linek.split("\t")
        if float(tabi[2]) >= sigVal:
            tabid, tabgene, tabscore = tabi[0].split("_")[0], tabi[1], float(tabi[2])
            if tabgene not in gene_dict:
                gene_dict[tabgene] = [[], []]
                gene_dict[tabgene][0].append(tabid)
                gene_dict[tabgene][1].append(float(tabscore))
            else:
                if tabscore >= multiHitPerc * gene_dict[tabgene][1][0]:
                    if tabscore > gene_dict[tabgene][1][0]:
                        gene_dict[tabgene][0].insert(0,tabid)
                        gene_dict[tabgene][1].insert(0,float(tabscore))
                    else:
                        gene_dict[tabgene][0].append(tabid)
                        gene_dict[tabgene][1].append(float(tabscore))
hmm_file.close()
out_file = open(outFile, "w")
out_file.write("Gene\t" + annN + "\tmaxScore\thitNumber\n")
allIDs = []
for item in gene_dict:
    gene = item
    priIDs = []
    hN = 0
    score = gene_dict[item][1][0]
    for IDind in range(len(gene_dict[item][0])):
        if gene_dict[item][1][IDind] >= multiHitPerc * score and gene_dict[item][0][IDind] not in priIDs:
            priIDs.append(gene_dict[item][0][IDind])
            if gene_dict[item][0][IDind] not in allIDs:
                allIDs.append(gene_dict[item][0][IDind])
    IDs = ";".join(priIDs)
    hN = len(priIDs)
    out_file.write(gene + "\t" + IDs + "\t" + str(score) + "\t" + str(hN) + "\n")
out_file.close()
if args.validIDs:
    out_file2 = open(outFile2, "w")
    for ids in allIDs:
        out_file2.write(str(ids)+"\n")
    out_file2.close()
print(len(allIDs))
    
