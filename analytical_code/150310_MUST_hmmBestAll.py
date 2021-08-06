#!/usr/bin/env python

import os
import sys
import argparse
import math

parser = argparse.ArgumentParser(description='Select significant annotations from HMM-output.')
parser.add_argument('koFile', help='KEGG output files from consolidate_hmmscan_results.pl')
parser.add_argument('mcFile', help='metaCyc output files from consolidate_hmmscan_results.pl')
parser.add_argument('spFile', help='Swiss-Prot output files from consolidate_hmmscan_results.pl')
parser.add_argument('pfFile', help='Pfam output files from consolidate_hmmscan_results.pl')
parser.add_argument('tiFile', help='TIGR Pfam output files from consolidate_hmmscan_results.pl')
parser.add_argument('-g','--numberOfGenes', type=int,help='number of genes used as input to hmmer, score cut-off is calculated as log2 of this')
parser.add_argument('-s','--scoreCutoff', type=float,default=20.0,help='lower cut-off for score, defaults to 20.0, however use of -g is recommended and -s will be ignored, if -g is used')

args = parser.parse_args()
koFile = args.koFile
mcFile = args.mcFile
spFile = args.spFile
pfFile = args.pfFile
tiFile = args.tiFile
annN = "ID"
if args.numberOfGenes:
    sigVal = math.log(args.numberOfGenes,2)
else:
    sigVal = args.scoreCutoff

outFile = "besthitsAllDB.tsv"

gene_dict = {}
hmm_file = open(koFile, "r")
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
            tabid, tabgene, tabscore = "KEGG:"+tabi[0].split("_")[0], tabi[1], float(tabi[2])
            if tabgene not in gene_dict:
                gene_dict[tabgene] = [[], []]
                gene_dict[tabgene][0].append(tabid)
                gene_dict[tabgene][1].append(float(tabscore))
            else:
                if tabscore >= gene_dict[tabgene][1][0]:
                    if tabscore > gene_dict[tabgene][1][0]:
                        gene_dict[tabgene][0].insert(0,tabid)
                        gene_dict[tabgene][1].insert(0,float(tabscore))
                    else:
                        gene_dict[tabgene][0].append(tabid)
                        gene_dict[tabgene][1].append(float(tabscore))
hmm_file.close()

hmm_file = open(mcFile, "r")
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
            tabid, tabgene, tabscore = "metaCyc:"+tabi[0].split("_")[0], tabi[1], float(tabi[2])
            if tabgene not in gene_dict:
                gene_dict[tabgene] = [[], []]
                gene_dict[tabgene][0].append(tabid)
                gene_dict[tabgene][1].append(float(tabscore))
            else:
                if tabscore >= gene_dict[tabgene][1][0]:
                    if tabscore > gene_dict[tabgene][1][0]:
                        gene_dict[tabgene][0].insert(0,tabid)
                        gene_dict[tabgene][1].insert(0,float(tabscore))
                    else:
                        gene_dict[tabgene][0].append(tabid)
                        gene_dict[tabgene][1].append(float(tabscore))
hmm_file.close()

hmm_file = open(spFile, "r")
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
            tabid, tabgene, tabscore = "swissProt:"+tabi[0].split("_")[0], tabi[1], float(tabi[2])
            if tabgene not in gene_dict:
                gene_dict[tabgene] = [[], []]
                gene_dict[tabgene][0].append(tabid)
                gene_dict[tabgene][1].append(float(tabscore))
            else:
                if tabscore >= gene_dict[tabgene][1][0]:
                    if tabscore > gene_dict[tabgene][1][0]:
                        gene_dict[tabgene][0].insert(0,tabid)
                        gene_dict[tabgene][1].insert(0,float(tabscore))
                    else:
                        gene_dict[tabgene][0].append(tabid)
                        gene_dict[tabgene][1].append(float(tabscore))
hmm_file.close()

hmm_file = open(pfFile, "r")
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
            tabid, tabgene, tabscore = "Pfam:"+tabi[0], tabi[1], float(tabi[2])
            if tabgene not in gene_dict:
                gene_dict[tabgene] = [[], []]
                gene_dict[tabgene][0].append(tabid)
                gene_dict[tabgene][1].append(float(tabscore))
            else:
                if tabscore >= gene_dict[tabgene][1][0]:
                    if tabscore > gene_dict[tabgene][1][0]:
                        gene_dict[tabgene][0].insert(0,tabid)
                        gene_dict[tabgene][1].insert(0,float(tabscore))
                    else:
                        gene_dict[tabgene][0].append(tabid)
                        gene_dict[tabgene][1].append(float(tabscore))
hmm_file.close()

hmm_file = open(tiFile, "r")
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
            tabid, tabgene, tabscore = "TIGR:"+tabi[0].split("_")[0], tabi[1], float(tabi[2])
            if tabgene not in gene_dict:
                gene_dict[tabgene] = [[], []]
                gene_dict[tabgene][0].append(tabid)
                gene_dict[tabgene][1].append(float(tabscore))
            else:
                if tabscore >= gene_dict[tabgene][1][0]:
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
        if gene_dict[item][1][IDind] >= score and gene_dict[item][0][IDind] not in priIDs:
            priIDs.append(gene_dict[item][0][IDind])
            if gene_dict[item][0][IDind] not in allIDs:
                allIDs.append(gene_dict[item][0][IDind])
    IDs = ";".join(priIDs)
    hN = len(priIDs)
    out_file.write(gene + "\t" + IDs + "\t" + str(score) + "\t" + str(hN) + "\n")
out_file.close()
print(len(allIDs))
    
