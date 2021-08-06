#!/usr/bin/env python
#_*_coding:utf-8_*_
import argparse
import re
from collections import Counter
import numpy
import itertools

def readFasta(file):
	with open(file) as f:
		records = f.read()
	if re.search('>', records) == None:
		print('The input file seems not in fasta format.')
		sys.exit(1)
	records = records.split('>')[1:]
	myFasta = []
	for fasta in records:
		array = fasta.split('\n')
		name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '-', ''.join(array[1:]).upper())
		myFasta.append([name, sequence])
	return myFasta

def DPC(fastas):
	AA = 'ACDEFGHIKLMNPQRSTVWY'
	encodings = []
	diPeptides = [aa1 + aa2 for aa1 in AA for aa2 in AA]
	header = ['#'] + diPeptides
	encodings.append(header)

	AADict = {}
	for i in range(len(AA)):
		AADict[AA[i]] = i
	
	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		tmpCode = [0] * 400
		for j in range(len(sequence) - 2 + 1):
			tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] +1
		if sum(tmpCode) != 0:
			tmpCode = [i/sum(tmpCode) for i in tmpCode]
		code = code + tmpCode
		encodings.append(code)
	return encodings

def savetsv(encodings, file = 'encoding.tsv'):
	with open(file, 'w') as f:
		if encodings == 0:
			f.write('Descriptor calculation failed.')
		else:
			for i in range(len(encodings[0]) -1):
				f.write(encodings[0][i] + '\t')
			f.write(encodings[0][-1] + '\n')
			for i in encodings[1:]:
				f.write(i[0] + '\t')
				for j in range(1, len(i) - 1):
					f.write(str(float(i[j])) + '\t')
				f.write(str(float(i[len(i)-1])) + '\n')
	return None


import collections 
def merge(lst1, lst2): 
	return [(sub + [lst2[i][-1]]) for i, sub in enumerate(lst1)] 
       
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Process fasta input with random forest virulence prediction model')
	parser.add_argument("--file", required=True, help="input fasta file")
	parser.add_argument("--out", dest='outFile',
						help="the generated descriptor file")
	args = parser.parse_args()
	
	fastas = readFasta(args.file)
	encodings = DPC(fastas)
	outFile = args.outFile if args.outFile != None else 'encoding.tsv'
	savetsv(encodings, outFile)




