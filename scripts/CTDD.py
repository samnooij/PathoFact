#!/usr/bin/env python
#_*_coding:utf-8_*_
import argparse
import re, math
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

def Count(aaSet, sequence):
	number = 0
	for aa in sequence:
		if aa in aaSet:
			number = number + 1
	cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
	cutoffNums = [i if i >=1 else 1 for i in cutoffNums]

	code = []
	for cutoff in cutoffNums:
		myCount = 0
		for i in range(len(sequence)):
			if sequence[i] in aaSet:
				myCount += 1
				if myCount == cutoff:
					code.append((i + 1) / len(sequence) * 100)
					break
		if myCount == 0:
			code.append(0)
	return code

def CTDD(fastas):
	group1 = {
		'hydrophobicity_PRAM900101': 'RKEDQN',
		'hydrophobicity_ARGP820101': 'QSTNGDE',
		'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
		'hydrophobicity_PONP930101': 'KPDESNQT',
		'hydrophobicity_CASG920101': 'KDEQPSRNTG',
		'hydrophobicity_ENGD860101': 'RDKENQHYP',
		'hydrophobicity_FASG890101': 'KERSQD',
		'normwaalsvolume': 'GASTPDC',
		'polarity':        'LIFWCMVY',
		'polarizability':  'GASDT',
		'charge':          'KR',
		'secondarystruct': 'EALMQKRH',
		'solventaccess':   'ALFCGIVW'
	}
	group2 = {
		'hydrophobicity_PRAM900101': 'GASTPHY',
		'hydrophobicity_ARGP820101': 'RAHCKMV',
		'hydrophobicity_ZIMJ680101': 'HMCKV',
		'hydrophobicity_PONP930101': 'GRHA',
		'hydrophobicity_CASG920101': 'AHYMLV',
		'hydrophobicity_ENGD860101': 'SGTAW',
		'hydrophobicity_FASG890101': 'NTPG',
		'normwaalsvolume': 'NVEQIL',
		'polarity':        'PATGS',
		'polarizability':  'CPNVEQIL',
		'charge':          'ANCQGHILMFPSTWYV',
		'secondarystruct': 'VIYCWFT',
		'solventaccess':   'RKQEND'
	}
	group3 = {
		'hydrophobicity_PRAM900101': 'CLVIMFW',
		'hydrophobicity_ARGP820101': 'LYPFIW',
		'hydrophobicity_ZIMJ680101': 'LPFYI',
		'hydrophobicity_PONP930101': 'YMFWLCVI',
		'hydrophobicity_CASG920101': 'FIWC',
		'hydrophobicity_ENGD860101': 'CVLIMF',
		'hydrophobicity_FASG890101': 'AYHWVMFLIC',
		'normwaalsvolume': 'MHKFRYW',
		'polarity':        'HQRKNED',
		'polarizability':  'KMHFRYW',
		'charge':          'DE',
		'secondarystruct': 'GNPSD',
		'solventaccess':   'MSPTHY'
	}

	groups = [group1, group2, group3]
	property = (
	'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
	'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
	'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

	encodings = []
	header = ['#']
	for p in property:
		for g in ('1', '2', '3'):
			for d in ['0', '25', '50', '75', '100']:
				header.append(p + '.' + g + '.residue' + d)
	encodings.append(header)

	for i in fastas:
		name, sequence = i[0], re.sub('-', '', i[1])
		code = [name]
		for p in property:
			code = code + Count(group1[p], sequence) + Count(group2[p], sequence) + Count(group3[p], sequence)
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


	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Process fasta input with random forest virulence prediction model')
	parser.add_argument("--file", required=True, help="input fasta file")
	parser.add_argument("--out", dest='outFile',
						help="the generated descriptor file")
	args = parser.parse_args()
	
	fastas = readFasta(args.file)
	encodings = CTDD(fastas)
	outFile = args.outFile if args.outFile != None else 'encoding.tsv'
	savetsv(encodings, outFile)




