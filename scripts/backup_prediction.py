#!/usr/bin/env python

import Bio
from Bio import SeqIO
import re
import pandas as pd
import argparse
import numpy as np
from sklearn import model_selection
from sklearn.ensemble import RandomForestClassifier
import pickle

#from PyBioMed.PyProtein import AAComposition

AALetter = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

def CalculateAAComposition(ProteinSequence):
        LengthSequence = len(ProteinSequence)
        Result={}
        for i in AALetter:
                Result[i] = round(float(ProteinSequence.seq.count(i))/ LengthSequence * 100, 2)
        return Result

def CalculateDipeptideComposition(ProteinSequence):
        LengthSequence = len(ProteinSequence)
        Result = {}
        for i in AALetter:
                for j in AALetter:
                        Dipeptide = i + j
                        Result[Dipeptide] = round(float(ProteinSequence.seq.count(Dipeptide)) / (LengthSequence -1) * 100, 2)

        return Result

def CalculateAADipeptideComposition(ProteinSequence):
        Result = {}
        Result.update(CalculateAAComposition(ProteinSequence))
        Result.update(CalculateDipeptideComposition(ProteinSequence))
        return Result


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Process fasta input with random forest virulence prediction model')
	parser.add_argument('infile', metavar='infile', type=str,
		help='The input file (FASTA format)')
	parser.add_argument('outfile', metavar='outfile', type=str,
		help='The outpufile (TSV format)')
	args = parser.parse_args()
	
	df = pd.DataFrame([])
	
	for Sequence in SeqIO.parse(args.infile, "fasta"):	
		ProteinSequence = Sequence	
		DIP = CalculateAADipeptideComposition(ProteinSequence)
		data = DIP
		data = pd.DataFrame(data.items())
		data = data.transpose()

		data.columns=data.iloc[0]
		data = data.drop(data.index[[0]])
		data['Sequence'] = Sequence.id
		df =df.append(data)

	rfc = pickle.load(open("finalized_model_3f.sav", 'rb'))
	X_predict = df.drop('Sequence', axis = 1)
	rfc_predict = rfc.predict(X_predict)	
	rfc_probs = rfc.predict_proba(X_predict)[:,1]
	
	Prediction = pd.DataFrame([])
	dfToList = df['Sequence'].tolist()
	Prediction['Sequence'] = dfToList
	Prediction['Prediction'] = rfc_predict
	Prediction['Probability'] = rfc_probs
	Prediction.to_csv(args.outfile, sep='\t')	



