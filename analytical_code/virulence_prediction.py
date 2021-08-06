#!/usr/bin/env python

import re
import pandas as pd
import argparse
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import pickle


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Process fasta input with random forest virulence prediction model')
	parser.add_argument('infile', metavar='infile', type=str,
		help='The input file (FASTA format)')
	parser.add_argument('outfile', metavar='outfile', type=str,
		help='The outpufile (TSV format)')
	args = parser.parse_args()

	df = pd.read_csv(args.infile, sep='\t')
	X = df.drop('#', axis=1)
	
	rfc = pickle.load(open("scripts/Virulence_factor_model.sav", 'rb'))
	rfc_predict = rfc.predict(X)	
	rfc_probs = rfc.predict_proba(X)[:,1]
	
	Prediction = pd.DataFrame([])
	dfToList = df['#'].tolist()
	Prediction['Sequence'] = dfToList
	Prediction['Prediction'] = rfc_predict
	Prediction['Probability'] = rfc_probs
	Prediction.to_csv(args.outfile, sep='\t', header=False)	



