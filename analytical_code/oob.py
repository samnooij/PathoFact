#!/usr/bin/env python

import pandas as pd
import numpy as np
import pickle

# Open file with pd.read_csv
df = pd.read_csv("matrix_model_complete.tsv", sep='\t')
#print(df.head())

X = df.drop('class', axis=1)
X = X.drop('#', axis=1)
y = df['class']

#print(X)
#print(y)

from sklearn.model_selection import train_test_split

# implementing train-test-split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=66)


from sklearn.ensemble import RandomForestClassifier 
forest = RandomForestClassifier(n_estimators = 1600, max_depth=340, max_features='auto', oob_score = True)
forest.fit(X_train, y_train)
print('Score: ', forest.score(X_train, y_train))

print('oob_score: ',forest.oob_score_)


