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

from sklearn import model_selection
from sklearn.ensemble import RandomForestClassifier

#random forest model_selection
rfc = RandomForestClassifier(n_estimators=1600, max_depth=340, max_features='auto')
rfc.fit(X_train, y_train)

# predictions
rfc_predict = rfc.predict(X_test)

from sklearn.model_selection import cross_val_score
from sklearn.metrics import classification_report, confusion_matrix

rfc_cv_score= cross_val_score(rfc, X, y, cv=10, scoring='roc_auc')

print("=== Confusion Matrix ===")
print(confusion_matrix(y_test, rfc_predict))
print('\n')
print("=== Classification Report ===")
print(classification_report(y_test, rfc_predict))
print('\n')
print("=== All AUC Scores ===")
print(rfc_cv_score)
print('\n')
print("=== Mean AUC Score ===")
print("Mean AUC Score - Random Forest: ", rfc_cv_score.mean())

#filename = 'finalized_model_3f.sav'
#pickle.dump(rfc, open(filename, 'wb'))
# Tuning Hyperparameters
#from sklearn.model_selection import RandomizedSearchCV

# number of trees in random forest
#n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]

# number of features as every split
#max_features = ['auto', 'sqrt']

# max depth
#max_depth = [int(x) for x in np.linspace(100, 500, num=11)]
#max_depth.append(None)

# create random grid
#random_grid = {
#	'n_estimators': n_estimators,
#	'max_features': max_features,
#	'max_depth': max_depth
#	}

# Random search of parameters
#rfc_random = RandomizedSearchCV(estimator = rfc, param_distributions = random_grid, n_iter = 100, cv = 3, verbose=2, random_state=42, n_jobs = -1)

# Fit the model
#rfc_random.fit(X_train, y_train)

#print(results)
#print(rfc_random.best_params_)

