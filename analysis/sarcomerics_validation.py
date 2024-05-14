import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import os
from sklearn import preprocessing
from math import *
import re
import plotly.graph_objects as go
if('analysis' in os.getcwd()):
    os.chdir('../')
else:
    print(os.getcwd())
    datatypes = {'POS': 'int64',
             'REF': 'str',
             'ALT': 'str'}
# alpha missense data
ams_prediction = pd.read_csv('_in/sarcomeric_annotated_ams.csv',header=0, dtype=datatypes)

# VEP annotation data (other predictors)
vep_prediction = pd.read_csv('VEP_format/vep_results/vep_prediction_scores.csv',
                             header = 0, sep = ',', dtype=datatypes,na_values="-")
merged = pd.merge(ams_prediction,vep_prediction, on = ['CHROM','POS','ALT'], how='left')

# Change CADD scale to be in the 0-1 range

scaler = preprocessing.MinMaxScaler()
a = np.array(range(1,99)).reshape(-1,1)
scaler.fit(a)
merged['CADD_scaled'] = scaler.transform(merged['CADD_phred'].to_numpy().reshape(-1,1))
merged['CADD_scaled'].describe()

scaler = preprocessing.MinMaxScaler()
a = np.array(range(floor(min(merged['FATHMM_score'])),ceil(max(merged['FATHMM_score'])))).reshape(-1,1)
scaler.fit(a)
merged['FATHMM_scaled'] = scaler.transform(merged['FATHMM_score'].to_numpy().reshape(-1,1))
merged['FATHMM_scaled'].describe()

roc_df = merged[merged['Class'] != 0.5] # remove variants non labeled as 0 (B/LB) or 1 (P/LP)

predictors = ['am_pathogenicity','CADD_scaled','REVEL_score','Polyphen2_HVAR_score','MutationTaster_score','MutPred_score',
              'SIFT4G_score','FATHMM_scaled','DANN_score','MetaLR_score','MetaRNN_score','PrimateAI_score']
pred = ['Class'] + predictors
roc_data = roc_df[~roc_df.loc[:,pred].isnull().any(axis=1)].loc[:,pred]