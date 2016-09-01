import pickle
import numpy as np

df = pickle.load(open('eig_values','rb'))
def dates(date):
    d = date.year + date.month/12.0 + date.day/365.25
    return d
for c in df.columns:
    f = open(c,'w')
    for i,t in enumerate(df.index):
        f.write(str(dates(t)) + '\t' + str(df[c][i]) + '\t'+'0.0'+c.split('.')[0] + '\t' + c.split('.')[1] + '\t' + t)
    f.close()
