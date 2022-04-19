import pandas as pd
import numpy as np

# smartnet/smartnet_vscode/data/final_data_dock6_knot1.txt.bak
write_object = open('smartnet/smartnet_vscode/data/temp.txt','w')
with open('smartnet/smartnet_vscode/data/final_data_autodock_knot1.txt.bak') as read_object:
    for line in read_object:
        info = line.strip().split('\t')


df = pd.read_csv('smartnet/smartnet_vscode/data/temp.txt', sep='\t', header=None)
df.columns = ['SMILES', 'SEQ']
df2tab = lambda data: pd.get_dummies(data.SMILES).T.dot(data.SEQ.str.get_dummies(','))
tab = df2tab(df)

tab.to_csv('smartnet/smartnet_vscode/data/final_data_autodock_knot1.txt', sep='\t')
#print(tab)