#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 09:17:54 2023

@author: wenrchen
"""
##remove duplicated PTM sites from MetaMorphous results

import pandas as pd
import sys

argv=sys.argv

df=pd.read_csv(argv[1],sep='\t')

df=df.sort_values(by=["PEP_QValue"],ignore_index=True)

drop_list=[]
protein_ptm_dict={}

for i in range(df.shape[0]):
    protein_id=df.iloc[i]['Protein Accession']
    ptm_type=df.iloc[i]['Mod type']
    ptm_pos=df.iloc[i]['Mod start']
    
    if(protein_id not in protein_ptm_dict.keys()):
        protein_ptm_dict[protein_id]=[(ptm_type,ptm_pos)]
    else:
        ptm_list=protein_ptm_dict[protein_id]
        if((ptm_type,ptm_pos) in ptm_list):
            drop_list.append(i)
        else:
            protein_ptm_dict[protein_id].append((ptm_type,ptm_pos))
            
df=df.drop(drop_list)
print(df.shape[0])
df.to_csv(argv[2],sep='\t',index=None)

