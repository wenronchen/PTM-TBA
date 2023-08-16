#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 23:04:33 2023

@author: wenrchen
"""
##remove histone for mspathfinder

import pandas as pd
import sys


def extract_gene(tmp):
    gene_name=tmp.split("|")[-2]
    return gene_name

def extract_ensg(tmp):
    flag=tmp.find("ENSG")
    tmp=tmp[flag:]
    flag=tmp.find("|")
    
    return tmp[:flag]

df=pd.read_csv(sys.argv[1],sep='\t') ##proteoforms with PTM df


filtered_df=pd.DataFrame(columns=df.columns)

for i in range(df.shape[0]):

    protein_id=df.iloc[i]['ProteinName']

    
    gene_name=extract_gene(protein_id)
    if(gene_name.find("H1")==-1 and gene_name.find("H2")==-1 and gene_name.find("H3")==-1 and gene_name.find("H4")==-1):
        filtered_df=filtered_df.append(df.loc[i],ignore_index=True)

print(filtered_df.shape[0])
filtered_df.to_csv(sys.argv[2],sep='\t',index=None)