#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: wenrchen
"""

import pandas as pd
import sys


def extract_gene(tmp):
    gene_name=tmp.split("|")[-2]
    return gene_name

def extract_uniprot_gene(tmp):
    flag=tmp.find("GN=")
    tmp=tmp[flag+3:]
    flag=tmp.find(" ")
    
    return tmp[:flag]

def extract_ensg(tmp): ##extract gene in "ENSGXXXX" format.
    flag=tmp.find("ENSG")
    tmp=tmp[flag:]
    flag=tmp.find("_")
    
    return tmp[:flag]

df=pd.read_csv(sys.argv[1],sep='\t') ##proteoforms with PTM df

ensg_gene_df=pd.read_csv(sys.argv[3],sep='\t',header=None)
ensg_gene_map={}
for i in range(ensg_gene_df.shape[0]):
    ensg_gene_map[ensg_gene_df.iloc[i][0]]=ensg_gene_df.iloc[i][1]

filtered_df=pd.DataFrame(columns=df.columns)

for i in range(df.shape[0]):
    if(df.iloc[i]['#unexpected modifications']==0):
        filtered_df=filtered_df.append(df.loc[i],ignore_index=True)
    else:
        protein_id=df.iloc[i]['Protein accession']
        ensg=extract_ensg(protein_id)
        gene_name=ensg_gene_map[ensg]
        if(gene_name.find("H1")==-1 and gene_name.find("H2")==-1 and gene_name.find("H3")==-1 and gene_name.find("H4")==-1):
            filtered_df=filtered_df.append(df.loc[i],ignore_index=True)

filtered_df.to_csv(sys.argv[2],sep='\t',index=None)

