#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 15:41:29 2023

@author: wenrchen
"""

##remove duplicated PTM sites from MSPathFinder results

import pandas as pd
import sys

argv=sys.argv
file_num=int(argv[1])

df=pd.read_csv(argv[3],sep='\t')
for i in range(file_num-1):
    df_tmp=pd.read_csv(argv[i+4],sep='\t')
    df=df.append(df_tmp,ignore_index=True)
    
print(df.shape[0])

##remain unique ID for the same precursor mass and same Composition
df=df.sort_values(by=["Mass","EValue"],ignore_index=True)

drop_list=[]
scan_dict={}
for i in range(df.shape[0]):
    scan=df.iloc[i]['Scan']
    mass=df.iloc[i]['Mass']
    composition=df.iloc[i]['Composition']
    
    if((mass,composition) not in scan_dict.keys()):
        scan_dict[(mass,composition)]=(df.iloc[i]['EValue'],scan,i)
    else:
        if(scan!=scan_dict[(mass,composition)][1]):
            if(df.iloc[i]['EValue']<scan_dict[(mass,composition)][0]):
                drop_list.append(scan_dict[(mass,composition)][2])
                scan_dict[(mass,composition)]=(df.iloc[i]['EValue'],scan,i)
            else:
                drop_list.append(i)
print(len(drop_list))
df=df.drop(drop_list)

df=df.sort_values(by="EValue",ignore_index=True)

drop_list=[]
protein_ptm_dict={}

for i in range(df.shape[0]):
    protein_id=df.iloc[i]['ProteinName']
    ptm_type=df.iloc[i]['Mod type']
    ptm_pos=df.iloc[i]['Mod global pos']
    
    if(ptm_type.find("Carbam")!=-1):
        drop_list.append(i)
    
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



    