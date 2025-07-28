#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 10:31:42 2023

@author: wenrchen
"""

##extract ptm information (mass and localization) from the PSM table of MSFragger results

import pandas as pd
import sys
import warnings

warnings.simplefilter("ignore")

def contain_lower(seq): ##check whether the sequence contains lower cases
    if(seq!=0):
        for s in seq:
            if(s.islower()):
                return True
    return False

df=pd.read_csv(sys.argv[1],sep='\t') ##PSM with PTM df
df=df.fillna(0)

output=sys.argv[2]

df=df.sort_values(by="Expectation",ignore_index=True)
df_ms=pd.DataFrame(columns=df.columns)
ptm_value,ptm_start,ptm_end,ptm_type=[],[],[],[]

n_term_psm_cnt=0
oxidation_psm_cnt=0
for i in range(0,df.shape[0]):
    
    if(str(df.iloc[i]['Assigned Modifications']).find("N-term")!=-1):
        n_term_psm_cnt+=1
    if(str(df.iloc[i]['Assigned Modifications']).find("M")!=-1):
        oxidation_psm_cnt+=1
    
    start=df.iloc[i]['Protein Start']
    assign_mods=str(df.iloc[i]['Assigned Modifications']).split(",")
    for assign_mod in assign_mods:                              ##assign mods are variable PTMs
        if(assign_mod.find("N-term")!=-1):
            ptm_type.append("N-term acetylation")
            ptm_start.append(start)
            ptm_end.append(start)
            
            flag=assign_mod.find("(")
            ptm_value.append(float(assign_mod[flag+1:-1])) ##42.0106 Da
            df_ms=df_ms.append(df.loc[i],ignore_index=True)
        elif(assign_mod.find("M")!=-1):
            ptm_type.append("Oxidation")
            flag=assign_mod.find("M")
            ptm_start.append(start+int(assign_mod[:flag])-1)
            ptm_end.append(start+int(assign_mod[:flag])-1)
            
            flag=assign_mod.find("(")
            ptm_value.append(float(assign_mod[flag+1:-1])) ##15.9949 Da
            df_ms=df_ms.append(df.loc[i],ignore_index=True)
    if(abs(df.iloc[i]['Delta Mass'])>=0.1):
        obs_mods=str(df.iloc[i]['Observed Modifications']).split(";")
        for obs_mod in obs_mods:
            ptm_type.append(obs_mod)
            ptm_value.append(df.iloc[i]['Delta Mass'])
            
            local=df.iloc[i]['MSFragger Localization']
            if(contain_lower(local)):

                ptm_indice=[idx for idx, aa in enumerate(local) if aa.islower()]
                ptm_start.append(start+ptm_indice[0])
                ptm_end.append(start+ptm_indice[-1])
                
            else:
                ptm_start.append(start)
                ptm_end.append(df.iloc[i]['Protein End'])
            df_ms=df_ms.append(df.loc[i],ignore_index=True)

print(n_term_psm_cnt,oxidation_psm_cnt)
print(len(ptm_start),len(ptm_value),len(ptm_type))
print(df_ms.shape[0])
df_ms.insert(26,"Mod start",ptm_start)
df_ms.insert(27,"Mod end",ptm_end)
df_ms.insert(28,"Mod type",ptm_type)
df_ms.insert(29,"Mass shift",ptm_value)

psm_list=list(df_ms['Spectrum'].values)
print(len(list(set(psm_list))))

df_ms.to_csv(output,sep='\t',index=None)
            
            
        
                
                
                
    
