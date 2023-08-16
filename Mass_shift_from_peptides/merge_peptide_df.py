#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 09:48:12 2022

@author: wenrchen
"""

###find start and end position for peptides

import pandas as pd 
import sys

df=pd.read_csv(sys.argv[1],sep='\t')
df_mod=pd.read_csv(sys.argv[2],sep='\t')

pos_dict={}
for i in range(df.shape[0]):
    pos_dict[df.iloc[i]['id']]=(df.iloc[i]['Start position'],df.iloc[i]['End position'])


start_list=[]
end_list=[]
for i in range(df_mod.shape[0]):
    pep_id=df_mod.iloc[i]['Peptide ID']
    start_list.append(pos_dict[pep_id][0])
    end_list.append(pos_dict[pep_id][1])
print(df_mod.shape[0],len(start_list))

df_mod.insert(2,"Start position", start_list)
df_mod.insert(3,"End position", end_list)

df_mod.to_csv(sys.argv[3],sep='\t',index=None)

                        