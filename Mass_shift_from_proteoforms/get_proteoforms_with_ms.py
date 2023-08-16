#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 10:38:17 2023

@author: wenrchen
"""

##get proteoform df with mass shifts including n-terminal acetylation

import pandas as pd
import sys


df=pd.read_csv(sys.argv[1],sep='\t')

df=df.loc[(df['#unexpected modifications']!=0)| (df['#variable PTMs']==1)| ((df['Protein N-terminal form']!="NONE") & (df['Protein N-terminal form']!="NME") )]

# drop_list=[]
# for i in range(df.shape[0]):
#     if(df.iloc[i]['#unexpected modifications']==0):
#         if(df.iloc[i]['Protein N-terminal form']=="None"):
#             drop_list.append(i)
# df=df.drop(drop_list)
print(df.shape[0])

df.to_csv(sys.argv[2],sep='\t',index=None)