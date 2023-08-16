#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: wenrchen
"""

import pandas as pd
import sys

def extract_ptm(tmp):
    while(tmp.find("[")!=-1):
        flag1=tmp.find("[")
        flag2=tmp.find("]")
        tmp_ptm=tmp[flag1+1:flag2]

        tmp=tmp[:flag1]+tmp[flag2+1:]
        
        if(tmp_ptm.find('Carbamidomethylation')==-1 and tmp_ptm.find('Acetyl')==-1):
            
            if(tmp_ptm.find(".")!=-1):
                return float(tmp_ptm)
            else:
                return tmp_ptm
  
def get_ptm_range(tmp,first):
    #remove the letter before "." and the letter after "."

    #remove the letter between "[" and  "]"(include"[" and "]")
    while(tmp.find("[")!=-1):
        
        if(tmp.find("[Acetyl]")!=-1):
            flag=tmp.find("[Acetyl]")
            tmp=tmp[:flag-3]+tmp[flag-2]+tmp[flag+len("[Acetyl]"):]
    
        while(tmp.find("(C)[Carbamidomethylation]")!=-1):
            flag=tmp.find("(C)[Carbamidomethylation]")
            tmp=tmp[:flag]+tmp[flag+1]+tmp[flag+len("(C)[Carbamidomethylation]"):]
        
        flag1=tmp.find("[")
        flag2=tmp.find("]")
        tmp=tmp[:flag1]+tmp[flag2+1:]
        
    flag=tmp.find(".")
    tmp=tmp[flag+1:]
    flag=tmp.find(".")
    tmp=tmp[:flag]
    
    flag1=tmp.find("(")
    flag2=tmp.find(")")

    return (flag1+first,flag2+first-2)        

df=pd.read_csv(sys.argv[1],sep='\t') ##proteoforms with PTM df without histone 

ptm_value=[]
ptm_range=[]
for i in range(df.shape[0]):
    seq=df.iloc[i]['Proteoform']
    first=df.iloc[i]['First residue']
    ptm_value.append(extract_ptm(seq))
    ptm_range.append(get_ptm_range(seq, first))

df.insert(18,"Mass shift", ptm_value)
df.insert(19,"PTM Range",ptm_range)
df.to_csv(sys.argv[2],sep='\t',index=None)
