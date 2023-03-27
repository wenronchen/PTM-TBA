#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 13:25:06 2022

@author: wenrchen
"""

#find confident oxidations 
import pandas as pd
import sys

def get_ptm_seq(tmp):
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
    
    # while(tmp.find("(C)")!=-1):
    #     flag=tmp.find("(C)")
    #     tmp=tmp[:flag]+tmp[flag+1]+tmp[flag+3:]
        
    
    flag1=tmp.find("(")
    flag2=tmp.find(")")
    
    return tmp[flag1+1:flag2]

def compare_mass(m1, m2, t):##determine if the m1-m2, m1-m2+1.00235, m1-m2-1.00235 satisfy the error tolerance t
    if(abs(m1-m2)<=t):
        return True
    elif(abs(m1-m2-1.00235)<=t):
        return True
    elif(abs(m1-m2+1.00235)<=t):
        return True
    else:
        return False


df=pd.read_csv(sys.argv[1],sep='\t')
#ptm_value=float(sys.argv[2])

#potential_aa=["M","D","K","N","P","Y","R","C"] ##oxidation 15.9949
#potential_aa=["S","T","Y"] ##Phosphorylation 79.9663
#potential_aa=["C","I","L","K","R","H","D","E","N","Q","S","T"] # Methylation 14.0157
#potential_aa=["K","C","S","T"] ##Acetylation 42.0106
#potential_aa=["K","N","R"] ##Dimethylation, 28.0313
#potential_aa=["N","Q"] ##Deamination, 0.9840
#potential_aa=["D","E"] ##Cation:Na, 21.9819; Cation:K, 37.95582

potential_aa_dict={21.981943:["D","E"],37.955882:["D","E"],52.911464:["D","E"],53.919289:["D","E"],\
                   43.005814:["K","R","C","M","S","T","Y"],79.966331:["S","T","Y"],183.035399:["H","K","S","Y"],\
                       14.01565:["C","I","L","K","R","H","D","E","N","Q","S","T"],28.0313:["K","N","R","P"],\
                          42.010565:["K","C","S","T","Y","H","R"], 44.008456:["S"],\
                              15.994915:["M","W","H","D","K","N","P","Y","R","C","G","U","E","I","L","Q","S","T","V"]}

remained=[]
cnt=0
for i in range(df.shape[0]):
    for ptm_value in potential_aa_dict.keys():
        
        if(compare_mass(df.iloc[i]['Mass shift'],ptm_value,0.1)):
            potential_aa=potential_aa_dict[ptm_value]
            cnt+=1
            ptm_range=df.iloc[i]['PTM Range'].split(",")
            ptm_first=int(ptm_range[0][1:])
            ptm_last=int(ptm_range[1][:-1])
            
            ptm_seq=get_ptm_seq(df.iloc[i]['Proteoform'])
            for aa in potential_aa:
                if(aa in ptm_seq):
                    remained.append(i)
                    break
            
        # if((ptm_last-ptm_first)<3):
        #     remained.append(i)

print(cnt)
remained=list(set(remained))
print(len(remained))

output_df=df.loc[remained]
output_df.to_csv(sys.argv[2],sep='\t',index=None)

class_four_df=df.drop(remained)
print(class_four_df.shape[0])
class_four_df.to_csv(sys.argv[3],sep='\t',index=None)


    