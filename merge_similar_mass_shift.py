#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 10:20:14 2022

@author: wenrchen
"""
##merge proteoforms with similar sequence and similar mass shifts.

import pandas as pd
import sys

def extract_ptm(tmp):
    while(tmp.find("[")!=-1):
        flag1=tmp.find("[")
        flag2=tmp.find("]")
        tmp_ptm=tmp[flag1+1:flag2]
        #print(tmp_ptm)
        
        tmp=tmp[:flag1]+tmp[flag2+1:]
        
        if(tmp_ptm.find('Carbamidomethylation')==-1 and tmp_ptm.find('Acetyl')==-1):
            if(tmp_ptm.find(".")!=-1):
                return float(tmp_ptm)
            else:
                return tmp_ptm

# def is_float(s):
#     try:
#         float(s)
#         return True
#     except ValueError:
#         return False
    
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
    
    #print(tmp)
    flag1=tmp.find("(")
    flag2=tmp.find(")")
    
    return range(flag1+first,flag2+first)

def Is_range_subset(range1,range2):
    
    if(range1==range2):
        return range1==range2
    else:
        if((range1.stop-range1.start)<(range2.stop-range2.start)):
        
            return ((range1.start>=range2.start) and (range1.stop<=range2.stop))
        else:
            return ((range2.start>= range1.start) and (range2.stop <= range1.stop))

def Is_range_overlap(range1,range2):
    return (range1.start<=range2.stop) and (range1.stop>=range2.start)

def compare_mass(m1, m2, t):##determine if the m1-m2, m1-m2+1.00235, m1-m2-1.00235 satisfy the error tolerance t
    if(abs(m1-m2)<=t):
        return True
    elif(abs(m1-m2-1.00235)<=t):
        return True
    elif(abs(m1-m2+1.00235)<=t):
        return True
    else:
        return False

df=pd.read_csv(sys.argv[1],sep='\t') ##proteoforms with PTM df
error_tolerance=float(sys.argv[2]) ##1.00235 Da
output=sys.argv[3]
output_n_term=sys.argv[4]

n_term_dict={} 
proteoform_dict={}
n_term_cnt=0
drop_list=[]
drop_list_n_term=[]
for i in range(df.shape[0]):
    if(df.iloc[i]['#unexpected modifications']!=0):
        protein_id=df.iloc[i]['Protein accession']
        if(protein_id not in proteoform_dict.keys()):
            proteoform_dict[protein_id]=[i]
        else:
            proteoform_dict[protein_id].append(i)
    else:
        drop_list.append(i)
    if(df.iloc[i]['Proteoform'].find("[Acetyl]")!=-1):
        n_term_cnt+=1
        protein_id=df.iloc[i]['Protein accession']
        if(protein_id not in n_term_dict.keys()):
            n_term_dict[protein_id]=[i]
        else:
            n_term_dict[protein_id].append(i)
    else:
        drop_list_n_term.append(i)


for k,v in proteoform_dict.items(): ##remove dup unexpected mass shifts
    for i in range(len(v)-1):
        first1=df.iloc[v[i]]['First residue']
        last1=df.iloc[v[i]]['Last residue']
        seq1=df.iloc[v[i]]['Proteoform']
        ptm1=extract_ptm(seq1)
        range1=get_ptm_range(seq1,first1)
        #print(range1)
        for j in range(i+1,len(v)):
            first2=df.iloc[v[j]]['First residue']
            last2=df.iloc[v[j]]['Last residue']
            #if(first1==first2 and last1==last2):
            seq2=df.iloc[v[j]]['Proteoform']
            ptm2=extract_ptm(seq2)
            range2=get_ptm_range(seq2,first2)
            if(str(ptm1).find('.')!=-1 and str(ptm2).find('.')!=-1):
                if(compare_mass(ptm1,ptm2,error_tolerance)):
                    if(Is_range_overlap(range1,range2)):
                        if(df.iloc[v[i]]['E-value']<=df.iloc[v[j]]['E-value']):
                            drop_list.append(v[j])
                        else:
                            drop_list.append(v[i])

for k,v in n_term_dict.items(): ##remove dup n-term acetylation
    for i in range(len(v)-1):
        first1=df.iloc[v[i]]['First residue']
        
        for j in range(i+1,len(v)):
            first2=df.iloc[v[j]]['First residue']
            if(first1==first2):
                if(df.iloc[v[i]]['E-value']<=df.iloc[v[j]]['E-value']):
                    drop_list_n_term.append(v[j])
                else:
                    drop_list_n_term.append(v[i])
#print(drop_list)
print(n_term_cnt)

df_n_term=df.drop(drop_list_n_term)
print(df_n_term.shape[0])
df_n_term.to_csv(output_n_term,sep='\t',index=None)

df=df.drop(drop_list)
print(df.shape[0])
df.to_csv(output,sep='\t',index=None)







                    
                
            
        