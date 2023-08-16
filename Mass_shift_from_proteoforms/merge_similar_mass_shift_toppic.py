#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 10:20:14 2022

@author: wenrchen
"""
##merge proteoforms with similar sequence and similar mass shifts.

import pandas as pd
import sys


fix_mod_form1="(C)[Carbamidomethylation]"
fix_mod_form2="C[Carbamidomethylation]"

def extract_ptm(tmp,n_term_flag):
    if(n_term_flag==1):
        flag=tmp.find("[42.0106]")
        if(tmp[flag-1]==")"):
            tmp=tmp[:flag-3]+tmp[flag-2]+tmp[flag+9:]
        else:
            tmp=tmp[:flag]+tmp[flag+9:]
    while(tmp.find("[")!=-1):
        flag1=tmp.find("[")
        flag2=tmp.find("]")
        tmp_ptm=tmp[flag1+1:flag2]

        tmp=tmp[:flag1]+tmp[flag2+1:]
        
        if(tmp_ptm.find('Carbamidomethylation')==-1):
            
            if(tmp_ptm.find(".")!=-1):
                return float(tmp_ptm)
            else:
                return tmp_ptm
  
def get_ptm_range(tmp,first,n_term_flag,variable_flag):
    #remove the letter before "." and the letter after "."

    #remove the letter between "[" and  "]"(include"[" and "]")
    
    if(n_term_flag==1):
        flag=tmp.find("[42.0106]")
        if(tmp[flag-1]==")"):
            tmp=tmp[:flag-3]+tmp[flag-2]+tmp[flag+9:]
        else:
            tmp=tmp[:flag]+tmp[flag+9:]
    
    while(tmp.find("[")!=-1):
        
        # if(tmp.find("[Acetyl]")!=-1):
        #     flag=tmp.find("[Acetyl]")
        #     tmp=tmp[:flag-3]+tmp[flag-2]+tmp[flag+len("[Acetyl]"):]
        
        
        while(tmp.find(fix_mod_form1)!=-1):
            flag=tmp.find(fix_mod_form1)
            tmp=tmp[:flag]+tmp[flag+1]+tmp[flag+len(fix_mod_form1):]
        
        while(tmp.find(fix_mod_form2)!=-1):
            flag=tmp.find(fix_mod_form2)
            tmp=tmp[:flag+1]+tmp[flag+len(fix_mod_form2):]
        
        flag1=tmp.find("[")
        flag2=tmp.find("]")
        tmp=tmp[:flag1]+tmp[flag2+1:]
        
        if(tmp.find("[")==-1):
            
            if(variable_flag==1):
                if(tmp[flag1-1]!=")"):
                    tmp=tmp[:flag1-1]+"("+tmp[flag1-1]+")"+tmp[flag1:]

    flag=tmp.find(".")
    tmp=tmp[flag+1:]
    flag=tmp.find(".")
    tmp=tmp[:flag]
    
    flag1=tmp.find("(")
    flag2=tmp.find(")")

    return range(flag1+first,flag2+first-2)


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
error_tolerance=float(sys.argv[2]) ##0.1 Da
output_n_term=sys.argv[3]
output=sys.argv[4]

#"[Delta:H(4)C(2)O(-1)S(1)]":44.0085
replace_table={"[Acetyl]":42.0106,"[Phospho]":79.9663, "[Oxidation]":15.9949, "[Methyl]":14.0157, "[Cation:Na]":21.9819,\
              "[Cation:K]":37.9559, "[Dimethyl]":28.0313, "[Carbamyl]":43.0058, "[AEBS]":183.0354, "[Cation:Fe[III]]":52.9115,\
                 "[Cation:Fe[II]]":53.9193, "[Delta:H(4)C(2)O(-1)S(1)]":44.0085,"[Dioxidation]":31.989829}

def replace_subseq(seq):
    for ptm in replace_table.keys():
        seq=seq.replace(ptm,"["+str(replace_table[ptm])+"]")
    
    # new_seq=seq.replace("[Phospho]", "[79.9663]")
    # new_seq=new_seq.replace("[Oxidation]","[15.9949]")
    
    return seq
    
for i in range(df.shape[0]):
    seq=df.iloc[i]['Proteoform']
    new_seq=replace_subseq(seq)
    
    df.at[i,"Proteoform"]=new_seq
    



n_term_dict={} 
proteoform_dict={}
n_term_cnt=0
drop_list=[]
drop_list_n_term=[]
for i in range(df.shape[0]):
    if((df.iloc[i]['#unexpected modifications']!=0) or (df.iloc[i]['#variable PTMs']==1)):
        protein_id=df.iloc[i]['Protein accession']
        if(protein_id not in proteoform_dict.keys()):
            proteoform_dict[protein_id]=[i]
        else:
            proteoform_dict[protein_id].append(i)
    else:
        drop_list.append(i)
    if(df.iloc[i]['Protein N-terminal form'].find("ACETYLATION")!=-1):
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
        
        n_term_flag=0
        if(df.iloc[v[i]]['Protein N-terminal form'].find("ACETYLATION")!=-1):
            n_term_flag=1
        variable_flag=0
        if(df.iloc[v[i]]['#variable PTMs']==1):
            variable_flag=1
        
        ptm1=extract_ptm(seq1,n_term_flag)
        range1=get_ptm_range(seq1,first1,n_term_flag,variable_flag)
        #print(range1)
        for j in range(i+1,len(v)):
            first2=df.iloc[v[j]]['First residue']
            last2=df.iloc[v[j]]['Last residue']
            #if(first1==first2 and last1==last2):
            seq2=df.iloc[v[j]]['Proteoform']
            
            
            n_term_flag=0
            if(df.iloc[v[j]]['Protein N-terminal form'].find("ACETYLATION")!=-1):
                n_term_flag=1
            variable_flag=0
            if(df.iloc[v[j]]['#variable PTMs']==1):
                variable_flag=1
            
            ptm2=extract_ptm(seq2,n_term_flag)
            range2=get_ptm_range(seq2,first2,n_term_flag,variable_flag)
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







                    
                
            
        