#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 09:59:37 2022

@author: wenrchen
"""
##find matched peptide and modifications

import pandas as pd
import sys

import time
start_time = time.time()

def extract_gene(tmp):
    gene_name=tmp.split("|")[-2]
    return gene_name

# def extract_uniprot_gene(tmp):
#     flag=tmp.find("GN=")
#     tmp=tmp[flag+3:]
#     flag=tmp.find(" ")
    
#     return tmp[:flag]

def Is_range_overlap(range1,range2):
    return (range1.start<=range2.stop) and (range1.stop>=range2.start)

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

    return range(flag1+first,flag2+first)

def compare_mass(m1, m2, t):##determine if the m1-m2, m1-m2+1.00235, m1-m2-1.00235 satisfy the error tolerance t
    if(abs(m1-m2)<=t):
        return True
    elif(abs(m1-m2-1.00235)<=t):
        return True
    elif(abs(m1-m2+1.00235)<=t):
        return True
    else:
        return False

# def remove_dup_mod(evidence):
#     evidence_dict={}
#     for e in evidence:
#         if(e[0] not in evidence_dict.keys()):
#            evidence_dict[e[0]]=e
#         else:
#             if(evidence_dict[e[0]][4]>e[4]):
#                 evidence_dict[e[0]]=e
#     return evidence_dict

dfname=sys.argv[1]        
if(dfname.find("xlsx")==-1):
   df=pd.read_csv(dfname,sep='\t') ##proteoforms with PTM df
else:
    df=pd.read_excel(dfname,sheet_name=0)
    
ptm_dict={"Phospho":79.9663,"Acetyl":42.0106,"Oxidation":15.9949}

df_peptide=pd.read_csv(sys.argv[2],sep='\t') ##peptides with modifications and position

error_tolerance=float(sys.argv[3])
output=sys.argv[4]
n_term=int(sys.argv[5])   
 
df_peptide=df_peptide.fillna(0) 

peptide_dict={}

for i in range(df_peptide.shape[0]):
    p=df_peptide.iloc[i]['Protein'].split("|")[0]
    #p=df_peptide.iloc[i]['Protein Accession'].split("|")[0] ##for metamoepheus
    if(p not in peptide_dict.keys()):
        peptide_dict[p]=[i]
    else:
        peptide_dict[p].append(i)           


ptm_evidence_df=pd.DataFrame(columns=df.columns)
ptm_evidence=[]
ptm_mass=[]
ptm_e=[]
ptm_mod=[]
match_list=[]
ptm_value=[]
match_cnt=0
proteoform_cnt=0
match_proteoform=[]
for i in range(df.shape[0]):
    #first_residue=df.iloc[i]['First residue']
    
    # n_term_flag=0
    # if(df.iloc[i]['Protein N-terminal form'].find("ACETYLATION")!=-1):
    #     n_term_flag=1
    # variable_flag=0
    # if(df.iloc[i]['#variable PTMs']==1):
    #     variable_flag=1
    
    # seq=df.iloc[i]['Proteoform']
    
    # ptm=extract_ptm(seq,n_term_flag)
    # ptm_range=get_ptm_range(seq,first_residue,n_term_flag,variable_flag)
    # if(n_term==1):
    #     ptm=42.0106
    #     ptm_range=range(first_residue,first_residue)
    
    # #protein_id=df.iloc[i]['Protein accession'].split("|")[1]
    # protein_id=df.iloc[i]['Protein accession']
    ptm=ptm_dict[df.iloc[i]['Modifications'].split(" ")[0]]
    ptm_range=range(df.iloc[i]['Mod global pos'],df.iloc[i]['Mod global pos'])
    protein_id=df.iloc[i]['ProteinName'].split("|")[0]
    
    
    evidence_tmp=[]
    valid_tmp=[]
    if(protein_id in peptide_dict.keys()):
        for peptide_id in peptide_dict[protein_id]:

            mod_start=df_peptide.iloc[peptide_id]["Mod start"]
            mod_end=df_peptide.iloc[peptide_id]["Mod end"]
            mod_range=range(mod_start,mod_end)
            mod_mass=df_peptide.iloc[peptide_id]['Mass shift']
            
            if(Is_range_overlap(ptm_range, mod_range)):
                ptm_evidence_df=ptm_evidence_df.append(df.loc[i],ignore_index=True)
                if(compare_mass(ptm,mod_mass,error_tolerance)):
                    match_list.append("Match")
                    match_proteoform.append(i)
                else:
                    match_list.append("Not Match")
                
                mod_seq=str(df_peptide.iloc[peptide_id]['MSFragger Localization'])
                if(mod_seq!="0"):
                    ptm_evidence.append((mod_seq,mod_range))
                else:
                    seq_tmp=df_peptide.iloc[peptide_id]['Peptide']
                    ptm_evidence.append((seq_tmp,mod_range))
                
                ptm_e.append(df_peptide.iloc[peptide_id]['Expectation'])
                
                # mod_seq=str(df_peptide.iloc[peptide_id]['Essential Sequence'])
                # ptm_evidence.append((mod_seq,mod_range))
                # ptm_e.append(df_peptide.iloc[peptide_id]['PEP_QValue']) ##for metamorpheus
                
                ptm_mod.append(df_peptide.iloc[peptide_id]['Mod type'])
                ptm_mass.append(mod_mass)

            
print(ptm_evidence_df.shape[0],len(match_list))
print(len(list(set(match_proteoform))))
ptm_evidence_df.insert(20,"PTM evidence",ptm_evidence)
ptm_evidence_df.insert(21,"Delta mass",ptm_mass)
ptm_evidence_df.insert(22,"Peptide E-value",ptm_e)
ptm_evidence_df.insert(23,"Modification",ptm_mod)
ptm_evidence_df.insert(24,"Mass shift Match",match_list)
#ptm_evidence_df.insert(25,"Mass Shift",ptm_value)
ptm_evidence_df.to_csv(output,sep='\t',index=None)

print("--- %s seconds ---" % (time.time() - start_time))