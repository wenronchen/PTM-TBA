#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 10:10:05 2022

@author: wenrchen
"""
##find associated peptides for proteoforms and check modifications

import pandas as pd
import sys

df=pd.read_csv(sys.argv[1],sep='\t') ##proteoforms with PTM df
df_peptide=pd.read_csv(sys.argv[2],sep='\t') ##peptides with modifications and position
df_evidence=pd.read_csv(sys.argv[3],sep='\t') ##evidence (localization) for modifications
ptm=sys.argv[4]
ptm_mass=float(sys.argv[5])
error_tolerance=float(sys.argv[6])
output=sys.argv[7]

def extract_gene(tmp):
    gene_name=tmp.split("|")[-2]
    return gene_name

def extract_ptm(tmp):
    while(tmp.find("[")!=-1):
        flag1=tmp.find("[")
        flag2=tmp.find("]")
        tmp_ptm=tmp[flag1+1:flag2]
        #print(tmp_ptm)
        
        tmp=tmp[:flag1]+tmp[flag2+1:]
        
        if(tmp_ptm.find('Carbamidomethylation')==-1 and tmp_ptm.find('Acetyl')==-1):
            
            return float(tmp_ptm)
        

# def remove_dup_mod(evidence):
#     evidence_mod_dict={}
#     for e in evidence:
#         if(e[0] not in evidence_mod_dict.keys()):
#            evidence_mod_dict[e[0]]=e
#         else:
#             if(evidence_mod_dict[e[0]][4]>e[4]):
#                 evidence_mod_dict[e[0]]=e
#     return evidence_mod_dict

peptide_dict={}
for i in range(df_peptide.shape[0]):
    protein_ids=str(df_peptide.iloc[i]['Proteins']).split(";")
    for p in protein_ids:
        if(p not in peptide_dict.keys()):
            peptide_dict[p]=[i]
        else:
            peptide_dict[p].append(i)
                
evidence_dict={}
for i in range(df_evidence.shape[0]):
    if(df_evidence.iloc[i]['Modifications']!="Unmodified"):
        if(df_evidence.iloc[i]['Mod. peptide ID'] not in evidence_dict.keys()):
            evidence_dict[df_evidence.iloc[i]['Mod. peptide ID']]=[i]
        else:
            evidence_dict[df_evidence.iloc[i]['Mod. peptide ID']].append(i)



ptm_evidence_df=pd.DataFrame(columns=df.columns)
ptm_evidence=[]
match_list=[]
ptm_value=[]
for i in range(df.shape[0]):
    first_residue=df.iloc[i]['First residue']
    last_residue=df.iloc[i]['Last residue']
    
    seq=df.iloc[i]['Proteoform']
    ptm_mass_shift=extract_ptm(seq)

    
    protein_id=df.iloc[i]['Protein accession']
    evidence_tmp=[]
    if(protein_id in peptide_dict.keys()):
        for peptide_id in peptide_dict[protein_id]:
            if(df_peptide.iloc[peptide_id]['Modifications'].find(ptm)!=-1):
                #evidence_ids=df_peptide.iloc[peptide_id]['Evidence IDs'].split(";")[:-1]
                evidence_ids=evidence_dict[peptide_id]
                for e in evidence_ids:
                    if(df_evidence.iloc[int(e)]['Modifications'].find(ptm)!=-1):
                        start_pos=int(df_peptide.iloc[peptide_id]['Start position'])
                        end_pos=int(df_peptide.iloc[peptide_id]['End position'])
                        if(start_pos>=first_residue and end_pos<=last_residue):
                            modified_seq=df_evidence.iloc[int(e)]['Modified sequence']
                            evidence_tmp.append((modified_seq,start_pos,end_pos))
                            break
    if(len(evidence_tmp)!=0):
        ptm_evidence_df=ptm_evidence_df.append(df.loc[i],ignore_index=True)
        
        for evidence in evidence_tmp:
            if(abs(ptm_mass-ptm_mass_shift)<=error_tolerance):
                match_list.append(evidence)
                print(i)
                break
        if(len(match_list)==len(ptm_evidence)):
            match_list.append("No Match")

        ptm_value.append(ptm_mass_shift)
        ptm_evidence.append(evidence_tmp)

print(ptm_evidence_df.shape[0],len(ptm_evidence))

ptm_evidence_df.insert(20,"PTM evidence",ptm_evidence)
ptm_evidence_df.insert(21,"Mass Match",match_list)
ptm_evidence_df.insert(22,"Mass Shift",ptm_value)
ptm_evidence_df.to_csv(output,sep='\t',index=None)
