#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 09:59:37 2022

@author: wenrchen
"""
##find matched peptide and modifications

import pandas as pd
import sys



def extract_gene(tmp):
    gene_name=tmp.split("|")[-2]
    return gene_name

def extract_uniprot_gene(tmp):
    flag=tmp.find("GN=")
    tmp=tmp[flag+3:]
    flag=tmp.find(" ")
    
    return tmp[:flag]

def Is_range_overlap(range1,range2):
    return (range1.start<=range2.stop) and (range1.stop>=range2.start)

def extract_ptm(tmp):
    while(tmp.find("[")!=-1):
        flag1=tmp.find("[")
        flag2=tmp.find("]")
        tmp_ptm=tmp[flag1+1:flag2]
        #print(tmp_ptm)
        
        tmp=tmp[:flag1]+tmp[flag2+1:]
        
        if(tmp_ptm.find('Carbamidomethylation')==-1 and tmp_ptm.find('Acetyl')==-1):
            
            return float(tmp_ptm)

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

def compare_mass(m1, m2, t):##determine if the m1-m2, m1-m2+1.00235, m1-m2-1.00235 satisfy the error tolerance t
    if(abs(m1-m2)<=t):
        return True
    elif(abs(m1-m2-1.00235)<=t):
        return True
    elif(abs(m1-m2+1.00235)<=t):
        return True
    else:
        return False

def remove_dup_mod(evidence):
    evidence_dict={}
    for e in evidence:
        if(e[0] not in evidence_dict.keys()):
           evidence_dict[e[0]]=e
        else:
            if(evidence_dict[e[0]][4]>e[4]):
                evidence_dict[e[0]]=e
    return evidence_dict

dfname=sys.argv[1]        
if(dfname.find("xlsx")==-1):
   df=pd.read_csv(dfname,sep='\t') ##proteoforms with PTM df
else:
    df=pd.read_excel(dfname,sheet_name=0)

df_peptide=pd.read_csv(sys.argv[2],sep='\t') ##peptides with modifications and position

# print(df_peptide.iloc[0]['Is Unique'])
# df_peptide=df_peptide[df_peptide['Is Unique']==True]
# print(df_peptide.shape[0])

error_tolerance=float(sys.argv[3])
output=sys.argv[4]
n_term=int(sys.argv[5])   
 
df_peptide=df_peptide.fillna(0) 

peptide_dict={}

for i in range(df_peptide.shape[0]):
    p=df_peptide.iloc[i]['Protein']
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
    first_residue=df.iloc[i]['First residue']
    last_residue=df.iloc[i]['Last residue']
    
    seq=df.iloc[i]['Proteoform']
    
    ptm=extract_ptm(seq)
    ptm_range=get_ptm_range(seq,first_residue)
    if(n_term==1):
        ptm=42.0106
        ptm_range=range(first_residue,first_residue)
    
    protein_id=df.iloc[i]['Protein accession']
    
    evidence_tmp=[]
    valid_tmp=[]
    if(protein_id in peptide_dict.keys()):
        for peptide_id in peptide_dict[protein_id]:
            start_pos=int(df_peptide.iloc[peptide_id]['Protein Start'])
            end_pos=int(df_peptide.iloc[peptide_id]['Protein End'])
            
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
                ptm_mod.append(df_peptide.iloc[peptide_id]['Mod type'])
                ptm_mass.append(mod_mass)
                    
                
                
            
            
            
            
            # if(df_peptide.iloc[peptide_id]['Assigned Modifications']!=0):
                
            #     assign_mod=df_peptide.iloc[peptide_id]['Assigned Modifications'].split(",")
            #     for m in assign_mod:
            #         if(m.find("M(15.9949)")!=-1):
            #             if(start_pos>=first_residue and end_pos<=last_residue):
            #                 mod_seq=df_peptide.iloc[peptide_id]['Peptide']
            #                 e=df_peptide.iloc[peptide_id]['Expectation']
            #                 evidence_tmp.append((mod_seq,start_pos,end_pos,15.9949,e,m))
            #                 break
            #         # if(m.find("(0.9840)")!=-1):
            #         #     if(start_pos>=first_residue and end_pos<=last_residue):
            #         #         mod_seq=df_peptide.iloc[peptide_id]['Peptide']
            #         #         e=df_peptide.iloc[peptide_id]['Expectation']
            #         #         evidence_tmp.append((mod_seq,start_pos,end_pos,0.9840,e,m))
            #         #         break
            #         if(m.find("(42.0106)")!=-1):
            #             if(start_pos>=first_residue and end_pos<=last_residue):
            #                 mod_seq=df_peptide.iloc[peptide_id]['Peptide']
            #                 e=df_peptide.iloc[peptide_id]['Expectation']
            #                 evidence_tmp.append((mod_seq,start_pos,end_pos,42.0106,e,m))
            #                 break
                    # if(m.find("C(14.0157)")!=-1):
                    #     if(start_pos>=first_residue and end_pos<=last_residue):
                    #         if(df_peptide.iloc[peptide_id]['MSFragger Localization']!=0):
                    #             mod_seq=df_peptide.iloc[peptide_id]['MSFragger Localization']
                    #         else:
                    #             mod_seq=df_peptide.iloc[peptide_id]['Peptide']
                    #         e=df_peptide.iloc[peptide_id]['Expectation']
                        
                    #         evidence_tmp.append((mod_seq,start_pos,end_pos,14.0157,e,m))
                    #         break
            # if(df_peptide.iloc[peptide_id]['Observed Modifications']!=0):
            #     if(start_pos>=first_residue and end_pos<=last_residue):
            #         mod_seq=df_peptide.iloc[peptide_id]['MSFragger Localization']
            #         if(mod_seq==0):
            #             mod_seq=df_peptide.iloc[peptide_id]['Peptide']
            #         e=df_peptide.iloc[peptide_id]['Expectation']
            #         delta=df_peptide.iloc[peptide_id]['Delta Mass']
            #         mods=df_peptide.iloc[peptide_id]['Observed Modifications']
            #         evidence_tmp.append((mod_seq,start_pos,end_pos,delta,e,mods))
            # elif(df_peptide.iloc[peptide_id]['MSFragger Localization']!=0):
            #     if(start_pos>=first_residue and end_pos<=last_residue):
            #         mod_seq=df_peptide.iloc[peptide_id]['MSFragger Localization']
            #         e=df_peptide.iloc[peptide_id]['Expectation']
            #         delta=df_peptide.iloc[peptide_id]['Delta Mass']
            #         valid_tmp.append((mod_seq,start_pos,end_pos,delta,e))
                
    # evidence_tmp=list(set(evidence_tmp))   
        
    # if(len(evidence_tmp)!=0):
    #     proteoform_cnt+=1
    #     evidence_dict=remove_dup_mod(evidence_tmp)
        
        # if(abs(ptm+57)<error_tolerance):
        #     ptm_evidence_df=ptm_evidence_df.append(df.loc[i],ignore_index=True)
            
        #     match_list.append("C57 Match")
            
        #     ptm_evidence.append("()")
        #     ptm_mass.append(-57)
        #     ptm_e.append(0)
        #     ptm_mod.append("Carbamidomethylation")
            
        #     ptm_value.append(ptm)
        # else:
            
            
    #     for k in evidence_dict.keys():
    #         ptm_evidence_df=ptm_evidence_df.append(df.loc[i],ignore_index=True)
    #         tmp=evidence_dict[k]
            
    #         ptm_evidence.append((tmp[0],tmp[1],tmp[2]))
    #         ptm_mass.append(tmp[3])
    #         ptm_e.append(tmp[4])
    #         ptm_mod.append(tmp[5])
    #         if(compare_mass(tmp[3],ptm,error_tolerance)):
    #             match_list.append("Match")
    #             match_proteoform.append(i)
    #             #print(i)
    #             match_cnt+=1
    #         else:
    #             match_list.append("No Match")
    #         ptm_value.append(ptm)
    # if(i not in match_proteoform):
    #     if(len(valid_tmp)!=0):
    #         valid_dict=remove_dup_mod(valid_tmp)
    #         for k in valid_dict.keys():
    #             ptm_evidence_df=ptm_evidence_df.append(df.loc[i],ignore_index=True)
    #             tmp=valid_dict[k]
                
    #             ptm_evidence.append((tmp[0],tmp[1],tmp[2]))
    #             ptm_mass.append(tmp[3])
    #             ptm_e.append(tmp[4])
    #             ptm_mod.append("None")
    #             if(compare_mass(tmp[3],ptm,error_tolerance)):
    #                 match_list.append("Delta Mass Match")
    #                 match_proteoform.append(i)
    #                 #print(i)
    #                 match_cnt+=1
    #             else:
    #                 match_list.append("No Match")
    #             ptm_value.append(ptm)
            
print(ptm_evidence_df.shape[0],len(match_list))
print(len(list(set(match_proteoform))))
ptm_evidence_df.insert(20,"PTM evidence",ptm_evidence)
ptm_evidence_df.insert(21,"Delta mass",ptm_mass)
ptm_evidence_df.insert(22,"Peptide E-value",ptm_e)
ptm_evidence_df.insert(23,"Modification",ptm_mod)
ptm_evidence_df.insert(24,"Mass shift Match",match_list)
#ptm_evidence_df.insert(25,"Mass Shift",ptm_value)
ptm_evidence_df.to_csv(output,sep='\t',index=None)