#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 13:25:06 2022

@author: wenrchen
"""

#find confident oxidations 
import pandas as pd
import sys

fix_mod_form1="(C)[Carbamidomethylation]"
fix_mod_form2="C[Carbamidomethylation]"

def get_ptm_seq(tmp,n_term_flag,variable_flag):
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


def read_ptm_list(filename):
    with open(filename, 'r') as file:
    # Read all the lines from the file into a list
        lines = file.readlines()
        print(list(lines[0].split(",")[2]))

    # Iterate through the list and print each line
        out_dict={}
        ptm_mass_dict={}
    
        for line in lines:
            ptm=line.split(",")
            
            out_dict[ptm[0]]=list(ptm[2])
            ptm_mass_dict[ptm[0]]=float(ptm[1])
    file.close()
    
    return out_dict,ptm_mass_dict

df=pd.read_csv(sys.argv[1],sep='\t')

df=df.fillna(0)

output_common=sys.argv[2]
output_class_four=sys.argv[3]
error_t=float(sys.argv[4])
ptm_list_file=sys.argv[5]
##PTMs with high frequency from bottom-up MS, 44.008456:["S"],31.989829:["P","R","K","M","F","W","Y","C","E","I","L","V"]
# potential_aa_dict={21.981943:["D","E"],37.955882:["D","E"],52.911464:["D","E"],53.919289:["D","E"],\
#                    43.005814:["K","R","C","M","S","T","Y"],79.966331:["S","T","Y"],183.035399:["H","K","S","Y"],\
#                        14.01565:["C","I","L","K","R","H","D","E","N","Q","S","T"],28.0313:["K","N","R","P"],\
#                           42.010565:["K","C","S","T","Y","H","R"], 44.008456:["S"],\
#                               15.994915:["M","W","H","D","K","N","P","Y","R","C","G","U","E","I","L","Q","S","T","V"]\
#                                   }
potential_aa_dict, ptm_mass_dict=read_ptm_list(ptm_list_file)
ptm_cnt_dict={}
remained=[]
cnt=0
oxi_cnt,ace_cnt,pho_cnt=0,0,0
for i in range(df.shape[0]):
    if(df.iloc[i]['Mass shift']!=0):
        
        n_term_flag=0
        if(df.iloc[i]['Protein N-terminal form'].find("ACETYLATION")!=-1):
            n_term_flag=1
        variable_flag=0
        if(df.iloc[i]['#variable PTMs']==1):
            variable_flag=1
        
        for ptm in potential_aa_dict.keys():
    
        # elif(ptm_value=="Acetyl"):
        #     ace_cnt+=1
        #     remained.append(i)
            ptm_value=ptm_mass_dict[ptm]
            if (df.iloc[i]['Mass shift'] in ptm_mass_dict.values()):
                mass_shift=ptm_mass_dict[df.iloc[i]['Mass shift']]
            else:
                mass_shift=float(df.iloc[i]['Mass shift'])
            if(compare_mass(mass_shift,ptm_value,error_t)):
                
                if(ptm_value==15.994915):
                    oxi_cnt+=1
                if(ptm_value==42.010565):
                    ace_cnt+=1
                if(ptm_value==79.966331):
                    pho_cnt+=1
                potential_aa=potential_aa_dict[ptm]
                cnt+=1
                ptm_range=df.iloc[i]['PTM Range'].split(",")
                ptm_first=int(ptm_range[0][1:])
                ptm_last=int(ptm_range[1][:-1])
                
                ptm_seq=get_ptm_seq(df.iloc[i]['Proteoform'],n_term_flag,variable_flag)
                for aa in potential_aa:
                    if(aa in ptm_seq):
                        
                        if(i not in remained):
                            if(ptm_value not in ptm_cnt_dict.keys()):
                                ptm_cnt_dict[ptm_value]=1
                            else:
                                ptm_cnt_dict[ptm_value]+=1
                        remained.append(i)
                        break


print(cnt)
print("Acetyl count:" +str(ace_cnt))
print("Oxidation count:" +str(oxi_cnt))
print("Phospho count:" +str(pho_cnt))
remained=list(set(remained))
print(len(remained))
# print(ptm_cnt_dict)

output_df=df.loc[remained]
output_df.to_csv(output_common,sep='\t',index=None)


class_four_df=df.drop(remained)
print(class_four_df.shape[0])
class_four_df.to_csv(output_class_four,sep='\t',index=None)


    
