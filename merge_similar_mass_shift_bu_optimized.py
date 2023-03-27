#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 15:49:11 2023

@author: wenrchen
"""
##remove duplicated mass shifts from those extracted from PSM table (optimized)
import pandas as pd
import sys

def Is_range_subset(range1,range2):
    
    if(range1==range2):
        return range1==range2
    else:
        if((range1.stop-range1.start)<(range2.stop-range2.start)):
        
            return ((range1.start>=range2.start) and (range1.stop<=range2.stop))
        else:
            return ((range2.start>= range1.start) and (range2.stop <= range1.stop))
        


df=pd.read_csv(sys.argv[1],sep='\t')

df=df.sort_values(by="Expectation", ignore_index=True)
error_tolerance=float(sys.argv[2])
output=sys.argv[3]

protein_ms_dict={} ## key=protein value=mass shifts on this protein
for i in range(df.shape[0]):
    protein_id=df.iloc[i]['Protein ID']
    if(protein_id not in protein_ms_dict.keys()):
        protein_ms_dict[protein_id]=[i]
    else:
        protein_ms_dict[protein_id].append(i)

drop_list=[]
for key,value in protein_ms_dict.items():
    ms_list=[]
    ms_dict={}
    for v in value:
        ms=df.iloc[v]['Mass shift']
        flag=0
        for m in ms_dict.keys():
            if(abs(m-ms)<=0.1):
                ms_dict[m].append(v)
                flag=1
        if(flag==0):
            ms_dict[ms]=[v]
    
    for m in ms_dict.keys():
        ms_df=df.loc[ms_dict[m],:] ##select rows with specified protein and similar mass shifts
        ms_df.insert(0,"original_index",ms_dict[m])
        ms_df=ms_df.sort_values(by=["Mod start"],ignore_index=True)
        
        i=0
        flag=1
        while (i<ms_df.shape[0] and flag<ms_df.shape[0]):
            while(i<ms_df.shape[0] and flag<ms_df.shape[0] and ms_df.iloc[flag]['Mod start']<=ms_df.iloc[i]['Mod end']):
                if(ms_df.iloc[flag]['Mod end']<=ms_df.iloc[i]['Mod end']):
                    drop_list.append(ms_df.iloc[i]["original_index"])
                    i=flag
                    flag=i+1
                else:
                    # i_len=ms_df.iloc[i]['Mod end']-ms_df.iloc[i]['Mod start']
                    # flag_len=ms_df.iloc[flag]['Mod end']-ms_df.iloc[flag]['Mod start']
                    # if(i_len<=flag_len):
                    #     drop_list.append(ms_df.iloc[i]["original_index"])
                    #     flag+=1
                    # else:
                    #     drop_list.append(ms_df.iloc[i]["original_index"])
                    #     i=flag
                    #     flag=i+1
                    drop_list.append(ms_df.iloc[flag]["original_index"])
                    #i=flag
                    flag=flag+1
                    
            if(i<ms_df.shape[0] and flag<ms_df.shape[0] and ms_df.iloc[flag]['Mod start']>ms_df.iloc[i]['Mod end']):
                i=i+1
                flag=i+1
            

df=df.drop(drop_list)
print(df.shape[0])
df.to_csv(output,sep='\t',index=None)