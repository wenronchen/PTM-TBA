#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 12:34:28 2023

@author: wenrchen
"""


###Extract PTM information for output from MSPathFinder

import pandas as pd
import sys

df=pd.read_csv(sys.argv[1],sep='\t')
print(df.shape[0])
df_ms=pd.DataFrame(columns=df.columns)

df=df.sort_values(by=["Scan"],ignore_index=True)


##remain unique ID for the same scan
drop_list=[]
scan_dict={}
for i in range(df.shape[0]):
    scan=df.iloc[i]['Scan']
    if(scan not in scan_dict.keys()):
        scan_dict[scan]=(df.iloc[i]['EValue'],i)
    else:
        if(df.iloc[i]['EValue']<scan_dict[scan][0]):
            drop_list.append(scan_dict[scan][1])
            scan_dict[scan]=(df.iloc[i]['EValue'],i)
        else:
            drop_list.append(i)
print(len(drop_list))
df=df.drop(drop_list)

df=df.fillna(0)
df=df[df['Modifications']!=0]
print(df.shape[0])

df=df.reset_index()


ptm_type=[]
ptm_pos=[]
ptm_global_pos=[]

acet_cnt,phospho_cnt,oxi_cnt=0,0,0
for i in range(df.shape[0]):
    if(float(df.iloc[i]['EValue'])<=0.01):
        mods=df.iloc[i]['Modifications'].split(",")
        for mod in mods:
            mod_info=mod.split(" ")
            #print(mod_info)
            if(mod_info[0]=="Acetyl"):
                acet_cnt+=1
            if(mod_info[0]=="Oxidation"):
                oxi_cnt+=1
            if(mod_info[0]=="Phospho"):
                phospho_cnt+=1
            
            if(mod_info[0].find("Carbam")==-1):
                ptm_type.append(mod_info[0])
                ptm_pos.append(int(mod_info[1]))
                ptm_global_pos.append(int(mod_info[1])+df.iloc[i]['Start']-1)
                
                df_ms=df_ms.append(df.loc[i],ignore_index=True)
            
print("Acetyl count:" +str(acet_cnt))
print("Oxidation count:" +str(oxi_cnt))
print("Phospho count:" +str(phospho_cnt))
print(len(ptm_type))

df_ms.insert(2,"Mod type",ptm_type)
df_ms.insert(3,"Mod pos",ptm_pos)
df_ms.insert(4,"Mod global pos",ptm_global_pos)

df_ms.insert(0,"File name",sys.argv[3])

df_ms.to_csv(sys.argv[2],sep='\t',index=None)
        
        
        
    
        
    
        
        
    


