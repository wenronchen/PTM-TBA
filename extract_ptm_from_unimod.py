#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 19:17:34 2023

@author: wenrchen
"""
import pandas as pd
import sys


argv=sys.argv

content = []
with open(argv[1], 'r', encoding='UTF-8') as f:
    for line in f.readlines():
        content.append(line.rstrip('\n'))

f = open(argv[2], "w")
#print("Name","Mass","Residues","Position","UnimodID", file = f, sep="\t") 

ptm_name=[]
l=0
while(l<len(content)):
    if(content[l].find("[Term]")!=-1):
        l+=1
        if(content[l].find("id: UNIMOD:0")==-1):
            unimod_id=int(content[l].split(":")[2])
            l+=1
            name=content[l].split(":")[1][1:]
            
            residue="" 
            while(content[l].find("is_a: UNIMOD:0 ! unimod root node")==-1):
                l+=1
                if(content[l].find("delta_mono_mass")!=-1):
                    mass_flag=content[l].find("\"")
                    mass=float(content[l][mass_flag+1:-1])
                if(content[l].find("_site")!=-1):
                    site_flag=content[l].find("\"")
                    site=content[l][site_flag+1:-1]
                    if(site.find("term")==-1):
                        residue=residue+site
            if(name not in ptm_name):
                ptm_name.append(name)
                print(name, mass, residue, "any", unimod_id, file = f, sep=",")
                
    l+=1     
                
f.close()           
                

        
        
        
        

# ID_list = [ss for ss in content if "id: " in ss] # find the "id:"
# # print('Total IDs are: {} '.format(len(ID_list)))

# Name_list = [ss for ss in content if "name: " in ss] # find the "id:"
# # print('Total Names are: {} '.format(len(Name_list)))




# index_sp = []
# for line in range(len(content)):
#     if "id: " in content[line]:
#         index_sp.append(line)
# # print(index_sp)

# id_list = []
# for i in range(len(index_sp)):
#     id_content = content[index_sp[i]]
#     id_list.append(id_content)
    
# index_sp = []
# for line in range(len(content)):
#     if "name: " in content[line]:
#         index_sp.append(line)
        
# name_list = []
# for i in range(len(index_sp)):
#     name_content = content[index_sp[i]]
#     name_list.append(name_content)
    
# data = {'ID': id_list,
#        'Name': name_list}
# id_name_list_df = pd.DataFrame(data) 

# # delete redundant information and make it better
# for i in range(len(id_name_list_df)):
#     temp_name = id_name_list_df['Name'].iloc[i]
#     id_name_list_df['Name'].iloc[i] = temp_name[6:]
#     temp_id = id_name_list_df['ID'].iloc[i]
#     id_name_list_df['ID'].iloc[i] = temp_id[11:]
# id_name_list_df


# index_sp = []
# for line in range(len(content)):
#     if "xref: delta_mono_mass" in content[line]:
#         index_sp.append(line)
        
# # extract the mass
# mass_list = []
# for i in range(len(index_sp)):
#     mass_content = content[index_sp[i]]
#     mass_list.append(mass_content)

# data = {'Unimod Mass': mass_list}
# mass_list_df = pd.DataFrame(data) 

# # extract the values from the mass list
# for i in range(len(mass_list_df)):
#     temp_name = mass_list_df['Unimod Mass'].iloc[i]
#     mass_list_df['Unimod Mass'].iloc[i] = temp_name[23:-1]
    
# # set a dummy value for id=0, ex. = 9999
# df2 = pd.DataFrame(data= ['9999'], columns = mass_list_df.columns)

# # reset the index
# mass_list_df_new = pd.concat([df2,mass_list_df], axis = 0, ignore_index = True)  
# mass_list_df_new

# # add the mass value into the name list
# id_name_mass_list_df = pd.concat([id_name_list_df, mass_list_df_new], axis = 1) 
# id_name_mass_list_df

# save it to csv file
# id_name_mass_list_df.to_csv("C:/Users/kunza/Tulane_project/Data/Unimod_id_name_mass_results_all.csv",index=False)