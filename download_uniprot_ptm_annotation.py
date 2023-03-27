#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 15:33:23 2022

@author: wenrchen
"""
#extract PTM annotation for specific genes from UniProt

import urllib.parse
import urllib.request
import sys

output=sys.argv[1]


# url = 'https://www.uniprot.org/uploadlists/'

# params = {
#     'from': 'ACC+ID',
#     'format': 'tab',
#     "organism": '9606',
#     "columns": "id,genes,feature(MODIFIED RESIDUE)",
#     "compress": "no",
#     "query": '',# NOTE: not sure if this one works
#     }

url = f"https://www.uniprot.org/uniprot/?query=organism:{9606}&columns=id,genes(PREFERRED),feature(MODIFIED%20RESIDUE)&format=tab"

# data = urllib.parse.urlencode(params)
# data = data.encode('utf-8')
req = urllib.request.Request(url)
with urllib.request.urlopen(req) as f:
    response = f.read()
    annot_filename = output
    with open(annot_filename, 'w') as as_file:
        decoded_response = response.decode("utf-8")
        if decoded_response[:5] == "Entry":
            as_file.write(decoded_response)
#print(response.decode('utf-8'))