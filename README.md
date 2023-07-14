# PTM-TBA

A toolkit for extracting PTM information and verifying PTM information with top-down MS, bottom-up MS, and UniProt annotations

## PTM-TBA tutorials
### Download

### Overview
* Preprocess mass shifts from peptide identifications
* Preprocess mass shifts from proteoforms
* Verify PTMs/mass shifts from proteoforms using peptides
* Extract PTM annotations from knowledge bases (UniProt and dbPTM)
* Verify PTMs/mass shifts from proteoforms using annotations

### 1. Preprocess mass shifts from peptide identifications
#### 1.1 Database search for bottom-up spectra
Requirements
* .mzML file of bottom-up spectra
* protein database in fasta format
* installed FragPipe

Perform the open search mode following the instructions [here](https://fragpipe.nesvilab.org/docs/tutorial_open.html).
#### 1.2 Extract mass shifts from peptide identifications
Input parameter:
* .tsv file containing identifications with mass shift information, e.g. psm.tsv in the output files of MS-Fragger
* output file name

Output:
* .tsv file containing mass shift information, including mass, possible range of localization, and annotation

Run the command:
```sh
python3 extract_ptm_from_psm.py psm.tsv psm_with_ms.tsv
```
#### 1.3 Remove duplicated mass shifts with greedy algorithm
Input parameter:
* .tsv file containing mass shift information, e.g. psm_with_ms.tsv
* error tolerance (Da), e.g. 0.1
* output file name

Output:
* .tsv file containing mass shift information without duplicates

Run the command:
```sh
python3 remove_dup_mass_shift_bu.py psm_with_ms.tsv 0.1 psm_with_ms_no_dup.tsv
```
#### (Optional) 1.4 Substitute FragPipe with MetaMorpheus 

### 2. Preprocess mass shifts from proteoform identifications
#### 2.1 Database search for top-down spectra
Requirements
* .mzML file of top-down spectra
* protein database in fasta format
* installed TopPIC Suite

Perform the open search following the instructions [here](https://www.toppic.org/software/toppic/tutorial.html).
#### 2.2 Get proteoforms with mass shifts
Input parameter:
* .tsv file containing identifications with mass shift information, e.g. toppic_proteoform_single.tsv in the output files of TopPIC
* output file name

Output:
* .tsv file containing proteoforms with mass shifts

Run the command:
```sh
python3 get_proteoforms_with_ms.py toppic_proteoform_single.tsv proteoform_with_ms.tsv
```
#### (Optional) 2.3 Remove proteoforms from histone proteins
Input parameter:
* .tsv file containing proteoforms with mass shifts
* output file name

Output:
* .tsv file containing proteoforms without those from histone proteins

Run the command:
```sh
python3 remove_histone.py proteoform_with_ms.tsv proteoform_with_ms_no_histone.tsv
```
#### 2.4 Remove duplicated mass shifts and split mass shifts into N-terminal acetylation and other PTMs
Input parameter:
* .tsv file containing proteoforms with mass shifts
* error tolerance (Da)
* output file name for N-terminal acetylation
* output file name for other PTMs

Output:
* .tsv file containing N-terminal acetylation
* .tsv file containing other PTMs

Run the command:
```sh
python3 remove_dup_mass_shift_top_down.py proteoform_with_ms_no_histone.tsv 0.1 mass_shift_n_term.tsv other_mass_shift.tsv 
```


### 1. Quick Start
##### 1.1  Run the examples of SW480 data sets in paper. 

    sh  PTM-TBA_one_step.sh -t 0.1
    
##### 1.2  Run with your own data sets and error tolerance t for your purpose.

    sh PTM-TBA_one_step.sh -td_input /path/to/proteoform.tsv -bu_input /path/to/peptides.tsv -anno uniprot_preferred_gene_modified_residue.tsv -t 0.1
    
### 2. Step-by-step pipeline
##### 2.1 Preprocessing mass shifts identified from top-down MS data
```sh
# Extract mass shifts from proteoform identifications
python3 extract_ptm_from_proteoform.py proteoform.tsv proteoform_with_ms.tsv
# (Optional) Remove proteoforms from histone proteins
python3 remove_histone.py proteoform_with_ms.tsv proteoform_with_ms_no_histone.tsv
# Remove duplicated mass shifts and split mass shifts with N-terminal acetylation and other PTMs, input with error tolerance 0.1 Da
python3 remove_dup_mass_shift_top_down.py proteoform_with_ms.tsv 0.1 other_mass_shift.tsv mass_shift_n_term.tsv
# Extract mass shift information including mass and range 
python3 extract_ptm_to_df.py mass_shift.tsv mass_shift_with_info.tsv
# Assign PTM type to mass shifts from high-frequency ones with a user-defined error tolerance (0.1 Da)
python3 find_confident_ptm.py mass_shift_with_info.tsv mass_shift_assigned_with_high_frequency.tsv mass_shift_not_identified.tsv 0.1

```
##### 2.2 Preprocessing mass shifts idenfied from bottom-up MS data
```sh
# Extract mass shifts from peptide identifications
python3 extract_ptm_from_psm.py psm.tsv psm_with_ms.tsv
# Remove duplicated mass shifts with greedy algorithm
python3 remove_dup_mass_shift_bu.py psm_with_ms.tsv 0.1 psm_with_ms_no_dup.tsv

```
##### 2.3  Verifying mass shifts from top-down MS with mass shifts from bottom-up MS with a user-defined error tolerance
```sh
python3 find_open_search_evidence.py input1=mass_shift_top_down.tsv input2=mass_shift_bottom_up.tsv t=0.1 output=mass_shift_td_matched_with_bu.tsv n-term-mode=0

```
##### 2.4 Verifying mass shifts from top-down MS with UniProt annotations with a user-defined error tolerance
```sh
python3 find_uniprot_evidence.py input1=mass_shift_top-down.tsv input2=uniprot_preferred_gene_modified_residue.tsv t=0.1 output=mass_shift_td_matched_with_anno.tsv
```



