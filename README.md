# PTM-TBA

A toolkit for extracting PTM information and verifying PTM information with top-down MS, bottom-up MS, and UniProt annotations

## PTM-TBA tutorials
### Download
You can download the Python code and test examples [here](https://github.com/wenronchen/PTM-TBA/archive/refs/heads/master.zip). 
### Overview
* [Bottom-up MS database search](https://github.com/wenronchen/PTM-TBA/blob/master/README.md#1-bottom-up-ms-database-search)
   - Database search using MS-Fragger
   - Database search using MetaMorpheus
   - Database search using MaxQuant
* Preprocess mass shifts from peptide identifications
* [Top-down MS database search](https://github.com/wenronchen/PTM-TBA/blob/master/README.md#2-preprocess-mass-shifts-from-proteoform-identifications)
    - Database search using TopPIC
    - Database search using MSPathFinder
* Preprocess mass shifts from proteoform identifications
* [Verification of mass shifts from proteoforms using mass shifts from peptides](https://github.com/wenronchen/PTM-TBA/blob/master/README.md#3-verifying-mass-shifts-from-proteoforms-with-mass-shifts-from-peptides)
    - MS-Fragger + TopPIC
    - MS-Fragger + MSPathFinder
    - MetaMorpheus + TopPIC
    - MaxQuant + TopPIC
* [PTM annotations extraction](https://github.com/wenronchen/PTM-TBA/blob/master/README.md#4-preprocess-ptm-annotations-from-knowledge-bases-uniprot-and-dbptm)
    - UniProt
    - dbPTM
* [Verification of mass shifts from proteoforms using PTM annotations](https://github.com/wenronchen/PTM-TBA/blob/master/README.md#5--verifying-mass-shifts-from-proteoforms-with-ptm-annotations)
   - UniProt + TopPIC
   - dbPTM + TopPIC
* [Extract PTM information from UNIMOD](https://github.com/wenronchen/PTM-TBA/blob/master/README.md#6-extract-ptm-information-from-unimod) 


### 1. Bottom-up MS database search
#### 1.1 Database search with MS-Fragger
Requirements
* .mzML file of bottom-up spectra
* protein database in fasta format
* installed FragPipe

Perform the open search mode following the instructions [here](https://fragpipe.nesvilab.org/docs/tutorial_open.html).

#### 1.2 Database search with MetaMorpheus
Download and install MetaMorpheus following the instructions [here](https://github.com/smith-chem-wisc/MetaMorpheus/wiki/Getting-Started#test-conda-installation-linux-macos-windows). Download the parameters for the five Tasks executed in MetaMorpheus.

Run MetaMorpheus via the command line:
```sh
metamorpheus -t Task1-SearchTaskconfig.toml Task2-CalibrateTaskconfig.toml Task3-SearchTaskconfig.toml Task4-GPTMDTaskconfig.toml Task5-SearchTaskconfig.toml -s test.raw -d human.fasta
```

#### 1.3 Database search with MaxQuant
Download and install MaxQuant following the instructions [here](https://www.maxquant.org/). Set the PTM in which you are interested as the variable PTM when you perform the database search using MaxQuant.

### 2. Preprocess mass shifts from peptide identifications 
#### 2.1 Extract mass shifts from peptide identifications (MS-Fragger and MetaMorpheus)
Input parameter:
* .tsv file containing identifications with mass shift information, e.g. psm.tsv in the output files of MS-Fragger
* search engine name, e.g. "msfragger"
* output file name

Output:
* .tsv file containing mass shift information, including mass, possible range of localization, and annotation

Run the command:
```sh
python3 extract_ptm_from_psm.py psm.tsv "msfragger" psm_with_ms.tsv
```
#### 2.2 Remove duplicated mass shifts with greedy algorithm
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
#### 2.3 Get peptides with modifications and positions (MaxQuant).
Input:
* peptides.txt
* modificationSpecificPeptides.txt

Output:
* mod_peptides.txt

Run the command:
```sh
python3 merge_peptide_df.py txt/peptides.txt txt/modificationSpecificPeptides.txt mod_peptides.tsv
```
### 3. Top-down MS database search
#### 3.1 Database search using TopPIC
Requirements
* .mzML file of top-down spectra
* protein database in fasta format
* installed TopPIC Suite
* PTM list pre-defined (Optional)

Perform the open search following the instructions [here](https://www.toppic.org/software/toppic/tutorial.html).

#### 3.2 Database search using MSPathFinder
Download and install MSPathFinder and perform database search using it following the instructions [here](https://github.com/PNNL-Comp-Mass-Spec/Informed-Proteomics).

### 4. Preprocess mass shifts from proteoform identifications

#### 4.1 Get proteoforms with mass shifts
Input parameter:
* .tsv file containing identifications with mass shift information, e.g. toppic_proteoform_single.tsv in the output files of TopPIC
* output file name

Output:
* .tsv file containing proteoforms with mass shifts

Run the command:
```sh
python3 get_proteoforms_with_ms.py toppic_proteoform_single.tsv proteoform_with_ms.tsv
```
#### (Optional) 4.2 Remove proteoforms from histone proteins
Input parameter:
* .tsv file containing proteoforms with mass shifts
* output file name

Output:
* .tsv file containing proteoforms without those from histone proteins

Run the command:
```sh
python3 remove_histone.py proteoform_with_ms.tsv proteoform_with_ms_no_histone.tsv
```
#### 4.3 Remove duplicated mass shifts and split mass shifts into N-terminal acetylation and other PTMs
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
#### 4.4 Extract mass shift information
Input parameter:
* .tsv file containing proteoforms with mass shifts
* output file name

Output:
* .tsv file containing mass shift information, including mass, possible range of localization, and annotation

Run the command:
```sh
python3 extract_ptm_to_df.py mass_shift.tsv mass_shift_with_info.tsv 
```
#### 4.5 Assign PTM type to mass shifts from high-frequency ones with a user-defined error tolerance
Input parameter:
* .tsv file containing mass shift information
* output file name for assigned mass shifts
* output file name for unassigned mass shifts
* error tolerance 

Output:
* .tsv file containing assigned mass shifts
* .tsv file containing unassigned mass shifts

Run the command:
```sh
python3 find_confident_ptm.py mass_shift_with_info.tsv mass_shift_assigned_with_high_frequency.tsv mass_shift_not_identified.tsv 0.1 
```
### 5. Verifying mass shifts from proteoforms with mass shifts from peptides
#### 5.1 Verifying mass shifts using MS-Fragger and MetaMorpheus results
Input parameter:
* input1: .tsv file containing top-down mass shift information
* input2: .tsv file containing bottom-up mass shift information
* t: error tolerance
* output: output file name for the matched entries
* n-term-mode: 0 or 1, 0 represents matching for mass shifts except N-terminal acetylation and 1 represents matching for N-terminal acetylation 

Output:
* .tsv file containing matched entries

Run the command:
```sh
python3 find_open_search_evidence.py input1=mass_shift_top_down.tsv input2=mass_shift_bottom_up.tsv t=0.1 output=mass_shift_td_matched_with_bu.tsv n-term-mode=0
```
#### 5.2 Verifying mass shifts using MaxQuant results
Input:
* .tsv file containing top-down mass shift information
* mod_peptides.tsv from step 2.3
* evidence.txt in outputs of MaxQuant
* type of ptm, e.g. Phospho
* monoisotopic mass of ptm, e.g. 79.9663
* error tolerance
* output file name of PTM with verified evidences

Output:
* PTM_with_evidence.tsv

Run the command:
```sh
python3 find_ptm_evidence.py mass_shift_top-down.tsv mod_peptides.tsv txt/evidence.txt Phospho 79.9663 0.1 phospho_evidence.tsv
```

### 6. PTM annotations extraction
#### 6.1 UniProt
Download Entries "Genes" and "Modified residues" for the species Homo Sapiens to a list. (You can download the extracted list here.)
```sh
python3 download_uniprot_ptm_annotation.py output.tsv
```

#### 6.2 Preprocess PTM annotations from dbPTM
Download the annotation of experimental PTM sites in dbPTM for specified PTM types [here](https://awi.cuhk.edu.cn/dbPTM/download.php).

Extract PTM sites, e.g. phosphorylation, about the species Homo Sapiens using the following command:
```sh
grep "_HUMAN" Phosphorylation > human_anno/phosphorylation_human.tsv
```

### 7.  Verifying mass shifts from proteoforms with PTM annotations
#### 7.1 Verifying using UniProt annotations
Input parameter:
* input1: .tsv file containing top-down mass shift information
* input2: .tsv file containing PTM annotations
* t: error tolerance
* output: output file name for the matched entries 

Output:
* .tsv file containing matched entries

Run the command:
```sh
python3 find_uniprot_evidence.py input1=mass_shift_top-down.tsv input2=uniprot_preferred_gene_modified_residue.tsv t=0.1 output=mass_shift_td_matched_with_anno.tsv
```
#### 7.2 Verifying using dbPTM annotations 
Input parameter:
* input1: .tsv file containing top-down mass shift information
* input2: .tsv file containing PTM annotations
* output: output file name for the matched entries
* n-term-mode: 0 or 1, 0 represents matching for mass shifts except N-terminal acetylation and 1 represents matching for N-terminal acetylation 

Output:
* .tsv file containing matched entries

Run the command:
```sh
python3 find_dbPTM_evidence.py input1=mass_shift_top-down.tsv input2=phosphorylation_human.tsv output=mass_shift_td_matched_with_anno.tsv n-term-mode=0
```

### 6. Extract PTM information from UNIMOD
Extract PTM information as follows:
* Name: name of the modification (Unimod PSI-MS name)
* Mass: monoisotopic mass of modification
* Residues: amino acids that can be modified
* Position: positions in the protein where the modification can be attached
* UnimodID: unmimod id of the modification

Input:
* The entire contents of Unimod in OBO format, you can download it [here](https://www.unimod.org/downloads.html).
* output file name for the extracted PTMs
Output:
* Extracted PTM list with information in txt format.

Run the command:
```sh
python3 extract_ptm_from_unimod.py unimod.obo unimod_ptm_list.txt
```

