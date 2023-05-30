# PTM-TBA

A toolkit for extracting PTM information and verifing PTM information with top-down MS, bottom-up MS and UniProt annotations.

## PTM-TBA tutorials
### 1. Quick Start
##### 1.1  Run the examples of SW480 data sets in paper. 

    sh  PTM-TBA_one_step.sh -t 0.1
    
##### 1.2  Run with your own data sets.

    sh PTM-TBA_one_step.sh -td_input /path/to/proteoform.tsv -bu_input /path/to/peptides.tsv -anno uniprot_preferred_gene_modified_residue.tsv
    
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
# Assign PTM type to mass shifts from high-frequency ones with error tolerance 0.1 Da
python3 find_confident_ptm.py mass_shift_with_info.tsv 0.1 mass_shift_assigned_with_high_frequency.tsv mass_shift_not_identified.tsv

```
##### 2.2 Preprocessing mass shifts idenfied from bottom-up MS data
```sh
# Extract mass shifts from peptide identifications
python3 extract_ptm_from_psm.py psm.tsv psm_with_ms.tsv
# Remove duplicated mass shifts with greedy algorithm
python3 remove_dup_mass_shift_bu.py psm_with_ms.tsv 0.1 psm_with_ms_no_dup.tsv

```
##### 2.3  Verifying mass shifts from top-down MS with mass shifts from bottom-up MS
```sh



