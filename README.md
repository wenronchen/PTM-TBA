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
    python3 get_proteoform_with_ms.py proteoform.tsv proteoform_with_ms.tsv
    
    # (Optional) Remove proteoforms from histone proteins
    python3 remove_histone.py proteoform_with_ms.tsv proteoform_with_ms_no_histone.tsv 
    # (Optional) Remove duplicated mass shifts
    python3 merge_simialr_mass_shifts.py proteoform_with_ms_no_histone.tsv proteoform_with_ms_no_dup.tsv 

