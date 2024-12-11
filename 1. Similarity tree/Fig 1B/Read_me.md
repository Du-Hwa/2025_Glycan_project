# README
These folders contains the source code and associated data for creating the similarity tree depicted in Figure 1B.

## Folder Structure and Required Data

### Data Requirements
Each folder contains data specific to the corresponding step. However, some data files required for the analysis are not distributed in this repository due to licensing restrictions. Please download these files from the following sources:

- **[LysM-RLK.fasta.gz](https://doi.org/10.5281/zenodo.7017981)**
- **[LysM-RLP.fasta.gz](https://doi.org/10.5281/zenodo.7017981)**
- **[all.LRR-RLK.fasta.gz](https://doi.org/10.5281/zenodo.7017981)**
- **[RLP.fasta.gz](https://doi.org/10.5281/zenodo.7017981)**

**Important:**  
- Rename `all.LRR-RLK.fasta.gz` to `LRR-RLK.fasta.gz` and `RLP.fasta.gz` to `LRR-RLP.fasta.gz`.  
- Place these renamed files in the folder `01.Seq_from_2022_Nat_plants_Ngou_et_al`.

Additionally, download the following files from TAIR:  
- **[TAIR10_pep_20101214](https://www.arabidopsis.org/download/list?dir=Proteins%2FTAIR10_protein_lists)**  
  Place this file in the folder `02.missing_sequence_from_TAIR`.  

- **[TAIR10_pep_20110103_representative_gene_model](https://www.arabidopsis.org/download/list?dir=Proteins%2FTAIR10_protein_lists)**  
  Place this file in the folder `06.Tree iToL/datafiles`.

### Step-by-Step Workflow
1. At each step, the final output file must be manually moved to the next step's folder.  
   For convenience, some final files have already been pre-copied into subsequent folders.  

2. **Step 5: IQ-TREE**  
   - The phylogenetic tree was constructed using the web-based tool [IQ-TREE](http://www.iqtree.org/).  

3. **Final Tree and Decoration**  
   - The final tree was visualized and decorated using [iTOL](https://itol.embl.de/).  
   - Files for the decoration were generated using the script `06.Tree iToL/iTOL_annotationfiles/06-2.iTOL_annotation.ipynb`.  
   - These files are stored in the folder `06.Tree iToL/iTOL_annotationfiles` and were applied to complete the tree decoration.
