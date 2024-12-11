# README

This folder contains codes for RNA-Seq data analysis and visualization, designed to generate the figures in the manuscript.

## Folder Structure and Workflow

### 1. Preparation
#### Quality Check and Quantification
- The quality of raw reads was assessed using **FastQC (ver 0.11.9)**.  
- The reads were quantified against the **TAIR10 reference transcriptome (TAIR10.transcripts.fa)** using **Kallisto (ver 0.46.2)**.  

#### Input Data Placement
- Place the processed data in the `Sequencing_data/Kallisto` folder.  
- Ensure the folder structure follows the specifications in `Sequencing_data/Kallisto/studydesign.txt`.  
  For example:  
  - `Sequencing_data/Kallisto/RHAM_0min_1` should contain the following files:  
    - `abundance.h5`  
    - `run_info.json`  
    - `abundance.tsv`  
- A total of **9 subfolders** should be present under `Sequencing_data/Kallisto`.  

#### Differential Expression Analysis (DESeq2)
1. **RG-I Data**  
   - Use `01-1.DESeq_RG-I.ipynb` to perform DESeq2 analysis on the above data.  

2. **OGs Data**  
   - For OGs, RNA-Seq results from *Bjornson et al., Nat. Plants, 2021* were used.  
     Reference: [http://10.1038/s41477-021-00874-5](http://10.1038/s41477-021-00874-5)  
   - Data is accessed using **getDEE2**, so no additional downloads are necessary.  
   - Ensure `Sequencing_data/NP/NP_ex_design.csv` is available in the specified path.  
   - Use `01-2.DESeq_OGs.ipynb` to perform DESeq2 analysis for OGs.

---

### 2. Step-by-Step Execution
The analysis must be performed sequentially as each step uses output files from the previous step.

1. **02.Volcano_plot.ipynb**  
   - Generates the **Supplementary Figure 5C** (Volcano plot).  

2. **03.PCA_plot.ipynb**  
   - Generates the **Supplementary Figure 5B** (PCA plot).  

3. **04.Heatmap.ipynb**  
   - Generates **Figure 3E** (Heatmap).  

4. **05.Venndiagram.ipynb**  
   - Generates **Figure 3F** (Venn diagrams).  

5. **06.New_GO_bubble_chart.ipynb**  
   - Generates **Figure 3G** and **Supplementary Figure 5D** (Bubble plot).  

---

### Notes
- Follow the execution order carefully to ensure compatibility between steps.  
- Each notebook contains detailed documentation for its specific step.  
