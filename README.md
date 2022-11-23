This software is an [R Markdown](http://rmarkdown.rstudio.com) Notebook of an integrative pipeline that analyze proteomics data obtained from MS.
This pipeline take in input the raw files from [Proteome Discoverer](https://www.thermofisher.com/it/en/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/proteome-discoverer-software.html) (PD) or [MaxQuant](https://maxquant.org/) (MQ). In the first phases it imputes and normalizes the data obtaining the matrix that is used in the enrichment analysis of the proteins. The enrichment perform a complete analysis in [EnrichR](https://maayanlab.cloud/Enrichr/) obtaining the results from all the DBs in it, and a PPI network analysis via [STRING](https://string-db.org/), During the enrichment steps, the only protein analyzed are the up- or down-regulated, discriminated by the following parameters. These parameters are defined in the *options.xlsx* file that is required for the execution since it contains the parameters required for the execution.

## Default options for the up-/down-regulated discovery:
- Signal log2 expression threshold > -∞ (No Limit, represent by value "inf" in the cell)
- log2 Fold Change threshold = 0.75 (Up-regulated > 0.75, Down-regulated < -0.75)
- P value threshold < 0.05

## Before starting
Get the software by clone this repository in your computer. It is possible to download it from git in the terminal or directly from this website downloading the zip of the repository.  
The pipeline require 4 files:
- **options.xlsx**: contain all parameters needed. Need to be in the same folder of this notebook. Later is reported the description of the fields.
- **Input[...].xslx**: describe the samples of the experiment. Most important it identify the Condition of each sample as we see later.
- **PEP[...]**: raw file of peptides obtained from PD or MQ.
- **PROT[...]**: raw file of protein groups obtained from PD or MQ.

Before run the pipeline, we need to modify the *Input[...].xslx* file and the *options.xlsx* file.

### How create/modify the *Input[...].xslx*
- ```Condition``` column (REQUIRED): define the condition of each sample that divide the samples in groups. The conditions need to be the same of the the Contrast Design. In the example the two condition are DMSO and TRT.
- ```Color``` columns (OPTIONAL): define a color for the samples in the graphs. If not present use a default palette.
- ```MS_batch``` columns: define the groups of batch in the samples. REQUIRED FOR BATCH EFFECT CORRECTION.
- ```Sample``` column: define the names for the samples.
  - In case of ```PD``` files use the *Input* file obtained from PD, an example of PD Input file is present in *file_input_raw/Input_filter.xlsx*, the **Sample** column in optional, if is not present the software extract the names for the **"File Name"** column.
  - In case of ```MQ``` analysis the  this column is REQUIRED

### How modify the *options.xlsx*
- ```Title of Analysis```: title of the experiment. It will be the title of the web page report.
- ```Description```: description of the current experiment. It is the first paragraph of the report.
- ```Software Analyzer```: determine with software was use to identify peptides and proteins. **PD** for Protein Discoverer, **MQ** for MaxQuant.
- ```Path input```: folder that contain the 3 raw files from PD. The 3 files need to contain in the file name: *"Input"* for the file with the information about the samples, *"PEP"* in the filename of the peptide file data, *"PROT"* in the filename of the protein file. (defaults: "file_input_raw/")
- ```Signal log2 expr thr```: possible to change the signal log2 expression threshold. Default: no limit ("inf")
- ```Log2 FC thr```: possible to change the Fold Change threshold
- ```P.Value thr```: possible to change the p.value threshold (suggest to use the default 0.05)
- ```Batch Correction```: boolean value for the execution of the batch effect correction performed by [proBatch](https://www.bioconductor.org/packages/release/bioc/html/proBatch.html). (TRUE or FALSE)
- ```Contrast Design```: write the formulas of the contrast comparison you want to analyse with [Limma](https://bioinf.wehi.edu.au/limma/). **AT LEAST 1 FORMULA IS REQUIRED**, so the column *"Formule"* cannot be empty. Instead the column *"Name"* is optional, it assign a personalized name to the rule in the same row. Example of contrast design: *TRT-DMSO*, where TRT are the treatment samples group, and DMSO the wild type. TRT and DMSO are the same of the *condition* column required in the Input files as previously described.
- ```Control Boxplot proteins```: list of proteins used as control of the intensities. For each protein a boxplot is generated comparing the mean of the intensities group by condition.
- ```Execute Enrichment```: boolean value for the execution of the enrichment step. (TRUE or FALSE)
- ```Execute PPI network STRINGdb```: boolean value for the execution of the network analysis. (TRUE or FALSE)
- ```P.Value thr for enrichment```: possible to change the p.value threshold of the enriched terms (suggest to use the default 0.05)
- ```Overlap size thr for enrichment```: possible to change the Overlap size threshold. The overlap size is the number of DEPs discovered in the enriched terms. (suggest to use the default 5)
- ```Enrichment```: determine of filter the plot of the enrichment results. The enrichment can be filter by terms or by DBs. In the column *"Terms to search"* you can write the word that you want to search in the results of EnrichR (EX: MYC, C-MYC, Senescence,...). Instead, in the column *"DB to analyse"* you can write the DBs that you want to see in your plots (EX: ChEA_2016, KEGG_2021_Human, BioPlanet_2019, GO_Biological_Process_2021,...), using the same name that you can find in EnrichR. If the two columns are empty, no plots of enrichment are returned.

## Example case study - Proteomics of MCF7 cells
Proteomics of MCF7 cells (LC-MS/MS Orbitrap Fusion Tribrid). PRIDE: [PXD009417](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD009417). Its execution is performed using the provided *options.xlsx* file.
Clamer M, Tebaldi T, Lauria F, et al. Active Ribosome Profiling with RiboLace. Cell Rep. 2018;25(4):1097-1108.e5. 


## Details on the output files

- **report.html**: complete report of the analysis with all pics and results of the enrichment.
- **database_env_R.RData**: RData containing the essential variables for further analysis.
- **normalised_intensity_table.xlsx**: excel file with the matrix used in the analysis. In sheet *Normalized proteins matrix intensity* is written the normalized and imputed version of the original intensities of each protein in each sample. In sheet *Normalized peptides matrix intensity* is written the normalized and imputed version of the original intensities of each peptide in each sample. In sheet *Normalized proteins intensity group by* is written the average and the standard deviation of each condition described in the Input file for each protein. In sheet *Normalized peptides intensity group by* is written the average and the standard deviation of each condition described in the Input file for each peptide.
- **DEPs_table.xlsx**: summary of the result of Limma. In sheet *DEPs_proteins*, Column *id* is the name of the protein. In sheet *DEPs_peptides*, column *Accession* is the UniprotID of protein, column *Description* is the description of protein, *GeneName* is Gene Symbol, *Annotated Sequence* is the sequence of peptides, *Modifications* is the modification that affect the peptides, *Position in Master Proteins* is the UniprotID with the position on the sequence. For each contrast rule are present 5 columns:
  - *class*: is "+" if the protein is up-regulated following the parameters set in **options.xlsx**, "-" if the protein is down-regulated, "=" otherwise
  - *log2_FC*: is the log2 Fold Change of the contrast
  - *p_val*: is the p. value
  - *p_adj*: is the adjusted p. value
  - *log2_expr*: is the signal log2 expression

- **general_mds_proteins.pdf**: MDS of the samples by proteins.
- **general_mds_peptides.pdf**: MDS of the samples by peptides.
- **corrected_mds_proteins.pdf**: MDS of the samples by proteins after the bacth effect correction.
- **corrected_mds_peptides.pdf**: MDS of the samples by peptides after the bacth effect correction.
- **boxplot_control_intensity.pdf**: boxplot for the control of the intensities of selected proteins.
- **DEPs_proteins_count_barplot.pdf**: number of DEPs find for each contrast design based on proteins. Light grey are the down-regulated, dark grey are the up-regulated genes.
- **DEPs_peptides_count_barplot.pdf**: number of DEPs find for each contrast design based on peptides. Light grey are the down-regulated, dark grey are the up-regulated genes.
- **DEPs_proteins_mds.pdf**: MDS of the DEPs discovered based on proteins.
- **DEPs_peptides_mds.pdf**: MDS of the DEPs discovered based on peptides.
- **DEPs_enrichment_data.RData**: result of the enrichment made with EnrichR.
- **DEPs_enrichment_table.xlsx**: excel with only the enriched terms with an high significance (P.Value < 0.05, Overlap Size >= 4).
- **enrich** directory: plots of the enrichment data under the parameters described in the **option.xlsx** file
  - **DBs_enr.pdf**: dot matrix of the enrichment result divided in up- and down-regulated genes. Filtered by the datasets choose in the **options** file.
  - **terms_enr.pdf**: dot matrix of the enrichment result divided in up- and down-regulated genes. Filtered by the terms written in the **options** file.
  - **DBs_enr_all.pdf**: dot matrix of the enrichment result combining the up- and down- regulated genes in a single list. Filtered by the datasets choose in the **options** file.
  - **terms_enr_all.pdf**: dot matrix of the enrichment result combining the up- and down- regulated genes in a single list. Filtered by the terms written in the **options** file.

- **network** directory: plots of the up- and down- regulated genes based on STRINGdb. For each contrast design are create 2 files:
<!-- - **Prefilter**: histogram representing the number of genes in each community discovered. -->
  - **Postfilter**: histogram representing the number of genes in each community discovered after a filter on the genes and on the links.
  - **Graph**: visual representation of the network discovered in two layouts: *fr* and *kk*.

## Contacts
gabriele.tome@studenti.unitn.it

toma.tebaldi@unitn.it
