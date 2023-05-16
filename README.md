# ProTN
ProTN is a novel R Shiny web app of an integrative pipeline that can perform all the steps of proteomics and phospho-proteomic analysis downstream of peptide quantification. ProTN works with any MS proteomic experiments managed with Proteome Discoverer and MaxQuant, two of the most used bioinformatic platforms to obtain peptide abundances from raw MS spectra. The implementation of ProTN focuses on being user-friendly, quick in doing the analysis, and comprehensive in the generation of images and tables that help the biological interpretation of results.

All the required information can be found on the info page of the web app.

## How execute the app
  1. Open the file **app.R** in RStudio.
  2. Click **Run App**
The app should automatically started installing all the required package.

## Workflow ProTN
ProTN is an integrative pipeline that analyze DDA proteomics data obtained from MS. It perform a complete analysis of the raw files from Proteome Discoverer (PD) or MaxQuant (MQ), with their biological interpretation with enrichement and network analysis. ProTN executes a dual level analysis, at protein and peptide level.

![image](https://github.com/TebaldiLab/ProTN/assets/39188419/5ce589c6-64f4-4af4-b708-538b05e983c5)

## Workflow PhosProTN
PhosProTN is an integrative pipeline for phosphoproteomic analysis of DDA experimental data obtained from MS. It perform a complete analysis of the raw files from Proteome Discoverer (PD) or MaxQuant (MQ), with their biological interpretation, enrichement and network analysis. 
PhosProTN analyse the phosphoproteomic data at peptide level, with background the proteomic analysis.

![image](https://github.com/TebaldiLab/ProTN/assets/39188419/147a2fc3-233a-4237-b37c-261e3e3d4bb3)

## Contacts
gabriele.tome@unitn.it

toma.tebaldi@unitn.it
