# Cardo-et-al-2021
Scripts for analysis of raw sequencing data and differential gene expression for the SETBP1 paper from Lucia F. Cardo. 
The scripts start with raw data in .fastq files and uses STAR to align the sequences to a reference genome (obtained from Ensembl).


To use:

  - Edit param.txt file to specify the genome assembly, the read length, number of samples and versions of each module to load into Raven.

Each script is run sequentially for each group (e.g. all 01_* scripts can be run at the same time). Once they are done, the next set of scripts can be run.
Once the gene counts have been produced (all.featurecount.txt), these can be downloaded and analysed locally using the R scripts in DEG_analysis.
