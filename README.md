# ASE
Allele Specific Expression Pipeline
This ASE pipeline is developed by using publicly available RNA-seq data and open-source software. Using this pipeline, we are able to profile pervasive allelic imbalance across 42 tissues and 34 breeds from the Farm-GTEx-pig consortium at both SNP and gene levels without the need for parental genotypes or whole genome sequence data.
This repository contains analysis pipelines used by the PigGTEx Consortium, including: 
* The pipline for RNA-seq alignment, quantification, SNP calling, genotype imputation, and quality control is similar as the FarmGTEc-PigGTEx. (This pipeline is available [here](https://github.com/FarmGTEx/PigGTEx-Pipeline-v0/tree/53064b6528079b3fbb71b58cbb4415b82de69fd4/02_RNA-Seq) );
* Allele-specific expression (ASE) analysis pipeline cantains mapping bias remove, ASE sites/genes detectation using phaseR. The phased vcf (individual-level imputed genotype was used in this study which from Farm-GTEx), fastq data.
