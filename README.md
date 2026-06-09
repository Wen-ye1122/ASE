# ASE
Allele Specific Expression Pipeline. 
This ASE pipeline is developed by using publicly available RNA-seq data and open-source software. Using this pipeline, we are able to profile pervasive allelic imbalance across 42 tissues and 34 breeds from the Farm-GTEx-pig consortium at both SNP and gene levels without the need for parental genotypes or whole genome sequence data.
This repository contains analysis pipelines used by the PigGTEx Consortium, including: 
* The pipline for RNA-seq and DNA-seq alignment, quantification, SNP calling, and quality control is similar as the FarmGTEc-PigGTEx. (This pipeline is available [here](https://github.com/FarmGTEx/PigGTEx-Pipeline-v0) );
* Allele-specific expression (ASE) analysis pipeline cantains mapping bias remove, ASE sites/genes detectation using phaseR. Binomal test for ASE sites/gene.
* The following files are required: the phased vcf (individual-level imputed genotype, in this study the files are from Farm-GTEx), fastq data.
