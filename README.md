# hcontortus-pooled-wgs_2025

This repository contains scripts for the manuscript:
**Whole Genome Pool Sequencing Reveals Genomic Variation in Field Haemonchus contortus Isolates in Sweden**

## Abstract
[WILL PASTE PAPER ABSTRACT HERE]


## Samples & Data

The analysis includes pooled *H. contortus* samples from Swedish sheep farms. After quality control, the final dataset contained:

* **12 farms**
* **15 pooled samples** : each sample comprising a pool of at least 2,000 larvae; 2 biological replicate samples from 3 farms (A, F, and G) each

* **whole genome sequencing of 15 pooled samples** - Illumina NovaSeq X Plus 
* ***Haemonchus contortus*** reference genome (~283 Mbp) and **GFF3 annotation** file (WormBase release 2021-04) used are available in the [WormBase repository](https://parasite.wormbase.org/Haemonchus_contortus_prjeb506/Info/Index/)



## Availability of data generated 

The datasets generated and analysed during the current study will be deposited in:
**European Nucleotide Archive**

Accession number:
**PRJEB112099**  

Data list:
[WILL UPDATE HERE]



## Main analyses and scripts

Pipeline figure:[UPLOAD the final version later]


1. Read quality control (QC) and trimming ✅
2. Mapping, sorting, duplicate marking
3. Variant calling and filtering
4. Variant annotation: SNPs and indels
5. Between-population differentiation: fixation index (FST),Fisher's exact test (FET)
6. Gene ontology (GO) enrichment analysis
7. Allele frequency and principal component analysis (PCA): based on SNPs
8. Within-population diversity: nucleotide diversity (Pi), Watterson’s theta, Tajima’s D
9. Figure generation

## Suggested citation

If you use this repository, please cite this paper **Whole Genome Pool Sequencing Reveals Genomic Variation in Field Haemonchus contortus Isolates in Sweden** [LINK]


## License

This work is licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/), unless otherwise stated.