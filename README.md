# hcontortus-pooled-wgs_2025

This repository contains scripts for the paper:
**Whole Genome Pool Sequencing Reveals Genomic Variation in Field Haemonchus contortus Isolates in Sweden**

This work establishes the first genome-wide baseline for field-derived H. contortus populations across multiple Swedish sheep farms and shows that pooled-sample whole-genome sequencing (Pool-seq) is a useful method for describing broad population-genomic patterns in this species. Our findings indicate relatively low genomic variation between Swedish H. contortus populations and high genomic variation within them. This study provides a valuable foundation for future research on parasite adaptation, resistance, surveillance, and control.

## Suggested citation
If you use this repository, please cite this paper, DOI: [10.1186/s12864-026-13449-5](https://doi.org/10.1186/s12864-026-13449-5)

## Samples

The analysis includes pooled *H. contortus* samples from Swedish sheep farms.
1. **14 farms**
2. **17 pooled samples** : each sample comprising a pool of at least 2,000 larvae; 2 biological replicate samples from 3 farms (A, F, and G) each
    * 2 samples were excluded from downstream analyses due to poor mapping quality
3. **whole genome sequencing of pooled samples** - Illumina NovaSeq X Plus 
4. ***Haemonchus contortus*** reference genome (~283 Mbp) and **GFF3 annotation** file (WormBase release 2021-04) used are available in the [WormBase repository](https://parasite.wormbase.org/Haemonchus_contortus_prjeb506/Info/Index/)


### Sample Metadata

The table below shows the correspondence between farm sample IDs (A–N) in the paper and the old laboratory sample IDs (s1–s17).

| Sample ID | Laboratory Sample ID |
|-----------|----------------------|
| A1        | s1                   |
| F1        | s2                   |
| C         | s3                   |
| D         | s4                   |
| E         | s5                   |
| A2        | s6                   |
| G1        | s7                   |
| B         | s8                   |
| G2        | s9                   |
| F2        | s10                  |
| H         | s11                  |
| N         | s12                  |
| I         | s13                  |
| L         | s14                  |
| K         | s15                  |
| M         | s16                  |
| J         | s17                  |

* Note: Samples A1/A2, F1/F2, and G1/G2 are independent biological replicates from the same farms. Samples M (s16) and J (s17) were excluded from downstream analyses



## Availability of data generated 

* The SNP datasets generated and analysed during the current study are available in the European Nucleotide Archive under the study accession number PRJEB112099.

* Some other data generated and analysed during the current study are available in this Github repository: Data, Analysis_results



## Main analyses and scripts

1. Quality control (QC) and trimming  ✅
2. Mapping and post-mapping processing  ✅
3. Variant calling  ✅
4. Variant filtering  ✅
5. Variant annotation  ✅
6. Between-population differentiation: fixation index (FST),Fisher's exact test (FET)  ✅
7. Prepare for gene ontology (GO) enrichment analysis  ✅
8. Allele frequency and principal component analysis (PCA): based on SNPs ✅
    * [Interactive 3D PCA plot](https://github.com/nxl365/hcontortus-pooled-wgs_2025/blob/main/Analysis_results/8.%20PCA/PCA_3D_genome_interactive.html)
9. Within-population diversity: nucleotide diversity (Pi) ✅


## License

This work is licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/), unless otherwise stated.

