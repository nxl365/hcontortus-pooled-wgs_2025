# Summary of Significant SNPs with annotations
will get cleaned gene list for GO enrichment analysis


## 1. Extract useful info of annotated VCFs into a TSV table

```bash
cd '/cfs/klemming/projects/snic/naiss2025-23-66/scripts/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-07_sig_SNP'

02_extract_sig_snp_anno_tsv.sh

```

### script content

```bash
#!/bin/bash -l

module load bcftools

### arguments: need when you run the script
SAMPLE_ID_1=$1
SAMPLE_ID_2=$2

windowFst_cut=$3
everysnpFst_cut=$4


## configs: setting in the script

dir_base="/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-07_sig_SNP/${SAMPLE_ID_1}_VS_${SAMPLE_ID_2}/windowFst_10000_Cut${windowFst_cut}-everysnpFst_cut${everysnpFst_cut}"
cd $dir_base


input="${SAMPLE_ID_1}_${SAMPLE_ID_2}_sig_snp_anno_cut${windowFst_cut}_${everysnpFst_cut}.vcf"

output="${SAMPLE_ID_1}_${SAMPLE_ID_2}_sig_snp_anno_cut${windowFst_cut}_${everysnpFst_cut}_extracted.tsv"

# main code
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/ANN\n' $input \
| awk -F'\t' 'BEGIN{OFS="\t"}{       ## input and output delimiter: \t
    split($7,a,",");            ## splite the ANN col by comma into array a[], can be different transcript
    split(a[1],b,"|");          ## splite the first annotation/transcript by | to array b[], which will be allele|annotation....
    print $1,$2,$3,$4,$5,$6,b[4],b[2],b[3],b[10],b[11];
}' \
> $output

## the columns you extract: 
# $1..$6：CHROM、POS、REF、ALT、QUAL、FILTER
# $7   b[4]：Gene_Name
# $8   b[2]：Annotation （Consequence）--- missense, nonsense, splice..., might have multiple consequences and joined with `&`
# $9   b[3]：Annotation_Impact-- LOW, HIGT...; only one impact
# $10  b[10]：HGVSc
# $11  b[11]：HGVSp

```

### run the script

```bash
# e.g. for 3 pairs
bash 02_extract_sig_snp_anno_tsv.sh s4 s15 1 1

bash 02_extract_sig_snp_anno_tsv.sh s3 s4 1 1

bash 02_extract_sig_snp_anno_tsv.sh G s4 1 1

```


### input
```bash

s4_s15_sig_snp_anno_cut1_1.vcf

s3_s4_sig_snp_anno_cut1_1.vcf

G_s4_sig_snp_anno_cut1_1.vcf

```

### output
uploaded these outputs to `Analysis_results/7.prepare for GO enrichment`
```bash
s4_s15_sig_snp_anno_extracted.tsv

s3_s4_sig_snp_anno_cut1_1_extracted.tsv

G_s4_sig_snp_anno_cut1_1_extracted.tsv

```



## 2. Extract genes, summarize annotation of SNPs, get cleaned genes for GO enrichment

```bash
cd '/cfs/klemming/projects/snic/naiss2025-23-66/scripts/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-07_sig_SNP'

03_extract_genes_stats.sh
```

### script content
```bash

#!/bin/bash -l



## arguments
SAMPLE_ID_1=$1
SAMPLE_ID_2=$2

windowFst_cut=$3
everysnpFst_cut=$4



## configs
dir_base="/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/06-FST_removed_s16_s17/06-07_sig_SNP/${SAMPLE_ID_1}_VS_${SAMPLE_ID_2}/windowFst_10000_Cut${windowFst_cut}-everysnpFst_cut${everysnpFst_cut}"

cd $dir_base


file_base="${SAMPLE_ID_1}_${SAMPLE_ID_2}_sig_snp_anno_cut${windowFst_cut}_${everysnpFst_cut}"


input_vcf="${file_base}.vcf"

input_tsv="${file_base}_extracted.tsv"



########## 1. count the number of SNPs with annotations ##########
output_count_snp="${file_base}-count_snp.txt"

#### 1.1 total number of SNPs with annotations
echo "Total number of SNPs with annotations:" > $output_count_snp # create/overwrite file
grep -v "^#" $input_vcf | wc -l >> $output_count_snp # append more
echo >> $output_count_snp  # add an empty line



#### 1.2 annotated SNP counts per chromosome
echo "Number of SNPs per chromosome with annotations:" >> $output_count_snp # create/overwrite file
grep -v "^#" $input_vcf | cut -f1 | sort | uniq -c >> $output_count_snp # append more






########## 2. extract the list of genes involved and their SNP counts ##########

#### 2.1 get: unique list of all genes with SNPs
output_gene_list="${file_base}-uni_genes_all.txt"

cut -f7 $input_tsv | awk 'length($0)>0' | sort -u \
> $output_gene_list

## cut: cut the 7th col-- gene name;  
## awk: only keep the gene name if it is not empty; 
## sort -u: sort in alphabetical order and only keep unique gene name




#### 2.2 get: number of SNPs per gene (from most to least)
output_gene_list_order="${file_base}-uni_genes_counts.tsv"

cut -f7 $input_tsv | awk 'length($0)>0' | sort | uniq -c | sort -k1,1nr \
> $output_gene_list_order

## sort: sort in alphabetical order
## uniq -c: collapse duplicates and count how many times per gene appears
## sort -k1,1nr : sort by the 1st col (counts), numerically and in reverse order (most to least)





########## 3. Summarize the Annotation consequence, impact of SNPs ##########

#### 3.1 get the number of SNPs per annotation consequence
output_consequence="${file_base}-uni_consequences.tsv"

cut -f8 $input_tsv | tr '&' '\n' | sed '/^$/d' | sort | uniq -c | sort -k1,1nr \
> $output_consequence

## tr '&' '\n'： transfer & to \n, 
## since there are some snp with more than 1 consequence -- e.g. splice_region_variant&intron_variant; so put these multiple consequences to different lines

## sed '/^$/d' : delete if any line is empty
## uniq -c: collapse duplicates and count how many times per gene appears
## sort -k1,1nr : sort by the 1st col (counts), numerically and in reverse order (most to least)




#### 3.2 get the number of SNPs per annotation impact
## only one impact per SNP
output_impact="${file_base}-uni_impact.tsv"

cut -f9 $input_tsv | sed '/^$/d' | sort | uniq -c | sort -k1,1nr \
> $output_impact




#### 3.3 get the Impact → the corresponding 1st Consequence counts:
### the impact is only based on the 1st consequence when there are multiple consequences for one SNP


output_impact_conse="${file_base}-impact_consequence_counts.tsv"

awk -F'\t' '{
    split($8, cons, /[&]/);  # split consequences on &
    if (length(cons[1])>0) {
            print $9, cons[1];    ## only print 1st consequences and its corresponding impact
    }
}' $input_tsv \
| sort | uniq -c | sort -k2,2 -k1,1nr \
> $output_impact_conse




########## 4. cleaned genes lists for GO enrichment analysis ##########

## 4.1 get cleaned all SNP-mapped genes names (drop empty/intergenic combos, e.g. HCON1 - HCON2)

output_gene_list_clean="${file_base}-uni_genes_all_clean.txt"

grep -v '-' $output_gene_list | sort -u \
> $output_gene_list_clean


## 4.2 Extract the gene names with HIGH impact SNPs
output_impact_MH="${file_base}-uni_genes-impact_High.txt"

awk -F'\t' '$9 ~ /HIGH/' $input_tsv \
| cut -f7 | sort -u \
> $output_impact_MH


## 4.3 Extract the gene names with Moderate impact SNPs

output_impact_MH="${file_base}uni_genes-impact_Moderate.txt"

awk -F'\t' '$9 ~ /MODERATE/' $input_tsv \
| cut -f7 | sort -u \
> $output_impact_MH



### if extract gene names with Moderate or HIGH: use
# awk -F'\t' '$9 ~ /MODERATE|HIGH/'

```


### run the script
```bash
# e.g. for 3 pairs

bash 03_extract_genes_stats.sh s4 s15 1 1

bash 03_extract_genes_stats.sh s3 s4 1 1

bash 03_extract_genes_stats.sh G s4 1 1

```



### input
```bash
# vcf
s4_s15_sig_snp_anno_cut1_1.vcf

s3_s4_sig_snp_anno_cut1_1.vcf

G_s4_sig_snp_anno_cut1_1.vcf

# tsv
s4_s15_sig_snp_anno_extracted.tsv

s3_s4_sig_snp_anno_cut1_1_extracted.tsv

G_s4_sig_snp_anno_cut1_1_extracted.tsv

```

### output
e.g. s4 VS s15
```bash
s4_s15_sig_snp_anno_cut1_1-count_snp.txt
s4_s15_sig_snp_anno_cut1_1-impact_consequence_counts.tsv
s4_s15_sig_snp_anno_cut1_1-uni_consequences.tsv
s4_s15_sig_snp_anno_cut1_1-uni_genes-impact_High.txt
s4_s15_sig_snp_anno_cut1_1-uni_genes_Moderate_High.txt
s4_s15_sig_snp_anno_cut1_1-uni_genes_all.txt
s4_s15_sig_snp_anno_cut1_1-uni_genes_all_clean.txt  ## this is the cleaned gene list for GO enrichment analysis! uploaded to the Analysis_results
s4_s15_sig_snp_anno_cut1_1-uni_genes_counts.tsv
s4_s15_sig_snp_anno_cut1_1-uni_impact.tsv
s4_s15_sig_snp_anno_cut1_1uni_genes-impact_Moderate.txt

```

