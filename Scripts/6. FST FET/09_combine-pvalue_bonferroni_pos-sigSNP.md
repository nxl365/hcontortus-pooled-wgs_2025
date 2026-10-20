#  Replicate P-values were combined using Fisher's method; Bonferroni correction; positions of significantly differentiated SNPs
Only for the 3 most divergent farm pairs:
s4 VS s15,  s3 VS s4,  G (has 2 replicates s7, s9 ) VS s4



## 1. replicate samples' P-values from FET were combined using Fisher's method

### input 
```python
df_fet_pval_1= df_fet_pval[["CHR", "POS", "Aver_min_coverage",sample_pair_1, sample_pair_2]]
df_fet_pval_1
```
```bash
# this is how the df_fet_pval_1 looks like:
# s7 and s9 are 2 replicates from farm G
CHR	POS	Aver_min_coverage	s4_VS_s7	s4_VS_s9
0	chr1	3805	15.0	0.170683	1.000000
1	chr1	4124	20.0	1.000000	0.507250
2	chr1	4625	16.0	0.434478	0.019523
3	chr1	4635	17.0	0.031004	0.031004
4	chr1	4636	16.0	0.068424	0.068424
...	...	...	...	...	...
16260893 rows × 5 columns
```

### python code
```python
######### faster for a big dataset, use numpy and scipy to calculate manually

## 2 p values
p1=df_fet_pval_1[sample_pair_1]
p2=df_fet_pval_1[sample_pair_2]    

## combine p value using Fisher's method
stat = -2 * (np.log(p1) + np.log(p2))  ### Fisher’s test statistic (X²) from the two p-values
p_comb = chi2.sf(stat, df=4) # Convert X² into the final combined p-value (right-tail probability of chi-square with df=4), k = replicate, k=2， so df=2k=4

df_fet_pval_1[sample_pair_name] = p_comb
df_fet_pval_1

```
### output

```bash
## this is how the df_fet_pval_1 looks like:
## G_VS_s4 is the combined p-values
CHR	POS	Aver_min_coverage	s4_VS_s7	s4_VS_s9	G_VS_s4
0	chr1	3805	15.0	0.170683	1.000000	0.472441
1	chr1	4124	20.0	1.000000	0.507250	0.851547
2	chr1	4625	16.0	0.434478	0.019523	0.048941
3	chr1	4635	17.0	0.031004	0.031004	0.007639
4	chr1	4636	16.0	0.068424	0.068424	0.029796
```




## 2. Bonferroni correction
α/m (0.05/16260893) = 0.05 is divided by the total number of genome-wide SNPs tested

```python
## total number of genome-wide SNPs tested
genome_wide_snps_nr=len(df_fet_pval)
genome_wide_snps_nr  # 16260893


## Bonferroni corrected p value
bon_p_cutoff=0.05/genome_wide_snps_nr 
print(f"Bonferroni corrected p-value cutoff: {bon_p_cutoff}")
# Bonferroni corrected p-value cutoff: 3.074861878741838e-09


## then keep only significantly differentiated SNPs < bon_p_cutoff instead of 0.05

```


## 3. create the list of positions of significant differentiated SNPs

this is for the next step: using bcftools to extract these significantly differentiated SNPs with annotation in annotated VCF files

```python
sig_snp_list = sig_inter_snp[['CHR','POS']]

###  sig_snp_list: should have the same CHR values in VCF files
	CHR	                                    POS
0	hcontortus_chr1_Celeg_TT_arrow_pilon	12834921
1	hcontortus_chr2_Celeg_TT_arrow_pilon	16950617
2	hcontortus_chr2_Celeg_TT_arrow_pilon	16950681
3	hcontortus_chr2_Celeg_TT_arrow_pilon	16950694
4	hcontortus_chr2_Celeg_TT_arrow_pilon	16950710
......
```

### output 
the positions of significantly differentiated SNPs for the 3 most divergent farm pairs are uploaded to `Analysis_results/6. FST/06_pos-sigSNP`

```bash
G_s4_sig_snp_pos_list-windowFst_10000_Cut1-everysnpFst_cut1.txt
s3_s4_sig_snp_pos_list-windowFst_10000_Cut1-everysnpFst_cut1.txt
s4_s15_sig_snp_pos_list-windowFst_10000_Cut1-everysnpFst_cut1.txt

```



