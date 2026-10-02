# look for 3 known SNPs for BZ resistance

- **F167Y**：`TTC → TAC`，Phe → Tyr
- **E198A**：`GAA → GCA`，Glu → Ala
- **F200Y**：`TTC → TAC`，Phe → Tyr
- all in HCON_00005260 gene region on chr1

## 1. extract lines including 3 SNPs for BZ resistance in anotated VCF files

```bash
extract_BZ_resistance_snps.sh
```
```bash
#!/bin/bash -l

# Criteria for identifying the target beta-tubulin resistance SNPs:
# 1. The variant must be located on chromosome 1.
# 2. The variant must fall within the HCON_00005260 gene region (7,027,492-7,031,447 bp).
# 3. The SnpEff annotation must assign the variant to gene HCON_00005260.
# 4. The annotated protein change must match the target mutation:
#    p.Phe167Tyr, p.Glu198Ala, or p.Phe200Tyr.
# 5. The REF allele must be a single nucleotide.
# 6. The ALT allele responsible for the annotated consequence must also be a single nucleotide.
# Only variants meeting all of these criteria are retained as target resistance SNPs.


output="BZ_resistance_SNPs_by_sample.csv"

echo 'Sample_No,Sample_Name,Mutation,Status,CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,Sample_Field' > "$output"

for dir in $(printf '%s\n' Sample_XJ-4126-* | sort -V)
do
    [ -d "$dir" ] || continue

    sample_name="$dir"
    sample_no=$(echo "$sample_name" | sed -E 's/Sample_XJ-4126-([0-9]+).*/\1/')

    vcf="${dir}/${sample_name}_filtered_snp.snpeff.ann.vcf"

    if [ ! -f "$vcf" ]; then
        echo "WARNING: missing $vcf" >&2
        continue
    fi

    grep -E 'p\.Phe167Tyr|p\.Glu198Ala|p\.Phe200Tyr' "$vcf" | \
    awk -F '\t' \
        -v sample_no="$sample_no" \
        -v sample_name="$sample_name" \
        '
        function csv(x, y) {
            y=x
            gsub(/"/, "\"\"", y)
            return "\"" y "\""
        }

        function has_effect(info, effect,
                            fields, nfields, i, ann,
                            anns, nanns, j, a) {

            ann=""

            nfields=split(info, fields, ";")

            for (i=1; i<=nfields; i++) {
                if (substr(fields[i],1,4)=="ANN=") {
                    ann=substr(fields[i],5)
                    break
                }
            }

            if (ann=="")
                return 0

            nanns=split(ann, anns, ",")

            for (j=1; j<=nanns; j++) {

                split(anns[j], a, "\\|")

                # ANN fields:
                # a[1]  = allele
                # a[4]  = gene
                # a[10] = HGVS.c
                # a[11] = HGVS.p

                if (a[4]=="HCON_00005260" &&
                    a[11]==effect &&
                    length($4)==1 &&
                    length(a[1])==1) {

                    return 1
                }
            }

            return 0
        }

        function print_found(mut) {

            print \
            csv(sample_no) "," \
            csv(sample_name) "," \
            csv(mut) "," \
            csv("FOUND") "," \
            csv($1) "," \
            csv($2) "," \
            csv($3) "," \
            csv($4) "," \
            csv($5) "," \
            csv($6) "," \
            csv($7) "," \
            csv($8) "," \
            csv($9) "," \
            csv($10)
        }

        BEGIN {
            f167=0
            e198=0
            f200=0
        }

        $1=="hcontortus_chr1_Celeg_TT_arrow_pilon" &&
        $2>=7027492 &&
        $2<=7031447 {

            if (has_effect($8,"p.Phe167Tyr")) {
                print_found("F167Y")
                f167=1
            }

            if (has_effect($8,"p.Glu198Ala")) {
                print_found("E198A")
                e198=1
            }

            if (has_effect($8,"p.Phe200Tyr")) {
                print_found("F200Y")
                f200=1
            }
        }

        END {

            if (!f167)
                print csv(sample_no) "," csv(sample_name) "," csv("F167Y") "," csv("NOT_FOUND") ",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\""

            if (!e198)
                print csv(sample_no) "," csv(sample_name) "," csv("E198A") "," csv("NOT_FOUND") ",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\""

            if (!f200)
                print csv(sample_no) "," csv(sample_name) "," csv("F200Y") "," csv("NOT_FOUND") ",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\""
        }
        ' >> "$output"

done

echo "Done."
echo "Output: $output"
```


### run

```bash
bash extract_BZ_resistance_snps.sh
```

### input

one sample example
```bash
/cfs/klemming/projects/snic/naiss2025-23-66/analysis/variant_explore/illumina_HC_wgs_XJ4126/025_snpEFF_vcf_annota/025-02_snpEFF_annotation

|-- Sample_XJ-4126-10 # in each sample's folder
    |-- Sample_XJ-4126-10_filtered_snp.snpeff.ann.vcf
```
### output 
the output was uploaded to the folder `Analysis_results`
```bash
BZ_resistance_SNPs_by_sample.csv
```


## 2. calculate AF of SNPs for BZ resistance
For detected SNPs, the resistance-associated allele frequency was calculated from FreeBayes read counts as AO / (RO + AO) × 100, where AO and RO represent the numbers of reads supporting the alternative and reference alleles, respectively 

```bash
01_calculate_BZ_resis_SNP_AF.py
```
```python

import csv

input_file = "../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/025_snpEFF_vcf_annota/025-07_check_knownSNPs_drug_resis/BZ_resistance_SNPs_by_sample.csv"

output_file = "../../../../../analysis/variant_explore/illumina_HC_wgs_XJ4126/025_snpEFF_vcf_annota/025-07_check_knownSNPs_drug_resis/BZ_resistance_SNP_frequency.csv"


with open(input_file, newline="") as infile, open(output_file, "w", newline="") as outfile:

    reader = csv.DictReader(infile)

    fields = [
        "Sample_No",
        "Sample_Name",
        "Mutation",
        "CHR",
        "POS",
        "Gene",
        "REF",
        "ALT",
        "DP",
        "AO",
        "RO",
        "Frequency_percent"
    ]

    writer = csv.DictWriter(outfile, fieldnames=fields)
    writer.writeheader()

    for row in reader:

        # Only keep target SNPs that were found
        if row["Status"] != "FOUND":
            continue

        # Extract INFO fields
        info = {}

        for item in row["INFO"].split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                info[key] = value

        DP = int(info["DP"])
        AO = int(info["AO"])
        RO = int(info["RO"])

        # Extract gene ID from SnpEff ANN
        gene = ""

        ann = info.get("ANN", "")

        for annotation in ann.split(","):
            ann_fields = annotation.split("|")

            if len(ann_fields) > 3:
                if ann_fields[3] == "HCON_00005260":
                    gene = ann_fields[3]
                    break

        # Calculate resistance allele frequency
        frequency = AO / (AO + RO) * 100

        writer.writerow({
            "Sample_No": row["Sample_No"],
            "Sample_Name": row["Sample_Name"],
            "Mutation": row["Mutation"],
            "CHR": row["CHROM"],
            "POS": row["POS"],
            "Gene": gene,
            "REF": row["REF"],
            "ALT": row["ALT"],
            "DP": DP,
            "AO": AO,
            "RO": RO,
            "Frequency_percent": round(frequency, 2)
        })

print("Done!")
print("Output:", output_file)

```

