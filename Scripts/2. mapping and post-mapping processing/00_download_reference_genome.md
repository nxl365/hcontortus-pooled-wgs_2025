# Download the *H. contortus* reference genome

The *H. contortus* reference genome was downloaded from WormBase

```bash
cd .../raw_data/Hc_reference_genome

# Download the reference genome
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS19.genomic.fa.gz
```

## Note

`.fna` and `.fa` files are both FASTA-format files. The file extension is different, but the sequence content can be the same.

If a downstream analysis requires a `.fa` file, a `.fna` file can be renamed or copied as `.fa`.

```bash
cp file.fna file.fa  # keep both files

mv file.fna file.fa  # rename the file and keep only the .fa version
```

## Decompress the reference genome

```bash
# Decompress the reference genome from .fa.gz to .fa
gunzip -c file.fa.gz > file.fa

# Example used in this analysis
gunzip -c haemonchus_contortus.PRJEB506.WBPS19.genomic.fa.gz > haemonchus_contortus.PRJEB506.WBPS19.genomic.fa
```

## File size

```text
82M     haemonchus_contortus.PRJEB506.WBPS19.genomic.fa.gz
275M    haemonchus_contortus.PRJEB506.WBPS19.genomic.fa
```
