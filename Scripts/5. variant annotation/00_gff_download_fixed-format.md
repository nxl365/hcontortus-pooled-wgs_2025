# H.contortus gff.file —- downloaded and fixed format



## Download gff3 annotation file:

```bash
https://parasite.wormbase.org/Haemonchus_contortus_prjeb506/Info/Index/  

haemonchus_contortus.PRJEB506.WBPS19.annotations.gff3

```


## Fixed the error in format

When we built the H. contortus database for SnpEff, we encountered an error because the GFF3 file had only eight columns, while the standard format requires nine. We therefore preprocessed the file to fix this issue.  

### Create the python script:

```bash
00.preprocess_gff.ipynb
```

```python


file_path = "../../../../../data/haemonchus_contortus.PRJEB506.WBPS19.annotations.gff3"

out_gff = "../../../../../data/fixed-haemonchus_contortus.PRJEB506.WBPS19.annotations.gff3"  # Output GFF after filling empty 9th fields
```

```python
#### check how many lines have each column count: 9 or non-9 cloumns


col_counts = {}                                # dictionary: store how many lines have N columns
examples = []                                  # list: store a few examples of lines with non-9 columns

with open(file_path) as f:                     
    for i, line in enumerate(f, start=1):      # iterate lines, keeping track of line numbers
        if line.startswith('#') or not line.strip():             # Skip comment lines or blank lines
            continue
        parts = line.rstrip('\n').split('\t')     # remove any trailing newline, split the line into fields by tab \t
        n = len(parts)
        col_counts[n] = col_counts.get(n, 0) + 1  # set the associated number of key: col_counts.get(n, 0), give the value if the there is a current value associated with n, give 0 if not
        if n != 9 and len(examples) < 5:          # If it’s not 9 columns and we need no more than 5 examples
            examples.append((                     # Save (line number, column count, content)
                i,
                n,
                line.strip()
            ))

# print out how many lines have each column count
print("How many lines have each column count:")
for n in sorted(col_counts):
    print(f"  {n} columns: {col_counts[n]} lines")

# print out 5 examples of lines with non-9 columns
print("\nExamples of non-9-column lines:")
for line_num, n, content in examples:
    print(f"  Line {line_num} has {n} columns: {content}")


##########################################################
#### result: all lines in this gff file have 9 columns!! did not find any lines with non-9 columns, WHY???
```


```python
##### Check what is in the ninth column of the lines that caused errors in SnpEff running :
#" FATAL ERROR: Offending line (lineNum: 94893): #'hcontortus_chr4_Celeg_TT_arrow_pilon dust    low_complexity_region   1       4353    .       -              .       '
#java.lang.ArrayIndexOutOfBoundsException: Index 8 out of bounds for length 8" 

with open(file_path) as f:
    for i, line in enumerate(f, start=1):
        if "hcontortus_chr4_Celeg_TT_arrow_pilon" in line and "low_complexity_region" in line:
            print(f"Line {i}: {line.strip()}")    # error line content

            fields = line.strip('\n').split('\t')      # Split the line into fields by tab
            print(f"Line {i}: {len(fields)} columns")  # Print number of columns in the line

            for idx, fld in enumerate(fields, start=1):
                print(f"  Column {idx}: {repr(fld)}")   # Print each column's content
            
            col9 = fields[8]                            # Get the 9th column details
            print(f"Length of column 9: {len(col9)}")
            print(f"Unicode code point of column 9 :", [ord(c) for c in col9])  # Print code points of each character in column 9
            break  # Stop after the first match for demonstration


##########################################################
#### result: 
#### in the GFF file, the 9th column of error line is expected to contain attributes, but in this case, it is empty; since the gff file provide an tab \t after the 8th column, after we splite by tab, it is considered as an empty 9th column;
#### so in the gff file, there are only 8 columns in this error line
```



```python
#### replaces empty 9th col of some lines with a '.'
with open(file_path) as fin, open(out_gff, 'w') as fout:
    for line in fin:
        # copy comments and blank lines unchanged
        if line.startswith('#') or not line.strip():
            fout.write(line)
            continue

        # Split into 9 parts
        parts = line.strip('\n').split('\t')
        # If the 9th column is an empty string, replace it with '.'
        if parts[8] == '':
            parts[8] = '.'
        # Write the fixed line
        fout.write('\t'.join(parts) + '\n')

print(f"Done: Filled empty 9th columns in {out_gff}")
```



```python
#### checked the fixed gff: how many empty 9th in their line: should be zero now
count_empty_col9 = 0

with open(out_gff) as f:
    for i, line in enumerate(f, start=1):
        if "hcontortus_chr4_Celeg_TT_arrow_pilon" in line and "low_complexity_region" in line:
            #print(f"Line {i}: {line.strip()}")    # error line content

            fields = line.strip('\n').split('\t')      # Split the line into fields by tab
            #print(f"Line {i}: {len(fields)} columns")  # Print number of columns in the line

            #for idx, fld in enumerate(fields, start=1):
                #print(f"  Column {idx}: {repr(fld)}")   # Print each column's content
            
            # check the 9th column details
            #col9 = fields[8]                           
            #print(f"Length of column 9: {len(col9)}")
            #print(f"Unicode code point of column 9 :", [ord(c) for c in col9])  # Print code points of each character in column 9
            #break  # Stop after the first match for demonstration

            if fields[8] == '':
                count_empty_col9 += 1
            

print(f"How many lines with empty 9th column in total: {count_empty_col9}") 

```



## input
```bash
.../data/
`haemonchus_contortus.PRJEB506.WBPS19.annotations.gff3`
```

## output:
* This is uploaded in the `Data` folder.
```bash
.../data/
`fixed-haemonchus_contortus.PRJEB506.WBPS19.annotations.gff3`
```