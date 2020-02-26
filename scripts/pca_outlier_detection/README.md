# Automatic outlier detection

### Requirements
R installation

### Usage
Usage:

For latest options run scripts/classify_outliers.R -h

```
scripts/classify_outliers.R -f data/PCA.eigenvec -s data/sample_annot.tsv -e data/PCA.eigenval -o test
```

Where:
- -f eigenvectors without headers as output by plink
- -s optional sample annotation file with header and two columns. First column must contain the invidual ID's and second column is arbitrary sample annotation
- -e optionl list of eigenvalues for PCs on in each line
- -o output prefix
- -n number of pcs to take in to account in outlier detection
- -p if set then skips plotting of the PCAs
