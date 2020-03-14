# pca_kinship
PCA and Kinship pipeline

## prune.py

Script that iteratively prunes a plink data set until a `target` number of snps is reached. 
```
usage: prune.py [-h] [--info FILE VALUE | --extract EXTRACT] --bed BED
                --out-path OUT_PATH --prefix PREFIX [--pargs PARGS] [--force]
                [--target TARGET] [--ld SIZE STEP THRESHOLD STEP2]

LD pruning of a bed file

optional arguments:
  -h, --help            show this help message and exit
  --info FILE VALUE     Info score file and filter value
  --extract EXTRACT     Path to list of variants to include
  --bed BED             BED filepath
  --out-path OUT_PATH   folder in which to save the results
  --prefix PREFIX       Output prefix
  --pargs PARGS         extra plink args
  --force               Flag on whether to force run
  --target TARGET       Target number of snps after pruning
  --ld SIZE STEP THRESHOLD STEP2    size,step,threshold,threshold_step

```

The algorithm will keep reducing the r2 parameter at `STEP2` intervals until either the number of snps is below target or it runs out of steps. At that point, the step that produced the closest number of snps to the target is returned and `target` number of snps are randomly returned.

## ped.py

The pipeline takes as an input a list of pruned variants from the imputation panel and then proceeds to build a plink file with only the selected snps. After this, king is run with flags --related --degree 3 and --duplicate. This returns the list of degree 3 couples as well as list of duplicate couples. Next, king is run with the flags --build --degree 3. This reconstructs the pedigree structure of the samples in the bed file.

```
usage: ped.py [-h] [--extract EXTRACT] [--fam FAM] --pheno-file PHENO_FILE
              --bed BED --out-path OUT_PATH --prefix PREFIX [--force]
              [--release]

kinship analysis & pedigree

optional arguments:
  -h, --help            show this help message and exit
  --extract EXTRACT     Path to list of variants to include
  --fam FAM             Optional .fam file to subset individuals
  --pheno-file PHENO_FILE  Phenotype filepath. Needs to contain SEX column
  --bed BED             BED filepath
  --out-path OUT_PATH   folder in which to save the results
  --prefix PREFIX       Output prefix
  --force               Flag on whether to force run
  --release             Flag to structure output for release

```

First a subset is created if `FAM` is provided. Then KING with --related --degree 3 and --duplicate flags is run.
This step also produces a plot of the distribution of kinship values.

Then KING with --build option is run. This allows to create a new .fam file with updated sex and FIDs.

## pca.py
