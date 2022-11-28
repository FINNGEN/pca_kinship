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

The pipeline takes as an input a list of pruned variants from the imputation panel and then proceeds to build a plink file with only the selected snps. After this, king is run with flags --related --degree 3 and --duplicate. This returns the list of degree 3 couples as well as list of duplicate couples. Next, a pedigree fam file is built using the registry data (`--pheno-file`). `meta` is a tsv file with ID and a tag (in our case release number) that is used for plotting in order to find eventual anomalies in the connectivity within the batch.

```
usage: ped.py [-h] [--extract EXTRACT] [--fam FAM] --pheno-file PHENO_FILE
              [--meta META] --bed BED -o OUT_PATH --prefix PREFIX [--force]
              [--release] [--plot]
              [-log {critical,error,warn,warning,info,debug}]

kinship analysis & pedigree

optional arguments:
  -h, --help            show this help message and exit
  --extract EXTRACT     Path to list of variants to include
  --fam FAM             Optional .fam file to subset individuals
  --pheno-file PHENO_FILE
                        Phenotype filepath. Needs to contain SEX column
  --meta META           File with batch info fo samples
  --bed BED             BED filepath
  -o OUT_PATH, --out-path OUT_PATH
                        folder in which to save the results
  --prefix PREFIX       Output prefix
  --force               Flag on whether to force run
  --release             Flag to structure output for release
  --plot                Flag to plot
  -log {critical,error,warn,warning,info,debug}, --log {critical,error,warn,warning,info,debug}
                        Provide logging level. Example --log debug',
                        default='warning'

```

First a subset is created if `FAM` is provided. Then KING with --related --degree 3 and --duplicate flags is run.
This step also produces a plot of the distribution of kinship values.

Then KING with --build option is run. This allows to create a new .fam file with updated sex and FIDs.

## pca.py
This is script runs the pca analysis used in FG.

```
usage: pca.py [-h] --bed BED [--tg-bed TG_BED] --sample-info SAMPLE_INFO
              --meta META [--degree DEGREE] [--kin KIN]
              [--pca-components PCA_COMPONENTS] [--pc-filter PC_FILTER]
              [--finn-prob-filter FINN_PROB_FILTER] -o OUT_PATH [--name NAME]
              [-log {critical,error,warn,warning,info,debug}] [--force]
              [--test] [--cpus CPUS] [--release]

FinnGen PCA pipeline.

optional arguments:
  -h, --help            show this help message and exit
  --bed BED             Folder in which the merged plink file is stored
  --tg-bed TG_BED       Plink 1k file
  --sample-info SAMPLE_INFO
                        Path to csv file with sample,batch
  --meta META           Path to file with regionofbirth info
  --degree DEGREE       Degree for Kinship
  --kin KIN             File with king related individuals
  --pca-components PCA_COMPONENTS
                        Components needed for pca
  --pc-filter PC_FILTER
                        Number of pcs on which to perform the outlier
                        detection method
  --finn-prob-filter FINN_PROB_FILTER
                        Filter falue to decide whether a finngen sample is
                        part of the EUR or FIN centroid
  -o OUT_PATH, --out_path OUT_PATH
                        Folder in which to save the results
  --name NAME           Name to append to output files
  -log {critical,error,warn,warning,info,debug}, --log {critical,error,warn,warning,info,debug}
                        Provide logging level. Example --log debug',
                        default='warning'
  --force               Replaces files by force
  --test                Flag for quick pca_outlier method without plots.
  --cpus CPUS           Number of cpus to use (default available cpus)
  --release             Flag for data release

```

`sample-info` is a tsv file with sample ID and a tag that defines it (batch,cohort,release). It's used for plotting in order to try to identify batch/release effects in the PCA.

`meta` instead contains info about the region of birth. It's specific to Finland and it allows us to visualize pcas on the territory, where we expect to find a very significant correlation.