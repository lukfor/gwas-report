# gwas-report

> Create interactive html reports from GWAS summary statistics

**This script is now part of [genomic-utils](https://github.com/genepi/genomic-utils).**

## Installation

- Download gwas-report.jar from [releases](https://github.com/lukfor/gwas-report/releases)
- Validate installation with `java -jar gwas-report.jar version`

## Plot results from regenie

Create a html report from file `regenie.txt.gz` and store output in `regenie.report.html`:

```
java -jar gwas-report.jar report \
  --input regenie.txt.gz \
  --output regenie.report.html
```

If the columns in your file are different from regenie, you could set them manually:

```
java -jar gwas-report.jar report \
  --input regenie.txt.gz \
  --chr CHROM \
  --pos GENPOS \
  --pval LOG10P \
  --output regenie.report.html
```


## Annotate top hits

You could use [tabix-merge](https://github.com/lukfor/tabix-merge) to annotate results from regenie with rsids to use this information in summary table and manhattan plot:

```
java -jar gwas-report.jar report \
  --input regenie.annotated.txt.gz \
  --chr CHROM \
  --pos GENPOS \
  --pval LOG10P \
  --rsid RSID \
  --ref REF \
  --alt ALT \
  --beta BETA \
  --annotation TOP_HIT \
  --output regenie.report.html
```

You could use `--gene` to define which column contains the gene name and `--anotation GENE` to show it in the plot:

```
java -jar gwas-report.jar report \
  --input regenie.annotated.txt.gz \
  --chr CHROM \
  --pos GENPOS \
  --pval LOG10P \
  --rsid RSID \
  --ref REF \
  --alt ALT \
  --beta BETA \
  --gene NEAREST_GENE \
  --annotation TOP_HIT \
  --output regenie.report.html
```

## Create csv file with binned infos

Bins variants in `regenie.txt.gz` and store binned infos in `regenie.report.txt`. This file could be used to create plots for publications.

```
java -jar gwas-report.jar report \
  --input regenie.txt.gz \
  --chr CHROM \
  --pos GENPOS \
  --pval LOG10P \
  --format CSV \
  --output regenie.report.txt
```


## License

`gwas-report` is MIT Licensed.
