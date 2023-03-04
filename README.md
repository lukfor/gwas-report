# gwas-report

> Create interactive html reports from GWAS summary statistics


## Installation

- Download gwas-report.jar from [releases](https://github.com/lukfor/gwas-report/releases)
- Validate installation with `java -jar gwas-report.jar version`

## Plot results from regenie

Create a html report from file `regenie.txt.gz` and store output in `regenie.report.html`:

```
java -jar gwas-report.jar report \
  --input regenie.txt.gz \
  --chr CHROM \
  --pos GENPOS \
  --pval LOG10P \
  --output regenie.report.html
```

## Use Annotations

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
  --beta \
  --output regenie.report.html
```

## Create csv file with binned infos

Bins variants in `regenie.txt.gz` and store binned infos in `regenie.report.txt`:

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
