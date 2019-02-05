# `BCB420.2019.STRING`

#### (STRING data annotatation of human genes)

&nbsp;

###### Chantal Ho

----

## 1 About this package:

This package describes the workflow to download disease datasets from [the GEO database](https://www.ncbi.nlm.nih.gov/geo/), how to map the IDs to [HGNC](https://www.genenames.org/) symbols, and how to annotate the example gene set.

&nbsp;

#### In this project ...

```text
 --BCB420.2019.GEO/
   |__.gitignore
   |__.Rbuildignore
   |__BCB420.2019.GEO.Rproj
   |__DESCRIPTION
   |__dev/
      |__toBrowser.R               # display .md files in your browser
   |__inst/
      |__extdata/
         |__ilmn2hugo.RData         # ILMN ID to HGNC symbol mapping tool
      |__img/
         |__[...]                  # image sources for .md document
      |__scripts/
         |__recoverIDs.R           # utility to use biomaRt for ID mapping
   |__LICENSE
   |__NAMESPACE
   |__R/
   |__README.md                    # this file

```

&nbsp;

----

## Data download and cleanup

You can either download datasets manually from the websites:
1) From https://www.ncbi.nlm.nih.gov/gds/, search "disease and illumina" to get disease datasets with illumina ids
2) Click on the GSE id, which leads you to a webpage where you select Donwload Series Matrix TXT files
3) Put the file into the folder "../data"
4) Call getGEO

Or, download datasets directly using getGEO() on the GSE reference ID.

## 4 Mapping ENSEMBL IDs to HGNC symbols

STRING network nodes are Ensembl protein IDs. These can usually be mapped to HGNC symbols, but there might be ambiguities e.g. because alternatively spliced proteins might have different ENSP IDs that  map to the same HGNC symbol, or HGNC symbols have changed (they are frequently updated). To provide the best possible interpretation, we need to build a map of ENSP IDs to HGNC symbols. This requires care, because it is not guaranteed that all ENSP IDs can be mapped uniquely.** However, the usability of the dataset for annotation depends on the quality of this mapping.**

&nbsp;

#### Preparations: packages, functions, files

To begin, we need to make sure the required packages are installed:

**`readr`** provides functions to read data which are particularly suitable for
large datasets. They are much faster than the built-in read.csv() etc. But caution: these functions return "tibbles", not data frames. ([Know the difference](https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html).)
```R
if (! requireNamespace("readr")) {
  install.packages("readr")
}
```

**`biomaRt`** biomaRt is a Bioconductor package that implements the RESTful API of biomart,
the annotation framwork for model organism genomes at the EBI. It is a Bioconductor package, and as such it needs to be loaded via the **`BiocManager`**,
&nbsp;

```R
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
```
