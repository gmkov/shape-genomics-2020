# README

## Description

Data and scripts associated with the article "Genomics of altitude-associated wing shape in two tropical butterflies"  (2020)


## Contents

(will be updated when final scripts are uploaded, and then released to Zenodo for a DOI)

In the 2.pipelines/ section you can find thorough explanations, with chunks of code, for carrying out GWAS/permutations and for obtaining wing morphology data from images. The outputs of these pipelines, which are used for plotting and analyses, are provided in the 1.data/ section.

The 3.scripts/ directory contains scripts that can be run to produce the figures and results presented in the manuscript and supplementary materials.

```
├── 1.data
│   ├── 1.master.era.mel.wild.reared.csv
│   ├── 2.haplotagging.samples.all.csv
│   ├── 3.mothers.rearing.for.SI.csv
│   └── association
│       ├── era.479
│       │   ├── aspect.ratio.txt
│       │   ├── cov.txt
│       │   └── era.n479.covariates.csv
│       └── mel.187
│           ├── aspect.ratio.txt
│           ├── cov.txt
│           └── mel.n187.covariates.csv
├── 2.pipelines
│   ├── gwas.permutations.md
│   ├── pictures.to.wing.data.fiji.md
│   └── pictures.to.wing.data.fiji.pdf
├── 3.scripts
│   └── compute50SNP10SNP.windows.R
└── README.md

```

