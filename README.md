# README

## Description

Data and scripts associated with the article "Genomics of altitude-associated wing shape in two tropical butterflies" (2021) Molecular Ecology

## Contents

Most data files required to run the analyses are found within 1.data/. Large output files are hosted in dropbox (listed below in the tree but some not uploaded to github, others hosted by git LFS), but can be obtained by running the pipelines in 2.pipelines/.

Dropbox links for large files (optional):
[1.data/association/erato.gwas.shape.cov3.lrt0.gz.10snp.10snp.pos.ranks.short](https://www.dropbox.com/s/qkn777ap6069fzy/erato.gwas.shape.cov3.lrt0.gz.10snp.10snp.pos.ranks.short?dl=0)
[1.data/association/erato.gwas.shape.cov3.lrt0.gz.50snp.10snp.pos.ranks.short](https://www.dropbox.com/s/hfpi0tlex7fixy7/erato.gwas.shape.cov3.lrt0.gz.50snp.10snp.pos.ranks.short?dl=0)
[1.data/association/melpomene.gwas.shape.cov3.lrt0.gz.10snp.10snp.pos.ranks.short](https://www.dropbox.com/s/aif3c6y5jlfpj3m/melpomene.gwas.shape.cov3.doasso6.10snp.10snp.pos.ranks.short?dl=0)
[1.data/association/melpomene.gwas.shape.cov3.lrt0.gz.50snp.10snp.pos.ranks.short](https://www.dropbox.com/s/ge4msjzug9ljz20/melpomene.gwas.shape.cov3.doasso6.50snp.10snp.pos.ranks.short?dl=0)
[1.data/gene.set.anno/Heliconius_erato_demophoon_v1_-_proteins.interproscan.goterms.fasta.tsv](https://www.dropbox.com/s/s9dv8kyz4e8sq87/Heliconius_erato_demophoon_v1_-_proteins.interproscan.goterms.fasta.tsv?dl=0)
[1.data/gene.set.anno/Heliconius_melpomene_melpomene_Hmel2.5.proteins.interproscan.goterms.fasta.tsv](https://www.dropbox.com/s/yb392nodcvxo6ea/Heliconius_melpomene_melpomene_Hmel2.5.proteins.interproscan.goterms.fasta.tsv?dl=0)

GFF3 can be obtained from Lepbase ([link](http://download.lepbase.org/v4/features/))

In the 2.pipelines/ section you can find thorough explanations, with chunks of code, for carrying out GWAS/permutations and for obtaining wing morphology data from images. The outputs of these pipelines, which are used for plotting and analyses, are provided in the 1.data/ section.

The 3.scripts/ directory contains scripts that can be run to produce the figures and results presented in the manuscript and supplementary materials.


```
├── 1.data
│ ├── 1.master.era.mel.wild.reared.csv
│ ├── 2.haplotagging.samples.all.csv
│ ├── 3.mothers.rearing.for.SI.csv
│ ├── 4.mothers.haplotagging.formap.csv
│ ├── association
│ │ ├── era.479
│ │ │ ├── aspect.ratio.txt
│ │ │ ├── cov.txt
│ │ │ ├── cov3.top1.rows.lrtperm.txt
│ │ │ ├── era.n479.covariates.csv
│ │ │ ├── erato.gwas.shape.cov3.lrt0.gz.10snp.10snp.pos.ranks.short
│ │ │ ├── erato.gwas.shape.cov3.lrt0.gz.10snp.10snp.pos.ranks.short.outliers
│ │ │ ├── erato.gwas.shape.cov3.lrt0.gz.50snp.10snp.pos.ranks.short
│ │ │ └── erato.gwas.shape.cov3.lrt0.gz.50snp.10snp.pos.ranks.short.outliers
│ │ └── mel.187
│ │ ├── aspect.ratio.txt
│ │ ├── cov.txt
│ │ ├── cov3.top1.rows.lrtperm.txt
│ │ ├── mel.n187.covariates.csv
│ │ ├── melpomene.gwas.shape.cov3.doasso6.10snp.10snp.pos.ranks.short
│ │ ├── melpomene.gwas.shape.cov3.doasso6.10snp.10snp.pos.ranks.short.outliers
│ │ ├── melpomene.gwas.shape.cov3.doasso6.50snp.10snp.pos.ranks.short
│ │ └── melpomene.gwas.shape.cov3.doasso6.50snp.10snp.pos.ranks.short.outliers
│ ├── gene.set.anno
│ │ ├── Heliconius_erato_demophoon_v1_-_proteins.interproscan.goterms.fasta.tsv
│ │ └── Heliconius_melpomene_melpomene_Hmel2.5.proteins.interproscan.goterms.fasta.tsv
│ └── ld.decay
│ ├── era.vcf.thinned.22318kb.minr2.0.05.minuschr02.ld.summary
│ └── mel.vcf.thinned.17205kb.minr2.0.2.summary
├── 2.pipelines
│ ├── gwas.permutations.md
│ ├── pictures.to.wing.data.fiji.md
│ └── pictures.to.wing.data.fiji.pdf
├── 3.scripts
│ ├── 1.map.R
│ ├── 2.common.garden.plots.analyses.R
│ ├── 3.ngs.admix.plot.R
│ ├── 4.GO.R
│ ├── 5.ld.decay.R
│ └── extra
│ ├── compute50SNP10SNP.windows.R
│ └── perm.ranking.R
└── README.md
```




