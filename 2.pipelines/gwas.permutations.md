# GWAS (ANGSD) and permutations

This pipeline outlines the necessary steps to run the GWAS and permutations.
The code below is not runnable in its current state, as the size of VCF files will require running the pipeline on a computing cluster, which vary widely in specifications.

## Input

VCFs for 666 individuals from Meier et al 2020

For each species there are two important files (which can be found in ./1.data/association/) , which have as many rows as individuals in the sample: 

* cov.txt - covariates for association (wing area, sex code (0=male, 1=female), admixture proportion) 
* aspect.ratio.txt - aspect ratio (elongatedness)

Make sure VCF order of individuals matches covariates order (I always use ascending order of unit/cam IDs). Check order of vcf with bcftools query -l x.vcf

## 1. NGSadmix

Obtain admixture proportions to control for population structure in the GWAS. Install NGSadmix from [here](http://www.popgen.dk/software/index.php/NgsAdmix), 

```
module load htslib/1.2.1 bcftools-1.9-gcc-5.4.0-b2hdt5n samtools/1.3.1 vcftools

# by Joana Meier: 
# replace with your file prefix
file=erato
# Fix the missing genotypes (they need a PL)
zcat $file.vcf.gz | sed -e 's|./.:.:.:.:.:.|0/0:./.:0.33,0.33,0.34:0:10,10,10:0,0|g' | gzip > $file.corr.vcf.gz

# Generate a Beagle file
angsd -out $file -nThreads 8 \
  -vcf-gl $file.corr.vcf.gz -doGlf 2 \
  -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 

# Run NGSadmix
nice NGSadmix -likes $file.beagle.gz -K 2 -minMaf 0.05 -outfiles $file.maf0.05.K2 

# Get 1 column for GWAS, which is added to cov.txt file
cut -f 1 -d " " $file.maf0.05.K2.qopt > NGSadmix.prop

```

## 2. ANGSD association

Current version of ANGSD has fixed vcf reading [angsd version: 0.933-101-g7ea6e4b (htslib: 1.10.2-131-g0456cec) build(Aug 24 2020 10:35:22)]- however, you have to use the option to **-vcf-pl** to specify we have PL likelihoods in our vcfs.


```
module load python-3.5.2-gcc-5.4.0-rdp6q5l gcc
module load R

REF=~/Heliconius_erato_demophoon_v1_-_scaffolds.fa
FILE=erato

angsd -yQuant ./1.data/assocation/era.479/aspect.ratio.txt \
 -doAsso 6 -doPost 1 -out out/${FILE}.gwas.shape.cov \
 -doMajorMinor 4 -doMaf 1 -vcf-pl ./erato.vcf.gz \
 -P 8 -model 1 -cov ./1.data/assocation/era.479/cov.txt \
 -ref $REF

# Compute window median, min p-values (see short explanation below)
Rscript ./3.scripts/extra/compute50SNP10SNP.windows.R out/${FILE}.gwas.shape.cov.lrt0.gz

# output will be "erato.gwas.shape.cov.lrt0.gz.50snp.10snp.pos"

```

The script 3.scripts/extra/compute50SNP10SNP.windows.R takes the output from ANGSD, calculates p.values from likelihood ratio test statistic (chi2 distribution, 1df) and creates sliding windows of 50SNP and 10SNP steps (fixed) with the R package winscanr ([](https://github.com/tavareshugo/WindowScanR)[tavareshugogithub](https://github.com/tavareshugo)), to obtain min, max, and median values from two columns: LRT.pval (p.values per snp) and position (to obtain mininum and maximum position)


## 3. Genome-wide permutations

Create 200 random “phenotype files” (i.e. aspect.ratio.txt) by randomly permuting the observed aspect ratio values. I use R:

```
mel.cov <-  read.csv("./1.data/assocation/mel.187/mel.n187.covariates.csv")
era.cov <-  read.csv("./1.data/assocation/era.479/era.n479.covariates.csv")

mel.ar <- mel.cov$aspect.ratio
era.ar <- era.cov$aspect.ratio

era.rand.ar <-list()
for(i in 1:200){
  era.rand.ar[[i]] <- sample(era.ar) }

for(i in 1:200){
  write.table(era.rand.ar[[i]], paste0("permutations/era.n479/aspect.ratio.", str_pad(i, 3, pad = "0"),".txt" ), row.names = F,col.names = F ) }

mel.rand.ar <-list()
for(i in 1:200){
  mel.rand.ar[[i]] <- sample(mel.ar) }

for(i in 1:200){
  write.table(mel.rand.ar[[i]], paste0("permutations/mel.n187/aspect.ratio.", str_pad(i, 3, pad = "0"),".txt" ), row.names = F,col.names = F ) }

```

This will give rise to 200 files per species, called “permutations/mel.n187/aspect.ratio.[001-200].txt”. These will be used as input for ANGSD.

On a cluster. Array job running GWAS 200 times per species, one for each permuted phenotype file
R script computed window averages for 50SNP win 10SNP step (same as observed data)

```
# create script perm.assoc.sh
module load python-3.5.2-gcc-5.4.0-rdp6q5l gcc
REF=~/Heliconius_erato_demophoon_v1_-_scaffolds.fa
PERM=`sed -n -e "$SLURM_ARRAY_TASK_ID p" permutations/mel.n187/aspect.ratio.`

angsd -yQuant permutations/mel.n187/aspect.ratio.${PERM}.txt \
-doAsso 6 -doPost 1 -out out/perm_${PERM}.gwas.shape.cov \
-doMajorMinor 4 -doMaf 1 -vcf-pl ./erato.vcf.gz \
-P 8 -model 1 -cov ./1.data/assocation/era.479/cov.txt \
-ref $REF

# Compute window medians, min p-values, same as with observed data
Rscript ./3.scripts/extra/compute50SNP10SNP.windows.R out/perm_${PERM}.gwas.shape.cov.lrt0.gz

# to run array job
sbatch --array 1-200 perm.assoc.sh

```



### Rank observed p-values among empirical (permutated) p-values

First extract all median p-values of all windows and permutations.

Then the R script 3.scripts/extra/perm.ranking.R will import the median p-values from 200 permutations (200 values per outlier window), add the observed values (real p-values), and rank all p-values per window (observed +permuted, total=201). The observed value ranking will determine where in the null distribution it falls, we will only consider significant those above the 99th percentile (i.e. with a ranking of 1/201 or 2/201). 


```
`# extract the column (i.e. all windows/rows) with median p-values from each permutation (careful col names may be shifted)`
`for`` i ``in`` perm_``*.``gwas``.``shape``.``cov``.``lrt0``.``gz``.``50snp``.``10snp``.``pos`
`do`
`n``=``$``{``i``:``5``:``3``}`
`awk ``-``F ``'\t'`` ``'{ print $10}'`` $i ``|`` column ``-``t ``>`` lrt``.``median``.``all``.``rows``.``perm$n``.``txt`
`done`

`# concatenate all rows of all permutations`
`paste ``-``d ``'\t'``  lrt``.``median``.``all``.``rows``.``perm``*.``txt ``>`` cov3``.``all``.``rows``.``lrtperm``.``txt`

`# check correct number of columns, 200`
`awk ``-``F``'\t'`` ``'{print NF; exit}'`` cov3``.``all``.``rows``.``lrtperm``.``txt`

`# run r in a cluster with:`
`# input1- concatenated permutation results and`
`# input 2- observed gwas results`
`nice ``Rscript`` ``../../../../``scripts``/``perm``.``ranking``.``R  cov3``.``all``.``rows``.``lrtperm``.``txt erato.gwas.shape.cov.lrt0.gz.50snp.10snp.pos`` ``&>`` perm``.``ranking``.``out``&`


```




### For SI: extract lowest p-value per permutation

Grab minimum p-value per permutation (genome-wide), store. The 95th percentile will give a genome-wide threshold at p=0.05, presented in the SI Fig. S11.

Not very realistic for per site null distributions, as in every gwas of 25 million snps there will be some spuriously low p-values in each genome-wide permutation. 

```
# use p value min per window, column8, equivalent to finding minimum p-value genome wide (but faster)
for i in perm_*.gwas.shape.cov.lrt0.gz.50snp.10snp.pos; do awk 'NR==1 || $8 < min {line = $0; min = $8} END {print min}' $i >> min.LRT.pval.window.txt; done &

```

##  

## For SI Fig. S11: LD decay

First thin the full vcf snp dataset, so subsample snps that are at least XXbp apart → median size of our 50SNP windows which corresponds to:

* erato 1200bp
* melpomene 1032bp

```
`ERA_VCF``=erato``.``vcf``.``gz`
MEL_VCF=melpomene.vcf.gz

module load vcftools
vcftools --gzvcf $ERA_VCF --thin 1200 --out ./era.thinned1200  `--``recode `&
vcftools --gzvcf $MEL_VCF --thin 1032 --out ./mel.thinned1032  --recode &

```


Obtain r2 between snps with plink

```
# --ld-window-kb max size chr, so that comparisons always on diff chromosomes
# chr1scaff1 era 22318219
# chr1scaff1 mel 17205807
# min r2 values based on previous studies (Martin et al 2019, Van Belleghem et al 2017)
module load plink
`plink ``--``vcf`` ``era``.``thinned1200``.``recode``.``vcf``  ``--``r2 ``--``ld``-``window``-``kb ``22318`` ``--``ld``-``window``-``r2 ``0.05`` ``--``ld``-``window ``999999`` ``--``allow``-``extra``-``chr ``--``set``-``missing``-``var``-``ids ``@:#`` ``--``threads ``20`` ``--``out`` era``.``vcf``.``thinned``.``22318kb``.``minr2``.``0.05`
plink --vcf mel.thinned1032.recode.vcf  --r2 --ld-window-kb 17205 --ld-window-r2 0.2 --ld-window 999999 --allow-extra-chr --set-missing-var-ids @:# --threads 20 --out mel.vcf.thinned.17205kb.minr2.0.2

nohup nice bash ld.plink.era.sh &> ld.plink.era.sh.out & 
nohup nice bash ld.plink.mel.sh &> ld.plink.mel.sh.out & 

# remove lines that contain Herato02- large inversion skews ld decay
sed '/Herato02/d' era.vcf.thinned.22318kb.minr2.0.05.ld >era.vcf.thinned.22318kb.minr2.0.05.minuschr02.ld

# get out distance between snps
cat mel.vcf.thinned.17205kb.minr2.0.2.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > mel.vcf.thinned.17205kb.minr2.0.2.summary
cat era.vcf.thinned.22318kb.minr2.0.05.minuschr02.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > era.vcf.thinned.22318kb.minr2.0.05.minuschr02.ld.summary

```









