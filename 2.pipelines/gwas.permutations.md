# GWAS (ANGSD), permutations

This pipeline outlines the necessary steps to run the GWAS and permutations.
The code below is not runnable in its current state, as the size of VCF files will require running the pipeline on a computing cluster, which vary widely in specifications.

## Input

VCFs for 666 individuals from Meier et al 2020

For each species there are two important files (which can be found in ./1.data/assocation/), which have as many rows as individuals in the sample: 

* cov.txt - covariates for association (wing area, sex code (0=male, 1=female), admixture proportion) 
* aspect.ratio.txt - aspect ratio (elongatedness)

Make sure VCF order of individuals matches covariates order (I always use ascending order of unit/cam IDs). Check order of vcf with bcftools query -l x.vcf

## 1. NGSadmix

Obtain admixture proportions to control for population structure in the GWAS. Install NGSadmix from [here](http://www.popgen.dk/software/index.php/NgsAdmix), 

```
# by Joana Meier
module load htslib/1.2.1 bcftools-1.9-gcc-5.4.0-b2hdt5n samtools/1.3.1 vcftools

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

Current version has fixed vcf reading [angsd version: 0.933-101-g7ea6e4b (htslib: 1.10.2-131-g0456cec) build(Aug 24 2020 10:35:22)]- however, you have to use the option to **-vcf-pl ** to specify we have PL likelihoods in our vcfs.


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
Rscript ./3.scripts/compute50SNP10SNP.windows.R out/${FILE}.gwas.shape.cov.lrt0.gz

# output will be "erato.gwas.shape.cov.lrt0.gz.50snp.10snp.pos"
```

The script 3.scripts/compute50SNP10SNP.windows.R takes the output from ANGSD, calculates p.values from likelihood ratio test statistic (chi2 distribution, 1df) and creates sliding windows of 50SNP and 10SNP steps (fixed) with the R package winscanr ([github](https://github.com/tavareshugo/WindowScanR)), to obtain median min, max, and median values from two columns: LRT.pval (p.values per snp) and position (to obtain mininum and maximum position)


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

# Compute window medians, min p-values
Rscript ./3.scripts/compute50SNP10SNP.windows.R out/perm_${PERM}.gwas.shape.cov.lrt0.gz


# to run array job
sbatch --array 1-200 perm.assoc.sh

```



### Extract outlier windows from permutations

Only extract permutation p-values from a list of outlier windows/SNPs  from the observed data GWAS (top 1%, rather than working with the whole genome, as it would be too much data to handle locally. This is concatenated into a single file with the 200 permutations as columns, and the outlier windows as rows. 

Obtain in R row numbers that correspond to the windows with the top 1% associtions (i.e. p-values above 99th percentile)

```
era.shape.gwas <- read.table("./1.data/assocation/era.479/erato.gwas.shape.cov.lrt0.gz.50snp.10snp.pos", row.names = NULL, header = F )
# add correct row number (which in the cluster includes the title)
era.shape.gwas$row.number<- 2:(nrow(era.shape.gwas)+1); tail(era.shape.gwas$row.number)

# get top 1% associated windows (lowest median p-values)
era.shape.top1.lrt <- subset(subset(era.shape.gwas), LRT.pval_median < quantile(na.rm = TRUE,LRT.pval_median, prob = 1 - 99/100)); head(era.shape.top1.lrt)

#make list of rows that correspond to highest lrt. will use these to extract from permutations
rows.to.extract <- era.shape.top1.lrt$row.number; tail(rows.to.extract ); head(rows.to.extract)
write.table(rows.to.extract, "era.n479/shape.gwas.doasso6.cov3.50snp.10snp.top1.rows.txt", row.names = F, col.names = F)


```

Then extract from each permutation, the median p-value corresponding to those windows (rows)

```
# get top1% rows for lrt.pval_median (CAREFUL COLUMN NUMBER, here column 10) for each permutation
for i in perm_*.gwas.shape.cov.lrt0.gz.50snp.10snp.pos
do
n=${i:5:3}
awk 'NR == FNR{a[$0]; next};FNR in a {print $10}' shape.gwas.doasso6.cov3.50snp.10snp.top1.rows.txt $i | column -t > lrt.mean.perm$n.txt
done

# concatenate top rows of all permutations
# copy to keep original, add the rest
cp shape.gwas.doasso6.cov3.50snp.10snp.top1.rows.txt top1.rows.lrtperm.txt
`paste ``-``d ``'\t'`` `` top1``.``rows``.``lrtperm``.``txt`` lrt``.``mean``.``perm``*.``txt ``>`` top1``.``rows``.``lrtperm1``.``txt``  ``;`` mv  top1``.``rows``.``lrtperm1``.``txt `` cov3``.``top1``.``rows``.``lrtperm``.``txt`


```



### Rank observed p-value among empirical (permutated) p-values

Now in R again, import the median p-values from 200 permutaitons (200 values per outlier window), add the observed values (real p-values), and rank them. The observed value ranking will determine where in the null distribution it falls, we will only consider significant those above the 99th percentile (i.e. with a ranking of 1/200 or 2/200). 

```
perm.era <- read.table("cov3.top1.rows.lrtperm.txt", blank.lines.skip = T)

## give permutations sensible names
colnames(perm.era)<- c("row", 1:200); head(perm.era)

# collapse permutations long format
perm.era.long <-gather(perm.era, perm, LRT_mean,-row)  ; names(perm.era.long)

# create df from observed data to add to permutations
to.add <- data.frame(row=era.shape.top1.lrt[,27], perm="observed", LRT_mean=era.shape.top1.lrt[,7]);names(to.add)

# add
perm.obvs.era.long <- rbind(perm.era.long,to.add)

# rank perm/obvserved within windows (i.e. row)
perm.obvs.era.long.rank <- perm.obvs.era.long %>%
  dplyr::group_by(row) %>%
  dplyr::mutate(rank = order(order(LRT_mean, decreasing=TRUE)))
obvs.era.long.rank <- subset(perm.obvs.era.long.rank, perm=="observed")

# add ranking to main df for plotting
era.shape.gwas$outlier.ranking.permutation200 <- obvs.era.long.rank$rank[match(era.shape.gwas$row.number, obvs.era.long.rank$row)]

```



### Alternative: extract lowest p-value per permutation

Grab minimum p-value per permutation (genome-wide), store. The 95th percentile will give a genome-wide threshold at p=0.05, presented in the SI Fig. S11.

Not very realistic for per site null distributions, as in every gwas of 25 million snps there will be some spuriously low p-values. 

```
# use p value min per window, column8, equivalent to finding minimum p-value genome wide (but faster)
for i in perm_*.gwas.shape.cov.lrt0.gz.50snp.10snp.pos; do awk 'NR==1 || $8 < min {line = $0; min = $8} END {print min}' $i >> min.LRT.pval.window.txt; done &

```

##  



