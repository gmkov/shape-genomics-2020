##!Rscript --vanilla
# Load packages required
require("data.table")
require("devtools")
require("windowscanr")
library(dplyr)
library(tidyr)

# Read in the vcf file name
args<-commandArgs(TRUE)
perm<-args[1]
obvs <- args[2]

#load
perm <- fread(perm, header=T, sep = '\t', stringsAsFactors=F)
gwas <- fread(obvs, header=F, sep = '\t', stringsAsFactors=F)

# reanme gwas
colnames(gwas) <- c("row.per.scaff","scaff", "win_start","win_end", "win_mean", "LRT.pval_n","LRT.pval_mean","LRT.pval_min" ,"LRT.pval_max","LRT.pval_median" , "beta_n","beta_mean","beta_min","beta_max","beta_median","SE_n" ,"SE_mean" ,
 "SE_min","SE_max","SE_median" , "Position_n" ,"midPos" , "Position_min" , "Position_max", "Position_median"); head(gwas)

# rename perm
perm <- as.data.frame(perm)
colnames(perm)<- c(1:200); head(perm)

# check
nrow(perm); nrow(gwas)

# add correct row number (which in the cluster includes the title)
gwas$row.number<- 2:(nrow(gwas)+1); tail(gwas$row.number)
perm$row.number<- 2:(nrow(perm)+1); tail(perm$row.number)

# collapse permutations long format
perm.long <-gather(perm, perm, LRT.pval_median, -row.number) ; head(perm.long)

# create df from observed data to add to permutations
to.add <- data.frame(row.number=gwas[,c("row.number")], perm="observed", LRT.pval_median=gwas[,c("LRT.pval_median")]);head(to.add)

# add
perm.obvs.long ← rbind(perm.long, to.add)

# rank perm/obvserved within windows (i.e. row)
perm.obvs.long.rank <- perm.obvs.long %>%
 dplyr::group_by(row.number) %>%
 dplyr::mutate(rank = order(order(LRT.pval_median, decreasing=FALSE))); head(perm.obvs.long.rank)

obvs.long.rank <- subset(perm.obvs.long.rank, perm=="observed"); head(obvs.long.rank)

# add ranking to main df for plotting
gwas$outlier.ranking.permutation200 <- obvs.long.rank$rank[match(gwas$row.number, obvs.long.rank$row.number)]
unique(gwas$outlier.ranking.permutation200)

# Write out results of ranking
write.table(gwas,file=paste0(obvs,".ranks"), sep = '\t',quote=FALSE)
