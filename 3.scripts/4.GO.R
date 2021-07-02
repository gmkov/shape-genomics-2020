##### GMK April 2019 #####
#### packages ####
rm(list=ls())
dev.off()

library(ggplot2)
library(dplyr)
library(tidyr)
library(qqman)
library(cowplot)
library(data.table)
library(stringr)
library(devtools)
library(ggman)
options(scipen = 999)
library(RcppCNPy)
library(GenomicRanges)
library(GenomicFeatures)
library(topGO)

##### functions ##### 



############################## 0. data prep #########################
# use 10 snp no overlaps to look at densities
era.shape.gwas <- read.table("1.data/association/era.479/erato.gwas.shape.cov3.lrt0.gz.10snp.10snp.pos.ranks.short.outliers", header = T)
mel.shape.gwas <- read.table("1.data/association/mel.187/melpomene.gwas.shape.cov3.doasso6.10snp.10snp.pos.ranks.short.outliers", header = T)

# references available at LepBase
ref.scaff.era <- read.table("./Heliconius_erato_demophoon_v1_-_scaffolds.fa.fai", row.names = NULL)
ref.scaff.mel <- read.table("./Hmel2.5.scaffolds.fa.fai", row.names = NULL)


############################## 1. get outlier ranges #########################
######## 1. era, prep outlier blocks ranges  ########
#### keep only outlier ranges with more than 10 outlier windows - i.e. clusters #####

# make granges, keep metadata, xxx blocks
head(era.shape.gwas)
era.all.summ.ranges <-  makeGRangesFromDataFrame(era.shape.gwas, keep.extra.columns = T, 
                                                 start.field="Position_min",end.field=c("Position_max"),
                         seqnames.field=c("scaff"),
                         ignore.strand =T)


######## 2. mel, prep outlier blocks ranges  ########
##### 2.2 keep only outlier ranges with more than 10 outlier windows - i.e. clusters #####

# make granges, keep metadata, xxx blocks
head(mel.shape.gwas)
mel.all.summ.ranges <-  makeGRangesFromDataFrame(mel.shape.gwas, keep.extra.columns = T, 
                                                 start.field="Position_min",end.field=c("Position_max"),
                                                 seqnames.field=c("scaff"),
                                                 ignore.strand =T)

############################## 2. extract genes (all and outlier) #########################
########### 2.1 prep gff ranges, extract ALL genes  ########
## get genes from gff with genomic ranges
# create txdb from gff https://www.biostars.org/p/167818/
txdb.era <- makeTxDbFromGFF("1.data/gene.set.anno/Heliconius_erato_demophoon_v1.gff3" , format="gff3")
txdb.mel <- makeTxDbFromGFF("1.data/gene.set.anno/Heliconius_melpomene_melpomene_Hmel2.5.gff3" , format="gff3")

# checks
exonsBy(txdb.era, by="gene")
mcols(GenomicFeatures::cds(txdb.era, columns=c("TXID", "TXNAME", "GENEID"))); GenomicFeatures::transcripts(txdb.era)
nrow(mcols(GenomicFeatures::genes(txdb.era, columns=c("TXID", "TXNAME", "GENEID")))) # total 13676 genes

# obtain gene lengths for all genes of each species, called "width"
era.gene.lengths <-as.data.frame(GenomicFeatures::transcripts(txdb.era))
mel.gene.lengths <-as.data.frame(GenomicFeatures::transcripts(txdb.mel))

########### 2.2 erato outlier genes   ########
length(GenomicFeatures::genes(txdb.era)) 
length(GenomicFeatures::genes(txdb.mel)) 


# use txname to overlap with ranges, as this has model instead of TU
overlap.genes.era <- subsetByOverlaps( GenomicFeatures::genes(txdb.era, columns=c("TXNAME")), era.all.summ.ranges, type="any")
head(as.data.frame(overlap.genes.era)); str(as.data.frame(overlap.genes.era$TXNAME))

# use values (as these will match real names) # 4k genes
genes.era <-  as.data.frame(overlap.genes.era$TXNAME)$value; length(genes.era); length(unique(genes.era))

########### 2.3 mel outlier genes   ########
length(GenomicFeatures::genes(txdb.mel)) # 20k genes

# use txname to overlap with ranges, as this has model instead of TU
overlap.genes.mel <- subsetByOverlaps( GenomicFeatures::genes(txdb.mel, columns=c("TXNAME")), mel.all.summ.ranges,type="any")

head(as.data.frame(overlap.genes.mel)); str(as.data.frame(overlap.genes.mel$TXNAME))
# use values (as these will match real names)
genes.mel <-  as.data.frame(overlap.genes.mel$TXNAME)$value; length(genes.mel)

############################## 3. go terms enrich #########################
############### 3.1 erato goterms  ####
library(topGO)

##### new (2021) interproscan goterm annotations #####
era.inter.new <- fread("1.data/gene.set.anno/Heliconius_erato_demophoon_v1_-_proteins.interproscan.goterms.fasta.tsv",header = F, sep = '\t' , fill = T); head(era.inter.new)

era.inter.new <-era.inter.new[,c(1,5)]
names(era.inter.new) <- c("gene.name", "go.term"); era.inter.new <-  as.data.frame(era.inter.new)
# change pipe for ", "
era.inter.new <- era.inter.new %>% mutate(go.term = str_replace_all(go.term, '\\|', ", ")); head(era.inter.new)

# remove entries without go terms
era.inter.new <- subset(era.inter.new, go.term!=""&go.term!="-"); head(era.inter.new)

# concatenate all go terms
era.inter.new.go <-  era.inter.new %>%
  group_by(gene.name)  %>%
  mutate(go.term = paste0(unique(go.term), collapse = ", ")) ; head(era.inter.new.go)
era.inter.new.go <- era.inter.new.go[!duplicated(era.inter.new.go$gene.name), ]; head(era.inter.new.go)

# subset string of go terms to keep only one
era.inter.new.go$go.term <- sapply(strsplit(era.inter.new.go$go.term, ",", fixed = TRUE), function(x)
  paste(unique(x), collapse = ",")); head(era.inter.new.go)

write.table(era.inter.new.go, "1.data/gene.set.anno/era.inter.new.unique.anno.txt", sep = '\t', quote = F, row.names = F, col.names = F)

##### topgo era outliers all ########
# load new interproscan erato GO annotation
era.geneID2GO <- topGO::readMappings(file="1.data/gene.set.anno/era.inter.new.unique.anno.txt")

# check how montgomery does it
era.geneID2GO <- topGO::readMappings(file="1.data/gene.set.anno/Heliconius_erato_demophoon_v1_-_proteins.interproscan.goterms.fasta.tsv", sep = "\t", header=F)
era.geneUniverse <- names(era.geneID2GO) 

## gene / topgodata prep #####
# grab all genes that overlap from outliers from above
era.genes <-  tibble(genes=genes.era); names(genes.era)<- genes.era; head(era.genes); nrow(unique(era.genes))

# unique, they were already unique but just in case
era.genesOfInterest <- as.character(unique(era.genes$genes)); length(era.genesOfInterest)
# indicate in the unverse whether gene present in list or not 
era.geneList <- factor(as.integer(era.geneUniverse %in% era.genesOfInterest))
names(era.geneList) <- era.geneUniverse; head(era.geneList)
sum(as.numeric(as.character(era.geneList))) #392 outlier genes have gene ontology

# select BP biological process ontology
era.myGOdata <- new("topGOdata", description="outliers", ontology="BP", allGenes=era.geneList,  
                    annot = annFUN.gene2GO, gene2GO = era.geneID2GO)
# important: feasible genes means that they have biological process ontology, they may have molecular function but not BP

## GOenrich #####
# use weight to account for hierarchy. or elim
result.weight.Fisher <- runTest(era.myGOdata, algorithm="weight01", statistic="fisher",node=10 )
hist(score(result.weight.Fisher), 50, xlab = "p-values")

era.allRes <- GenTable(era.myGOdata, weightFisher=result.weight.Fisher ,
                   orderBy = "weightFisher"  , 
                   ranksOf ="result.weight.Fisher", topNodes = 50)

# showSigOfNodes(era.myGOdata, score(result.weight.Fisher), firstSigNodes = 5, useInfo ='all')
printGraph(era.myGOdata, result.weight.Fisher, firstSigNodes = 5, fn.prefix = "1.data/gene.set.anno/tGO.era", 
           useInfo = "def", pdfSW = TRUE)

era.allRes$Term <- gsub(" [a-z]*\\.\\.\\.$", "", era.allRes$Term); era.allRes$Term <- gsub("\\.\\.\\.$", "", era.allRes$Term); era.allRes$Term
era.allRes$Term <- paste( era.allRes$Term, era.allRes$GO.ID,sep=", "); era.allRes$Term
era.allRes$Term <- factor(era.allRes$Term, levels=rev(era.allRes$Term)); 
era.allRes$weightFisher <- as.numeric(era.allRes$weightFisher)
era.allRes <- subset(era.allRes, weightFisher<0.05) # important, keep only significant

############### 3.2 melpomene goterms #####
##### new (2021) interproscan goterm annotations #####
mel.inter.new <- fread("1.data/gene.set.anno/Heliconius_melpomene_melpomene_Hmel2.5.proteins.interproscan.goterms.fasta.tsv",header = F, sep = '\t' , fill = T); head(mel.inter.new)
mel.inter.new <-mel.inter.new[,c(1,5)]; head(mel.inter.new)
names(mel.inter.new) <- c("gene.name", "go.term"); mel.inter.new <-  as.data.frame(mel.inter.new)
# change pipe for ", "
mel.inter.new <- mel.inter.new %>% mutate(go.term = str_replace_all(go.term, '\\|', ", ")); head(mel.inter.new)

# remove entries without go terms
mel.inter.new <- subset(mel.inter.new, go.term!=""&go.term!="-"); head(mel.inter.new)

# concatenate all go terms
mel.inter.new.go <-  mel.inter.new %>%
  group_by(gene.name)  %>%
  mutate(go.term = paste0(unique(go.term), collapse = ", ")) ; head(mel.inter.new.go)
mel.inter.new.go <- mel.inter.new.go[!duplicated(mel.inter.new.go$gene.name), ]; head(mel.inter.new.go)

# subset string of go terms to keep only one
mel.inter.new.go$go.term <- sapply(strsplit(mel.inter.new.go$go.term, ",", fixed = TRUE), function(x)
  paste(unique(x), collapse = ",")); head(mel.inter.new.go); nrow(mel.inter.new.go)

write.table(mel.inter.new.go, "1.data/gene.set.anno/mel.inter.new.unique.anno.txt", sep = '\t', quote = F, row.names = F, col.names=F)

##### topgo mel outliers all ########
# load all melto annotation
mel.geneID2GO <- topGO::readMappings(file="1.data/gene.set.anno/mel.inter.new.unique.anno.txt")
mel.geneUniverse <- names(mel.geneID2GO) ; head(mel.geneID2GO); length(mel.geneID2GO)

# outliers
mel.genes <-  tibble(genes=genes.mel); names(genes.mel)<- genes.mel; nrow(mel.genes)

# unique!!!
mel.genesOfInterest <- as.character(unique(mel.genes$genes)) ; length(mel.genesOfInterest)
mel.geneList <- factor(as.integer(mel.geneUniverse %in% mel.genesOfInterest))
sum(as.numeric(as.character(mel.geneList))) #261 outlier genes have gene ontology

length(mel.geneList) # total genes
names(mel.geneList) <- mel.geneUniverse 
mel.myGOdata <- new("topGOdata", description="outliers", ontology="BP", allGenes=mel.geneList,  annot = annFUN.gene2GO, gene2GO = mel.geneID2GO)

# use weight to account for himelrchy
result.weight.Fisher <- runTest(mel.myGOdata, algorithm="weight01", statistic="fisher", nodes=5)
mel.allRes <- GenTable(mel.myGOdata, weightFisher=result.weight.Fisher ,
                       orderBy = "weightFisher"  , 
                       ranksOf ="result.weight.Fisher", topNodes = 50)
mel.allRes <- subset(mel.allRes, weightFisher<0.05)
mel.allRes$Term <- gsub(" [a-z]*\\.\\.\\.$", "", mel.allRes$Term); mel.allRes$Term <- gsub("\\.\\.\\.$", "", mel.allRes$Term); 
mel.allRes$Term <- paste( mel.allRes$Term, mel.allRes$GO.ID,sep=", ")
mel.allRes$Term <- factor(mel.allRes$Term, levels=rev(mel.allRes$Term)); mel.allRes$weightFisher <- as.numeric(mel.allRes$weightFisher)

############### 3.3. era mel plot GOenrich #####
## era
era.allRes$in.melpomene <- era.allRes$Term %in%  mel.allRes$Term

era.go.plot <-ggplot(subset(era.allRes, weightFisher<0.05), aes(x=Term, y=-log10(weightFisher), fill=in.melpomene)) +
  geom_bar(stat="identity",size=.9,width = .8 )+
  #stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("Enrichment -log10(pvalue)") +
  scale_color_manual(values=c("transparent", "red"))+
  scale_fill_manual(values=c("darkgrey", "red"))+
  #ggtitle(expression(italic('H. erato')~ 'shared outlier windows'))+
  # scale_y_continuous(expand = c(0,0), breaks = round(seq(0,2.5, by = 1), 1)) +
  # scale_y_continuous(expand = c(0,0), breaks = round(seq(0, 4, by = 2), 1), limits = c(0,4)) +
  theme_bw()+ 
  scale_fill_manual(values=c("darkgrey", "red"))+
  # scale_fill_viridis_c()+ labs(fill="p-value")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length.y.right =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        # axis.text.y.left = element_blank(),
        # axis.text.y= element_text(face=ifelse(in.melpomene==TRUE,"bold","italic")),
        #panel.background = element_rect(fill="grey90", color="transparent"), 
        axis.text = element_text(size=10),  legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(1,0.5, 0.5,0.5), "cm"),
        axis.title=element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
        strip.background =element_rect(fill="transparent", color="transparent")) +coord_flip();era.go.plot
## mel
mel.allRes$in.erato <- mel.allRes$Term %in% era.allRes$Term; head(mel.allRes)

mel.go.plot <- ggplot(mel.allRes, aes(x=Term, y=-log10(weightFisher), fill=in.erato)) +
  geom_bar(stat="identity",size=.9,width = .8 )+
 # stat_summary(geom = "bar", fun = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("Enrichment") +
  scale_color_manual(values=c("transparent", "red"))+
  scale_fill_manual(values=c("darkgrey", "red"))+
  #ggtitle(expression(italic('H. melpomene')~ 'shared outlier windows'))+
  scale_y_continuous(expand = c(0,0), breaks = round(seq(0, 5, by = 2), 1), limits = c(0,5)) +
  theme_bw()+
  #scale_fill_viridis_c()+
  labs(fill="p-value")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),#axis.ticks.x=element_blank(),
        axis.ticks.length.y.right =unit(-0.15, "cm"), #axis.ticks.margin=unit(0.5, "cm"),
        # axis.text.y.left = element_blank(),
        #axis.text.y= element_text(face=ifelse(mel.allRes$in.erato==TRUE,"bold","italic")),
        #panel.background = element_rect(fill="grey90", color="transparent"), 
        axis.text = element_text(size=10), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(),# axis.text.x=element_blank(), #axis.text.y=element_blank(), #
        plot.margin = unit(c(1,0.5, 0.5,0.5), "cm"),
        axis.title=element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank(),# axis.title.x=element_blank(),#axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_blank(),
        strip.background =element_rect(fill="transparent", color="transparent")) +coord_flip(); mel.go.plot


plot_grid(era.go.plot, mel.go.plot, cols = 1, rel_heights = c(1.4,1), align = "hv")

############################## 4. SHARING SIMULATIONS go terms enrich #########################
########### 4.1 gene length prep ###########
# how long are genes in observed dataset of genes that have GO terms - era.geneList (from universe) that is present in gene set (==1)
subset(era.gene.lengths, tx_name %in% names(era.geneList[era.geneList==1]))$width


# how much sharing in observed dataset
era.mel.go.terms.sharing <- as.data.frame(tibble(perm.no="observed", no.outlier.terms.era=nrow(era.allRes), no.outlier.terms.mel=nrow(mel.allRes),
                                                     no.terms.shared.between.sp=sum(era.allRes$Term %in% mel.allRes$Term ),
                                                 no.era.outlier.genes=length(era.geneList[era.geneList==1]), 
                                                 mean.era.gene.length=mean(subset(era.gene.lengths, tx_name %in% names(era.geneList[era.geneList==1|era.geneList==0]))$width),
                                                 mean.era.outlier.gene.length=mean(subset(era.gene.lengths, tx_name %in% names(era.geneList[era.geneList==1]))$width),
                                                 no.mel.outlier.genes=length(mel.geneList[mel.geneList==1]), 
                                                 mean.mel.gene.length=mean(subset(mel.gene.lengths, tx_name %in% names(mel.geneList[mel.geneList==1|mel.geneList==0]))$width),
                                                 mean.mel.outlier.gene.length=mean(subset(mel.gene.lengths, tx_name %in% names(mel.geneList[mel.geneList==1]))$width),
                                                 go.terms.shared.era.mel=paste0(subset(era.allRes, Term %in% mel.allRes$Term )$Term, collapse = '|') )); era.mel.go.terms.sharing

########### 4.2 random outliers prep ###########
#### obtain positions of snp windows (all) & save ####
era.shape.gwas.long <- fread("1.data/association/era.479/erato.gwas.shape.cov3.lrt0.gz.10snp.10snp.pos.ranks.short")
summarise(group_by(era.shape.gwas.long,is.outlier), n=n()) #1.32% are outliers, 14665 oulier windows
era.shape.gwas.long <- era.shape.gwas.long[,c("scaff", "Position_min", "Position_max")]; head(era.shape.gwas.long)

mel.shape.gwas.long <- fread("1.data/association/mel.187/melpomene.gwas.shape.cov3.doasso6.10snp.10snp.pos.ranks.short")
summarise(group_by(mel.shape.gwas.long,is.outlier), n=n()) #1.0009% are outliers, 10573 oulier windows
mel.shape.gwas.long <- mel.shape.gwas.long[,c("scaff", "Position_min", "Position_max")]; head(mel.shape.gwas.long)


#### create random outlier datasets ####
era.shape.gwas.outlier.ran <- NA; mel.shape.gwas.outlier.ran <- NA
era.all.summ.ranges.ran <- NA; mel.all.summ.ranges.ran <- NA
overlap.genes.era.ran <- NA; overlap.genes.mel.ran <- NA
genes.era.ran <- NA; era.genesOfInterest.ran.list <- NA; era.geneList.ran.list <- NA
genes.mel.ran <- NA; mel.genesOfInterest.ran.list <- NA; mel.geneList.ran.list <- NA
era.myGOdata.ran.list <- NA; result.weight.Fisher.era.ran.list<-NA; era.allRes.ran.list <-NA
mel.myGOdata.ran.list <- NA; result.weight.Fisher.mel.ran.list<-NA; mel.allRes.ran.list <-NA

# to store
era.mel.go.terms.sharing.ran <- as.data.frame(tibble(perm.no=NA, no.outlier.terms.era=NA, no.outlier.terms.mel=NA,
                                                     no.terms.shared.between.sp= NA,
                                                     no.era.outlier.genes=NA, mean.era.gene.length=NA,mean.era.outlier.gene.length=NA,
                                                     no.mel.outlier.genes=NA,  mean.mel.gene.length=NA, mean.mel.outlier.gene.length=NA,
                                                     go.terms.shared.era.mel=NA)); era.mel.go.terms.sharing.ran

for (i in 1:1000){
  ###### random outliers #####
  # sample the same number of rows as outliers in observed dataset
  era.shape.gwas.outlier.ran <-  sample_n(era.shape.gwas.long, 14665, replace = F)
  # order by scaff, position
  era.shape.gwas.outlier.ran <- era.shape.gwas.outlier.ran[order( c(era.shape.gwas.outlier.ran[,1], era.shape.gwas.outlier.ran[,2])), ]
  # mel- sample the same number of rows as outliers in observed dataset
  mel.shape.gwas.outlier.ran <-  sample_n(mel.shape.gwas.long, 10573, replace = F)
  # order by scaff, position
  mel.shape.gwas.outlier.ran <- mel.shape.gwas.outlier.ran[order( c(mel.shape.gwas.outlier.ran[,1], mel.shape.gwas.outlier.ran[,2])), ]
  
  ###### create outlier ranges #####
  era.all.summ.ranges.ran <-  makeGRangesFromDataFrame(era.shape.gwas.outlier.ran, keep.extra.columns = F,
                                                            start.field="Position_min",end.field=c("Position_max"), seqnames.field=c("scaff"), ignore.strand =T)
  mel.all.summ.ranges.ran <-  makeGRangesFromDataFrame(mel.shape.gwas.outlier.ran, keep.extra.columns = F, 
                                                            start.field="Position_min",end.field=c("Position_max"), seqnames.field=c("scaff"), ignore.strand =T)
  ###### find gene overlaps with ranges #####
  overlap.genes.era.ran <- subsetByOverlaps( GenomicFeatures::genes(txdb.era, columns=c("TXNAME")), era.all.summ.ranges.ran, type="any")
  overlap.genes.mel.ran <- subsetByOverlaps( GenomicFeatures::genes(txdb.mel, columns=c("TXNAME")), mel.all.summ.ranges.ran, type="any")
  
  # store names, around 5500 genes overlap with these, whereas in real is 4698 || mel random overlap 3191 and real overlap 2500
  genes.era.ran <-  as.data.frame(overlap.genes.era.ran$TXNAME)$value; length(genes.era.ran); length(unique(genes.era.ran))
  genes.mel.ran <-  as.data.frame(overlap.genes.mel.ran$TXNAME)$value; length(genes.mel.ran); length(unique(genes.mel.ran))
  
  ###### store gene set  #####
  # grab all genes that overlap from outliers from above
  era.genes.ran <-  tibble(genes=genes.era.ran)
  mel.genes.ran <-  tibble(genes=genes.mel.ran)
  # unique, they were already unique but just in case
  era.genesOfInterest.ran.list <- as.character(unique(era.genes.ran$genes)); length(era.genesOfInterest.ran.list)
  mel.genesOfInterest.ran.list <- as.character(unique(mel.genes.ran$genes)); length(mel.genesOfInterest.ran.list)
  
  ###### create gene set with gene universe (those with goterms)  #####
  # indicate in the unverse whether gene present in list or not 
  # create named gene list from universe with 0s and 1s
  era.geneList.ran.list <- factor(as.integer(era.geneUniverse %in% era.genesOfInterest.ran.list))
  names(era.geneList.ran.list) <- era.geneUniverse
  
  mel.geneList.ran.list <- factor(as.integer(mel.geneUniverse %in% mel.genesOfInterest.ran.list))
  names(mel.geneList.ran.list) <- mel.geneUniverse
  
  ###### era - create GOdataset, run enrichment, store #####
  era.myGOdata.ran.list <- new("topGOdata", description="outliers", ontology="BP", allGenes=era.geneList.ran.list, annot = annFUN.gene2GO, gene2GO = era.geneID2GO)
  result.weight.Fisher.era.ran.list <- runTest(era.myGOdata.ran.list, algorithm="weight01", statistic="fisher",node=10 )
  era.allRes.ran.list <- GenTable(era.myGOdata.ran.list, weightFisher=result.weight.Fisher.era.ran.list , orderBy = "weightFisher"  , ranksOf ="result.weight.Fisher.era.ran.list", topNodes = 50)
  era.allRes.ran.list$Term <- gsub(" [a-z]*\\.\\.\\.$", "", era.allRes.ran.list$Term); era.allRes.ran.list$Term <- gsub("\\.\\.\\.$", "", era.allRes.ran.list$Term); era.allRes.ran.list$Term <- paste( era.allRes.ran.list$Term, era.allRes.ran.list$GO.ID,sep=", ")
  era.allRes.ran.list$Term <- factor(era.allRes.ran.list$Term, levels=rev(era.allRes.ran.list$Term)); 
  era.allRes.ran.list$weightFisher <- as.numeric(era.allRes.ran.list$weightFisher); era.allRes.ran.list <- subset(era.allRes.ran.list, weightFisher<0.05)
  

  ###### mel- create GOdataset, run enrichment, store #####
  mel.myGOdata.ran.list <- new("topGOdata", description="outliers", ontology="BP", allGenes=mel.geneList.ran.list, annot = annFUN.gene2GO, gene2GO = mel.geneID2GO)
  result.weight.Fisher.mel.ran.list <- runTest(mel.myGOdata.ran.list, algorithm="weight01", statistic="fisher",node=10 )
  mel.allRes.ran.list <- GenTable(mel.myGOdata.ran.list, weightFisher=result.weight.Fisher.mel.ran.list , orderBy = "weightFisher"  , ranksOf ="result.weight.Fisher.mel.ran.list", topNodes = 50)
  mel.allRes.ran.list$Term <- gsub(" [a-z]*\\.\\.\\.$", "", mel.allRes.ran.list$Term); mel.allRes.ran.list$Term <- gsub("\\.\\.\\.$", "", mel.allRes.ran.list$Term); mel.allRes.ran.list$Term <- paste( mel.allRes.ran.list$Term, mel.allRes.ran.list$GO.ID,sep=", ")
  mel.allRes.ran.list$Term <- factor(mel.allRes.ran.list$Term, levels=rev(mel.allRes.ran.list$Term)); 
  mel.allRes.ran.list$weightFisher <- as.numeric(mel.allRes.ran.list$weightFisher); mel.allRes.ran.list <- subset(mel.allRes.ran.list, weightFisher<0.05)
  
  ###### compare sharing, store  #####
  # store number of outlier terms, number shared between species per perm
  era.mel.go.terms.sharing.ran[i,] <- c(i, nrow(era.allRes.ran.list), nrow(mel.allRes.ran.list),sum(era.allRes.ran.list$Term %in% mel.allRes.ran.list$Term ),
                                        length(era.geneList.ran.list[era.geneList.ran.list==1]), 
                                        mean(subset(era.gene.lengths, tx_name %in% names(era.geneList.ran.list[era.geneList.ran.list==1|era.geneList.ran.list==0]))$width),
                                        mean(subset(era.gene.lengths, tx_name %in% names(era.geneList.ran.list[era.geneList.ran.list==1]))$width),
                                        length(mel.geneList.ran.list[mel.geneList.ran.list==1]), 
                                        mean(subset(mel.gene.lengths, tx_name %in% names(mel.geneList.ran.list[mel.geneList.ran.list==1|mel.geneList.ran.list==0]))$width),
                                        mean(subset(mel.gene.lengths, tx_name %in% names(mel.geneList.ran.list[mel.geneList.ran.list==1]))$width),
                                        paste0(subset(era.allRes.ran.list, Term %in% mel.allRes.ran.list$Term )$Term, collapse = '|')  )
  }

head(era.mel.go.terms.sharing)
# what percetange of the outlier GO terms are shared across species
era.mel.go.terms.sharing$perc.outlier.go.shared.era <- (as.numeric(as.character(era.mel.go.terms.sharing$no.terms.shared.between.sp))
                                                        /as.numeric(as.character(era.mel.go.terms.sharing$no.outlier.terms.era)))*100
era.mel.go.terms.sharing$perc.outlier.go.shared.mel <- (as.numeric(as.character(era.mel.go.terms.sharing$no.terms.shared.between.sp))
                                                        /as.numeric(as.character(era.mel.go.terms.sharing$no.outlier.terms.mel)))*100

era.mel.go.terms.sharing.ran$perc.outlier.go.shared.era <- (as.numeric(as.character(era.mel.go.terms.sharing.ran$no.terms.shared.between.sp))/
                                                              as.numeric(as.character(era.mel.go.terms.sharing.ran$no.outlier.terms.era)))*100
era.mel.go.terms.sharing.ran$perc.outlier.go.shared.mel <- (as.numeric(as.character(era.mel.go.terms.sharing.ran$no.terms.shared.between.sp))/
                                                              as.numeric(as.character(era.mel.go.terms.sharing.ran$no.outlier.terms.mel)))*100

######## collate results and plot ########
str(era.mel.go.terms.sharing); str(era.mel.go.terms.sharing.ran)
era.mel.go.terms.sharing.obs.ran <- rbind(era.mel.go.terms.sharing, era.mel.go.terms.sharing.ran)
era.mel.go.terms.sharing.obs.ran$perm.type <- if_else(era.mel.go.terms.sharing.obs.ran$perm.no=="observed","observed", "random.outliers" )

#write.csv(era.mel.go.terms.sharing.obs.ran, "1.data/gene.set.anno/era.mel.go.terms.sharing.obs.ran.1000.perm.csv", row.names = F)
#write.csv(era.mel.go.terms.sharing.obs.ran, "1.data/gene.set.anno/era.mel.go.terms.sharing.obs.ran2.1000.perm.csv", row.names = F)

era.mel.go.terms.sharing.obs.ran <- read.csv( "1.data/gene.set.anno/era.mel.go.terms.sharing.obs.ran.1000.perm.csv")

library(ggbeeswarm)
era.perm.plot <- ggplot(data=era.mel.go.terms.sharing.obs.ran, aes(x=perm.type, y=perc.outlier.go.shared.era))+
  geom_boxplot(outlier.colour = NA)+geom_beeswarm(alpha=.2)+ylab("Percentage enriched GO-terms shared across species") +
  scale_x_discrete(labels=c("observed" = "Observed GWAS\noutliers", "random.outliers" = "Simulated\noutliers"))+
  ylim(0,52) + theme_classic()+theme(
    axis.title.x=element_blank(),axis.title.y=element_text(face = "bold", size = 12),
    axis.text.x = element_text(size=12,face = "bold"),axis.text.y = element_text(size=10)
); era.perm.plot 

mel.perm.plot <- ggplot(data=era.mel.go.terms.sharing.obs.ran, aes(x=perm.type, y=perc.outlier.go.shared.mel))+
  geom_boxplot(outlier.colour = NA)+geom_beeswarm(alpha=.2)+ylab("Percentage enriched GO-terms shared across species") +
  scale_x_discrete(labels=c("observed" = "Observed GWAS\noutliers", "random.outliers" = "Simulated\noutliers"))+
  ylim(0,52) + theme_classic()+theme(
    axis.title.x=element_blank(),axis.title.y=element_text(face = "bold", size = 12),
    axis.text.x = element_text(size=12,face = "bold"),axis.text.y = element_text(size=10)
  ); mel.perm.plot 

plot_grid(era.perm.plot, mel.perm.plot, labels = c("H. erato", "H. melpomene"), label_fontface = "italic", label_x=0.1)

# gene length, add mean for all genes
era.perm.plot <- ggplot(data=era.mel.go.terms.sharing.obs.ran, aes(x=perm.type, y=mean.era.outlier.gene.length))+geom_boxplot()+geom_point()+
  annotate("text", x = 1.5, y = era.mel.go.terms.sharing.obs.ran$mean.era.gene.length, label = ">> mean gene length")+
  ylab("Mean gene length overlapping with outliers") +
  scale_x_discrete(labels=c("observed" = "Observed GWAS\noutliers", "random.outliers" = "Simulated\noutliers"))+
  ylim(15000,35000) + theme_classic()+theme(
    axis.title.x=element_blank(),axis.title.y=element_text(face = "bold", size = 12),
    axis.text.x = element_text(size=12,face = "bold"),axis.text.y = element_text(size=10)
  ); era.perm.plot 

mel.perm.plot <- ggplot(data=era.mel.go.terms.sharing.obs.ran, aes(x=perm.type, y=mean.mel.outlier.gene.length))+geom_boxplot()+geom_point()+
  annotate("text", x = 1.5, y = era.mel.go.terms.sharing.obs.ran$mean.mel.gene.length, label = ">> mean gene length")+
  ylab("Mean gene length overlapping with outliers") +
  scale_x_discrete(labels=c("observed" = "Observed GWAS\noutliers", "random.outliers" = "Simulated\noutliers"))+
  ylim(5000,18000) + theme_classic()+theme(
    axis.title.x=element_blank(),axis.title.y=element_text(face = "bold", size = 12),
    axis.text.x = element_text(size=12,face = "bold"),axis.text.y = element_text(size=10)
  ); mel.perm.plot 

plot_grid(era.perm.plot, mel.perm.plot, labels = c("H. erato", "H. melpomene"), label_fontface = "italic", label_x=0.2)

plot_grid(era.perm.plot, mel.perm.plot, labels = c("era", "mel"))

era.perm.plot <- ggplot(data=era.mel.go.terms.sharing.obs.ran, aes(x=perm.type, y=no.era.outlier.genes))+geom_boxplot()+geom_point(); era.perm.plot 
mel.perm.plot <-  ggplot(data=era.mel.go.terms.sharing.obs.ran, aes(x=perm.type, y=no.mel.outlier.genes))+geom_boxplot()+geom_point(); mel.perm.plot
plot_grid(era.perm.plot, mel.perm.plot, labels = c("era", "mel"))

era.perm.plot <- ggplot(data=era.mel.go.terms.sharing.obs.ran, aes(x=perm.type, y=no.outlier.terms.era))+geom_boxplot()+geom_point(); era.perm.plot 
mel.perm.plot <-  ggplot(data=era.mel.go.terms.sharing.obs.ran, aes(x=perm.type, y=no.outlier.terms.mel))+geom_boxplot()+geom_point(); mel.perm.plot
plot_grid(era.perm.plot, mel.perm.plot, labels = c("era", "mel"))


era.perm.plot <- ggplot(data=era.mel.go.terms.sharing.obs.ran, aes(x=perc.outlier.go.shared.era, y=no.outlier.terms.era, colour=perm.type))+
  geom_smooth(method="lm")+geom_point()+theme(legend.position = "none"); era.perm.plot 
mel.perm.plot <- ggplot(data=era.mel.go.terms.sharing.obs.ran, aes(x=perc.outlier.go.shared.mel, y=no.outlier.terms.mel, colour=perm.type))+
  geom_smooth(method="lm")+geom_point(); mel.perm.plot
plot_grid(era.perm.plot, mel.perm.plot, labels = c("era", "mel"), rel_widths = c(0.7,1))

## what go-terms tend to get enriched on both species
# from text to column r
go.terms.shared.df <- separate(era.mel.go.terms.sharing.obs.ran, go.terms.shared.era.mel,  paste(rep("go.term.", 11), 1:11, sep=""), sep = '\\|', remove = F)
head(go.terms.shared.df )
go.terms.shared.df.long <- gather(go.terms.shared.df[,c(1,12:22)], "go.term.no", "go.terms" , -perm.no)


go.terms.shared.df.long.summ <-  summarise(group_by(subset(go.terms.shared.df.long, perm.no!="observed"&go.terms!=""), go.terms),
          n.times.shared.sims=n())

go.terms.shared.df.long.summ$shared.in.observed <- if_else( go.terms.shared.df.long.summ$go.terms %in% 
                                                             subset(go.terms.shared.df.long, perm.no=="observed"&go.terms!="")$go.terms , "yes", "no")

# percentage of times these categories were enriched AND shared across species
go.terms.shared.df.long.summ$perc.sims.shared.across.sp <- (go.terms.shared.df.long.summ$n.times.shared.sims/1000)*100

#### save results #####
# all data with go terms
#write.csv(go.terms.shared.df, "1.data/gene.set.anno/era.mel.go.terms.shared.df.obs.ran.1000.perm.csv", row.names = F)

# summary of go-terms enriched more often and shared across species
write.csv(go.terms.shared.df.long.summ, "1.data/gene.set.anno/era.mel.go.terms.shared.df.long.summ.obs.ran.1000.perm.csv", row.names = F)
era.allRes$perc.sims.shared.across.sp <- go.terms.shared.df.long.summ$perc.sims.shared.across.sp[match(
  era.allRes$Term, go.terms.shared.df.long.summ$go.terms)]
mel.allRes$perc.sims.shared.across.sp <- go.terms.shared.df.long.summ$perc.sims.shared.across.sp[match(
  mel.allRes$Term, go.terms.shared.df.long.summ$go.terms)]

# split go.terms
era.allRes$Term.only <- sapply(str_split(era.allRes$Term, ","), `[[`, 1)
mel.allRes$Term.only <- sapply(str_split(mel.allRes$Term, ","), `[[`, 1)

## add percentage times categories were enriched AND shared across species in sumlations
head(era.allRes); head(mel.allRes)
nrow(era.allRes); nrow(mel.allRes)
write.csv(mel.allRes, "1.data/gene.set.anno/mel.allRes.short.fisher.csv", row.names = F)
write.csv(era.allRes, "1.data/gene.set.anno/era.allRes.short.fisher.csv", row.names = F)

