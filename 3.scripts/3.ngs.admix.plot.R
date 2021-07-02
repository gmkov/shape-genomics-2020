#### size batch cam coll oct 18 ####

rm(list=ls())
dev.off()
library(ggplot2)
library(dplyr)
library(tidyr)
library(qqman)
library(cowplot)
library(data.table)
library(stringr)
library(tools)
library(ggpubr)
library(EnvStats)
#install.packages('stringr')
#### load data ####
haplo <- read.csv("~/Dropbox/PhD/29.shape.paper/data/joana/haplotagging.samples.all.csv",stringsAsFactors = FALSE)

# order according to population
head(haplo)
ord<-order(haplo[haplo$species=="erato"&haplo$aspect.ratio!="",][,"TransectPosition"])
ord<-order(haplo[haplo$species=="erato"&haplo$aspect.ratio!="",][,"CollectionSite"])

barplot(t(haplo[haplo$species=="erato"&haplo$aspect.ratio!="",]$ngs.admix.joana.k2)[,ord],col=1:3,space=0,border=NA,xlab="Individuals",ylab="Demo2 Admixture proportions for K=3")
text(tapply(1:nrow(haplo[haplo$species=="erato"&haplo$aspect.ratio!="",]),haplo[haplo$species=="erato"&haplo$aspect.ratio!="",][ord,1],mean),-0.05,unique(haplo[haplo$species=="erato"&haplo$aspect.ratio!="",][ord,1]),xpd=T)
abline(v=cumsum(sapply(unique(haplo[haplo$species=="erato"&haplo$aspect.ratio!="",][ord,1]),function(x){sum(haplo[haplo$species=="erato"&haplo$aspect.ratio!="",][ord,1]==x)})),col=1,lwd=1.2)


# order by transect position
haplo <- read.csv("~/Dropbox/git/shape-genomics-2020/1.data/2.haplotagging.samples.all.csv",stringsAsFactors = FALSE)
names(haplo)
haplo.era.order <- haplo[haplo$species=="erato"&haplo$aspect.ratio!=""&haplo$TransectPosition!="",]
haplo.era.order$TransectPosition <- as.numeric(as.character(haplo.era.order$TransectPosition))
haplo.era.order <-haplo.era.order[order(haplo.era.order$TransectPosition),]
haplo.era.order$CollectionSite <- factor(haplo.era.order$CollectionSite , 
                                         levels = unique(haplo.era.order$CollectionSite[order(haplo.era.order$TransectPosition)]))
haplo.era.order$id <- factor(haplo.era.order$id,  levels = unique(haplo.era.order$id[order(haplo.era.order$TransectPosition)]))

era.ngs.summ <- dplyr::summarise(dplyr::group_by(haplo.era.order, subsp, CollectionSite),
                                 n=n(), 
                                 admix=mean(ngs.admix.k2),
                                 TransectPosition=as.numeric(as.character(unique(TransectPosition))))

ggplot(era.ngs.summ  , aes(fill=subsp, y=admix, x=CollectionSite)) + 
  geom_bar(position="dodge", stat="identity")+theme(axis.text.x = element_text(angle=45))

era.admix.plot <-ggplot(subset(era.ngs.summ,n>2)  , aes(colour=subsp, y=admix, x=TransectPosition)) + 
  geom_point(aes(size=n))+theme(axis.text.x = element_text(angle=45))+theme_classic()+
  xlim(0,85)+
  guides(colour=guide_legend(title="Subspecies"), size=guide_legend(title="Sample size"))+
  ylab("Admixture proportion")+ xlab("Position along transect (km)")+
  scale_colour_manual(values=c("black", "#009E73", "#E69F00"), limits=c("notabilis", "hybrid", "lativitta"))+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size=12, face="bold", colour = "black"),
        legend.position = c("right"), axis.text = element_text(size=12, face="bold", colour = "black"),
        #axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.line = element_line( colour = "black"),
        legend.background = element_rect(fill = "transparent")); era.admix.plot

## mel 
haplo.mel.order <- haplo[haplo$species=="melpomene"&haplo$aspect.ratio!=""&haplo$TransectPosition!="",]
mel.ngs.summ <- dplyr::summarise(dplyr::group_by(haplo.mel.order, subsp, CollectionSite),
                                 n=n(), 
                                 admix=mean(ngs.admix.joana.k2),
                                 TransectPosition=as.numeric(as.character(unique(TransectPosition))))


ggplot(mel.ngs.summ  , aes(fill=subsp, y=admix, x=CollectionSite)) + 
  geom_bar(position="stack", stat="identity")+theme(axis.text.x = element_text(angle=45))

mel.ngs.summ
mel.admix.plot <- ggplot(subset(mel.ngs.summ,n>2)  , aes(colour=subsp, y=1-admix, x=TransectPosition)) + 
  geom_point(aes(size=n))+theme(axis.text.x = element_text(angle=45))+theme_classic()+
  ylab("Admixture proportion")+ xlab("Position along transect (km)")+
  guides(colour=guide_legend(title="Subspecies"), size=guide_legend(title="Sample size"))+
scale_colour_manual(values=c("black", "#009E73", "#E69F00"), limits=c("plesseni", "hybrid", "malleti"))+
  xlim(0,85)+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size=12, face="bold", colour = "black"),
        legend.position = c("right"), axis.text = element_text(size=12, face="bold", colour = "black"),
        #axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.line = element_line( colour = "black"),
        legend.background = element_rect(fill = "transparent")) ; mel.admix.plot

plot_grid(era.admix.plot, mel.admix.plot, ncol = 1, labels = "AUTO")
ggsave("~/Dropbox/PhD/29.shape.paper/figs/supplementary/admixture.prop.png", width = 7, height = 7)

