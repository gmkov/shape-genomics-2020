##### GMK April 2021 #####
#### ld decay plotting ####

rm(list=ls())
dev.off()
######################### packages ####


library(devtools)
library("windowscanr")
#install_github("drveera/ggman")
library(ggplot2)
library(dplyr)
library(tidyr)
library(qqman)
library(cowplot)
library(data.table)
library(stringr)
library(devtools)
library(ggman)
library(cowplot)
library(EnvStats)
library(zoo)
library(ggpubr)
options(scipen = 999)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(RIdeogram)
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/PhD/29.shape.paper/data/data.clean/")

######################### 0. prep data ######
era.ld <-read.delim("1.data/ld.decay/era.vcf.thinned.22318kb.minr2.0.05.minuschr02.ld.summary",sep="",header=T,check.names=F,stringsAsFactors=F)
mel.ld <-read.delim("1.data/ld.decay/mel.vcf.thinned.17205kb.minr2.0.2.summary",sep="",header=T,check.names=F,stringsAsFactors=F)

colnames(era.ld) <- c("dist","rsq")
colnames(mel.ld) <- c("dist","rsq")


era.ld$distc <- cut(era.ld$dist,breaks=seq(from=min(era.ld$dist)-1,to=max(era.ld$dist)+1,by=100000))
era.ld1 <- era.ld %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
era.ld1 <- era.ld1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                       end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                       mid=start+((end-start)/2))

era.p <- ggplot(data=era.ld1,aes(x=start,y=mean))+
  geom_point(size=0.4,colour="grey20")+
  stat_smooth(method = "gam", colour="black")+
  #geom_line(data=era.ld1,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
  labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+  theme_bw()+
  geom_vline(xintercept = 1200, lty=2)+
  scale_x_continuous(breaks=c(0,1*10^6,2*10^6,3*10^6),labels=c("0","1","2","3"), limits =c(0,3000000) );era.p


mel.ld$distc <- cut(mel.ld$dist,breaks=seq(from=min(mel.ld$dist)-1,to=max(mel.ld$dist)+1,by=100000))
mel.ld1 <- mel.ld %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
mel.ld1 <- mel.ld1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                              end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                              mid=start+((end-start)/2))

mel.p <- ggplot(data=mel.ld1,aes(x=start,y=mean))+
  geom_point(size=0.4,colour="grey20")+
  stat_smooth(method = "gam", colour="black")+
  #geom_line(data=mel.ld1,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
  labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+  theme_bw()+
  geom_vline(xintercept = 1032, lty=2)+
  scale_x_continuous(breaks=c(0,1*10^6,2*10^6,3*10^6),labels=c("0","1","2","3"), limits =c(0,3000000) );mel.p



plot_grid(NA,era.p, mel.p, rel_heights = c(0.05, 1,1),ncol = 1, labels = c("", "A. H. erato", "B. H. melpomene"), 
          label_x = 0, label_y = 1.05, align = "v", hjust = 0)



