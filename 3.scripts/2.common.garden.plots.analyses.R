######## packages ########
rm(list=ls())
dev.off()
#devtools::install_github("kassambara/ggpubr")
library(ggplot2)
library(dplyr)
library(tidyr)
library(qqman)
library(cowplot)
library(data.table)
library(stringr)
library(tools)
library(ggpubr)
library(rptR)
library(EnvStats)
library(MASS)
library(lme4)
library(lmerTest)
library(relaimpo)
library(lmSupport)
library(ggbeeswarm)
library(EnvStats) # for adding n to boxplots
library(plyr)
library(survminer)

######## data ########
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/git/shape-genomics-2020/")
df <- read.csv("1.data/1.master.era.mel.wild.reared.csv", na.strings = c("", " ", NA))
df <- subset(df, sex!="")
df$altitude.type<- if_else(df$altitude>600,"hig", "low")
df.reared <-  subset(df, area!=""& type.reared=="reared"&!(is.na(mother.id))&mother.id!="")

### Fig. S1 aspect ratio is normally distributed ###
gghistogram(subset(df, !(sex==""&area!="")), x = "aspect.ratio",add = "mean", rug = TRUE,add.params=list(size=1.5, linetype="longdash"),
            color = "transparent", fill="sex",alpha=.5, bins=20, palette = c("black", "orange"), facet.by = c("type.reared","species"), scales="free_y",
            panel.labs.background = list(fill = "lightgrey", color = "lightgrey"), panel.labs.font.x=list(face="bold.italic",size=12)
            , panel.labs.font.y=list(face="bold",size=12), ylab = "Count", xlab="Wing Aspect Ratio",
            panel.labs = list(type.reared = c("Reared", "Wild"), species = c("H. erato", "H. melpomene")))

############################## 1. plots. mother-offspring regression / altitude boxplot ###############
##### 1.0 data prep ######
# mothers with wing phenotypes
mothers.with.pheno<- subset(df, is.mother.brood=="yes"&area!=""&!(is.na(area)))
mothers<- subset(df, is.mother.brood=="yes")

### era melpomene main fig 2
summ <- dplyr::summarise(dplyr::group_by(subset(df, area!=""& (type.reared=="reared")&!(is.na(mother.id))), mother.id,species,rearing.batch),
                         n=n(),
                         aspect.ratio.mean=mean(aspect.ratio, na.rm=TRUE),
                         aspect.ratio.sd=sd(aspect.ratio),
                         aspect.ratio.n=length(aspect.ratio),
                         aspect.ratio.se=sd(aspect.ratio)/sqrt(aspect.ratio.n),
                         area.mean=mean(area, na.rm=TRUE),
                         area.sd=sd(area),
                         area.n=length(area),
                         area.se=sd(area)/sqrt(area.n),
                         area=mean(area),
                         altitude=mean(altitude),
                         dev.time=mean(dev.time,na.rm=T),
                         no.fem=length(subset(sex, sex=="female")),
                         no.male=length(subset(sex, sex=="male")),
                         ratio=no.male/no.fem,
                         pupa.clean.mass=mean(pupa.clean.mass))
summ$ar.mother <-  mothers$aspect.ratio[match(summ$mother.id,mothers$mother.id)]
summ$area.mother <-  mothers$area[match(summ$mother.id,mothers$mother.id)]
summ$altitude.type <-if_else(summ$altitude>600,"hig","low")
summ$locality <- df$locality[match(summ$mother.id, df$mother.id)]
levels(summ$species) <- c("H. erato", "H. melpomene")

# subset so that all broods have at least 3 phenotyped offspring
df.more3 <- subset(df, area!=""& aspect.ratio!=""& mother.id!=""&!(is.na(mother.id))& type.reared=="reared" &mother.id %in% summ[summ$n>2,]$mother.id)
summ$species <- factor(summ$species, labels = c("H. erato", "H. melpomene"))
df.more3$species <- factor(df.more3$species, labels = c("H. erato", "H. melpomene"))
df.more3$mother.id <-  droplevels(df.more3$mother.id)
summ.more3<- subset(summ, n>2)
summ.more3$mother.unit.id <- mothers$cam.id[match(summ.more3$mother.id, mothers$mother.id)]

# check trait variances
sd(mothers$aspect.ratio); sd(summ$aspect.ratio.mean); sd(summ.sex[summ.sex$sex=="male",]$aspect.ratio.mean); sd(summ.sex[summ.sex$sex=="female",]$aspect.ratio.mean)
sd(df.more3[df.more3$sex=="male"&df.more3$species=="H. erato"&df.more3$type.reared=="reared",]$aspect.ratio); sd(df.more3[df.more3$sex=="female"&df.more3$type.reared=="reared"&df.more3$species=="H. erato",]$aspect.ratio)

mean(df.more3[df.more3$sex=="male",]$aspect.ratio); mean(df.more3[df.more3$sex=="female",]$aspect.ratio)
sd(df.more3[df.more3$sex=="male"&df.more3$species=="H. melpomene",]$aspect.ratio); sd(df.more3[df.more3$sex=="female"&df.more3$species=="H. melpomene",]$aspect.ratio)

# how many broods with mothers data and more than 2
nrow(subset(summ, n>2&ar.mother!="")); nrow(subset(summ, n>2))
nrow(subset(summ, species=="H. erato"& n>2 &ar.mother!=""));nrow(subset(summ, species=="H. erato"& n>2 ))
nrow(subset(summ, species=="H. melpomene"& n>2 &ar.mother!=""));nrow(subset(summ, species=="H. melpomene"& n>2))

# how many era mel
sum(subset(summ, n>2&species=="H. erato")$n)
sum(subset(summ, n>2&species=="H. melpomene")$n)

##### 1.1 fig.2 A C - AR mother-offspring regression #####

fig2.a.c <-  ggplot(subset(summ, n>2), aes(x= ar.mother, y=aspect.ratio.mean))+
  geom_abline(slope=1, intercept=, colour="grey",lty=2)+
  stat_smooth(method="lm", colour="black")+
  stat_regline_equation(aes(),label.y.npc = 0.85, label.x = 1.96,font.label = c("bold"), size=5)+
  #geom_vline(xintercept = mean(subset(summ, n>2&species=="H. erato"&!(is.na(ar.mother)))$ar.mother), colour="red")+
  #geom_hline(yintercept = mean(subset(summ, n>2&species=="H. erato"&!(is.na(ar.mother)))$aspect.ratio.mean), colour="red")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), digits = 1, label.x = 1.96, size=5)+
  #stat_cor( size=5)+
  geom_errorbar(width=0, aes(ymin=aspect.ratio.mean-aspect.ratio.se, ymax=aspect.ratio.mean+aspect.ratio.se), alpha=.9, colour="darkgrey") +
  geom_point(aes(size=n))+
  ylab("Mid-offspring aspect ratio")+xlab("Mother aspect ratio")+ labs(size = "Brood size")+
  facet_wrap(~species, ncol = 1, scales = "free")+
  scale_y_continuous(breaks=seq(1.8,2.3,0.1), limits = c(1.96,2.225))+
  scale_x_continuous(breaks=seq(1.8,2.3,0.1), limits = c(1.96,2.225))+
  theme_bw()+ theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
                    panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
                    panel.grid.major = element_blank(), axis.text=element_text(size=12),
                    panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold",vjust = 1),
                    legend.text.align = 0, legend.position = "none",strip.text = element_text(size=16,face="bold.italic", hjust = 0),
                    strip.background =element_rect(fill="transparent", color="transparent")
  ); fig2.a.c

##### 1.2 fig.S5 A - area mother-offspring regression #####
fig2.area.a <-  ggplot(subset(summ, n>2), aes(x=area.mother, y=area.mean))+
  geom_abline(slope=1, intercept=, colour="grey",lty=2)+
  stat_smooth(data=subset(subset(summ,n>2), species=="H. melpomene" ),method="lm", colour="black")+
  stat_regline_equation(font.label = c("bold"), size=5, label.y.npc = 0.9)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), digits = 1, label.x = 340, label.y = 550, size=5)+
  geom_errorbar(width=0, aes(ymin=area.mean-area.se, ymax=area.mean+area.se), alpha=.9, colour="darkgrey") +
  geom_point(aes(size=n))+
  ylab("Mid-offspring wing area")+xlab("Mother wing area")+ labs(size = "Brood size")+
  facet_wrap(~species, ncol = 1, scales = "free")+
  scale_y_continuous(breaks=seq(300,550,50), limits = c(340,550))+
  scale_x_continuous(breaks=seq(300,550,50), limits = c(340,550))+
  theme_bw()+ theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
                    panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
                    panel.grid.major = element_blank(), axis.text=element_text(size=12),
                    panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold",vjust = 1),
                    legend.text.align = 0, legend.position = "none",strip.text = element_text(size=16,face="bold.italic", hjust = 0),
                    strip.background =element_rect(fill="transparent", color="transparent")); fig2.area.a

##### 1.3 fig.S2 B - AR centered/scaled mother-offspring regressions #####
## could scale and centre phenotypes, to account for different variances between parents/offsrping

supp.ar.centred.era <-  ggplot(subset(summ, n>2&species=="H. erato"), aes(x= scale(ar.mother, center = T), y=scale(aspect.ratio.mean, center = T)))+
  geom_abline(slope=1, intercept=, colour="grey",lty=2)+
  stat_smooth(method="lm", colour="black")+
  stat_regline_equation(aes(),label.y.npc = 0.85,font.label = c("bold"), size=5)+
  #geom_vline(xintercept = mean(subset(summ, n>2&species=="H. erato"&!(is.na(ar.mother)))$ar.mother), colour="red")+
  #geom_hline(yintercept = mean(subset(summ, n>2&species=="H. erato"&!(is.na(ar.mother)))$aspect.ratio.mean), colour="red")+
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), digits = 1, label.x = 1.96, size=5)+
  stat_cor( size=5)+
  #geom_errorbar(width=0, aes(ymin=scale(aspect.ratio.mean, center = T)-scale(aspect.ratio.se, center = T), 
  #                           ymax=scale(aspect.ratio.mean, center = T)+scale(aspect.ratio.se, center = T)), alpha=.9, colour="darkgrey") +
  geom_point(aes(size=n))+
  ylab("Mid-offspring aspect ratio\n(centred, scaled)")+xlab("Mother aspect ratio (centred, scaled)")+ labs(size = "Brood size")+
  #facet_wrap(~species, ncol = 1, scales = "free")+
  #scale_y_continuous(breaks=seq(1.8,2.3,0.1), limits = c(1.96,2.225))+
  #scale_x_continuous(breaks=seq(1.8,2.3,0.1), limits = c(1.96,2.225))+
  theme_bw()+ theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
                    panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
                    panel.grid.major = element_blank(), axis.text=element_text(size=12),
                    panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold",vjust = 1),
                    legend.text.align = 0, legend.position = "none",strip.text = element_text(size=16,face="bold.italic", hjust = 0),
                    strip.background =element_rect(fill="transparent", color="transparent")
  ); supp.ar.centred.era 

supp.ar.centred.mel <-  ggplot(subset(summ, n>2&species=="H. melpomene"), aes(x= scale(ar.mother, center = T), y=scale(aspect.ratio.mean, center = T)))+
  geom_abline(slope=1, intercept=, colour="grey",lty=2)+
  stat_smooth(method="lm", colour="black")+
  stat_regline_equation(aes(),label.y.npc = 0.85,font.label = c("bold"), size=5)+
  #geom_vline(xintercept = mean(subset(summ, n>2&species=="H. erato"&!(is.na(ar.mother)))$ar.mother), colour="red")+
  #geom_hline(yintercept = mean(subset(summ, n>2&species=="H. erato"&!(is.na(ar.mother)))$aspect.ratio.mean), colour="red")+
  #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), digits = 1, label.x = 1.96, size=5)+
  stat_cor( size=5)+
  #geom_errorbar(width=0, aes(ymin=scale(aspect.ratio.mean, center = T)-scale(aspect.ratio.se, center = T), 
  #                           ymax=scale(aspect.ratio.mean, center = T)+scale(aspect.ratio.se, center = T)), alpha=.9, colour="darkgrey") +
  geom_point(aes(size=n))+
  ylab("")+xlab("Mother aspect ratio (centred, scaled)")+ labs(size = "Brood size")+
  #facet_wrap(~species, ncol = 1, scales = "free")+
  #scale_y_continuous(breaks=seq(1.8,2.3,0.1), limits = c(1.96,2.225))+
  #scale_x_continuous(breaks=seq(1.8,2.3,0.1), limits = c(1.96,2.225))+
  theme_bw()+ theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
                    panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
                    panel.grid.major = element_blank(), axis.text=element_text(size=12),
                    panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold",vjust = 1),
                    legend.text.align = 0, legend.position = "none",strip.text = element_text(size=16,face="bold.italic", hjust = 0),
                    strip.background =element_rect(fill="transparent", color="transparent")
  ); supp.ar.centred.mel

##### 1.4 fig.2 B D - AR high/low F1 boxplot #####

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.075, 1), symbols = c("****", "***", "**", "*", "+", "ns"))

dplyr::summarise(dplyr::group_by(subset(summ,n>2), species, altitude.type),ar=mean(aspect.ratio.mean))
t.test(aspect.ratio.mean~altitude.type, subset(summ,n>2& species=="H. erato"))
t.test(aspect.ratio.mean~altitude.type, subset(summ,n>2& species=="H. melpomene"))

my_comparisons1 <- list( c("hig", "low") )
fig2.b <-ggplot(subset(summ,n>2), aes(x=altitude.type, y=aspect.ratio.mean)) +
  geom_beeswarm(data=subset(df.more3,type.reared=="reared"), aes(x=altitude.type, y=aspect.ratio), alpha=.9, 
                inherit.aes = F, colour="lightgrey" , size=1, dodge.width = 1)+
  stat_n_text(y.pos = 1.9)+
  geom_boxplot(alpha=.5) + #geom_beeswarm(alpha=.6, aes(size=n), groupOnX = T,dodge.width = 1)+
  geom_jitter(aes(size=n), width = .05, alpha=.5)+
  ylab("Mid-offspring aspect ratio")+ labs(size = "Brood size")+
  scale_y_continuous(breaks=seq(1.8,2.3,0.1))+
  scale_x_discrete(limits=c("hig", "low"), labels = c("F1\nHigh", "F1\nLow"))+facet_wrap(~species, scales = "free", ncol = 1)+
  stat_compare_means(comparisons = my_comparisons1,method = "t.test", label = "p.signif",symnum.args=symnum.args, bracket.size = 0.25, 
                     step.increase=.01, tip.length = 0, vjust = 0.4, size=6 )+
  theme_classic()+theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text.x=element_text(size=14,face="bold"),axis.text.y=element_text(size=12),
        panel.grid.minor = element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "right",strip.text = element_blank(),
        strip.background =element_rect(fill="transparent", color="transparent")
  ); fig2.b

##### 1.5 fig.S5 B - area high/low boxplot #####
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.074, 1), symbols = c("****", "***", "**", "*", "+", "ns"))
dplyr::summarise(dplyr::group_by(subset(summ,n>2), species, altitude.type),ar=mean(aspect.ratio.mean))
t.test(area.mean~altitude.type, subset(summ,n>2& species=="H. erato"))
t.test(area.mean~altitude.type, subset(summ,n>2& species=="H. melpomene"))

my_comparisons1 <- list( c("hig", "low") )
fig2.area.c  <-ggplot(subset(summ,n>2), aes(x=altitude.type, y=area.mean)) +
  geom_beeswarm(data=subset(df.more3,type.reared=="reared"), aes(x=altitude.type, y=area), alpha=.9, 
                inherit.aes = F, colour="lightgrey" , size=1, dodge.width = 1)+
  stat_n_text(y.pos = 230)+
  geom_boxplot(alpha=.5) + #geom_beeswarm(alpha=.6, aes(size=n), groupOnX = T,dodge.width = 1)+
  geom_jitter(aes(size=n), width = .05, alpha=.5)+
  ylab("Mid-offspring aspect ratio")+ labs(size = "Brood size")+
  #scale_y_continuous(breaks=seq(1.8,2.3,0.1))+
  scale_x_discrete(limits=c("hig", "low"), labels = c("F1\nHigh", "F1\nLow"))+facet_wrap(~species, scales = "free", ncol = 1)+
  stat_compare_means(comparisons = my_comparisons1,method = "t.test", label = "p.signif",symnum.args=symnum.args, bracket.size = 0.25,
                     step.increase=.05, tip.length = 0, vjust = 0.2, size=6)+
  theme_classic()+theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text.x=element_text(size=14,face="bold"),axis.text.y=element_text(size=12),
        panel.grid.minor = element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank(),
        legend.text.align = 0, legend.position = "right",strip.text = element_blank(),
        strip.background =element_rect(fill="transparent", color="transparent")
  ); fig2.area.c

##### 1.6 fig.S2 A split by sex ######
summ.sex <- dplyr::summarise(group_by(subset(df, area!=""& type.reared=="reared"&!(is.na(mother.id))&mother.id!=""), mother.id,species,rearing.batch,sex),
                  n=n(),
                  aspect.ratio.mean=mean(aspect.ratio, na.rm=TRUE),
                  aspect.ratio.sd=sd(aspect.ratio),
                  aspect.ratio.n=length(aspect.ratio),
                  aspect.ratio.se=sd(aspect.ratio)/sqrt(aspect.ratio.n),
                  area=mean(area),
                  altitude=mean(altitude),
                  dev.time=mean(dev.time,na.rm=T),
                  no.fem=length(subset(sex, sex=="female")),
                  no.male=length(subset(sex, sex=="male")),
                  ratio=no.male/no.fem)
summ.sex$ar.mother <-  mothers$aspect.ratio[match(summ.sex$mother.id,mothers$mother.id)]
summ.sex$area.mother<-  mothers$area[match(summ.sex$mother.id,mothers$mother.id)]
summ.sex$altitude.type <-if_else(summ.sex$altitude>600,"high","low")
levels(summ.sex$species) <- c("H. erato", "H. melpomene")

fig.supp.2 <-  ggplot(subset(summ.sex, n>2), aes(x=ar.mother, y=aspect.ratio.mean, colour=sex, group=sex, fill=sex))+
  geom_abline(slope=1, intercept=, colour="grey",lty=2)+
  stat_smooth(method="lm")+
  stat_regline_equation(aes(),label.y.npc = 0.8,font.label = c("bold"), size=5)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 1.95, size=5, p.accuracy = 0.001, r.accuracy = 0.01)+
  geom_errorbar(width=0, aes(ymin=aspect.ratio.mean-aspect.ratio.se, ymax=aspect.ratio.mean+aspect.ratio.se), alpha=.7) +
  geom_point(aes(size=n))+
  #geom_point(size=3)+
  ylab("Mid-offspring aspect ratio")+xlab("Mother aspect ratio")+ labs(size = "Brood size")+
  facet_wrap(~species)+
  #geom_text(aes(label=mother.id), nudge_x = .01)+
  scale_color_manual(values = c("black", "orange"))+scale_fill_manual(values = c("black", "orange"))+
  #scale_colour_viridis_d(end = .7)+scale_fill_viridis_d(end = .7)+
  ylim(1.94,2.225)+xlim(1.94,2.225)+theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text=element_text(size=12),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold"),
        legend.text.align = 0, legend.position = "right",strip.text = element_text(size=14,face="bold.italic"),
        strip.background =element_rect(fill="transparent", color="transparent")
  ); fig.supp.2 

# fig. S2
ggarrange(fig.supp.2, ggarrange(supp.ar.centred.era, supp.ar.centred.mel, NULL, ncol = 3, labels = c("B", "",""), widths = c(1,1,.2)), # Second row with box and dot plots
          nrow = 2, labels = "A",align = "hv" ) 

##### 1.7 fig.S5 B area split by sex ######
fig2.area.b <-  ggplot(subset(summ.sex, n>2), aes(x=area.mother, y=area, colour=sex, group=sex, fill=sex))+
  geom_abline(slope=1, intercept=, colour="grey",lty=2)+
  stat_smooth(data=subset(subset(summ.sex, n>2), species=="H. melpomene" ), method="lm")+
  stat_regline_equation(aes(),label.y.npc = 0.7,font.label = c("bold"), size=5,label.x = 300)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 300, size=5, p.accuracy = 0.001, r.accuracy = 0.01)+
  #geom_errorbar(width=0, aes(ymin=area-area.se, ymax=area+area.se), alpha=.7) +
  geom_point(aes(size=n))+
  #geom_point(size=3)+
  ylab("Mid-offspring wing area")+xlab("Mother wing area")+ labs(size = "Brood size")+
  facet_wrap(~species, ncol = 1)+
  #geom_text(aes(label=mother.id), nudge_x = .01)+
  scale_color_manual(values = c("black", "orange"))+scale_fill_manual(values = c("black", "orange"))+
  #scale_colour_viridis_d(end = .7)+scale_fill_viridis_d(end = .7)+
  #ylim(1.94,2.225)+xlim(1.94,2.225)+theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text=element_text(size=12),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold"),
        legend.text.align = 0, legend.position = "none",#strip.text = element_text(size=14,face="bold.italic"),
        strip.text = element_blank(),
        strip.background =element_rect(fill="transparent", color="transparent"),axis.title.y=element_blank()
  ); fig2.area.b

##### 1.8 fig.S4 AR wild/reared high/low t.tests ####
# use same altitude range as mothers
hist(summ$altitude, breaks = 20) # very few mothers between 650 and 1000- so remove those from the wild sample
min(mothers$altitude); max(mothers$altitude)
# prep dataset, use min/max altitude mothers collected
df.wild <- subset(df, altitude>374&altitude<1522&(altitude>1000|altitude<650)&
                    sex!=""&type.reared=="wild"&is.na(mother.id))

df.wild.cg <-data_frame(indiv.brood=df.wild$cam.id,altitude=df.wild$altitude,altitude.type=df.wild$altitude.type, aspect.ratio=df.wild$aspect.ratio, sex=df.wild$sex, type.reared=df.wild$type.reared, species=df.wild$species); head(df.wild.cg)
levels(df$species) <- c("H. erato", "H. melpomene")
df.cg <-data_frame(indiv.brood=subset(summ.sex, n>2)$mother.id,altitude=subset(summ.sex, n>2)$altitude,altitude.type=subset(summ.sex, n>2)$altitude.type, aspect.ratio=subset(summ.sex, n>2)$aspect.ratio.mean, 
                   sex=subset(summ.sex, n>2)$sex, type.reared=rep("reared", nrow(subset(summ.sex, n>2))), species=subset(summ.sex, n>2)$species); head(df.cg)
df.wild.cg <-rbind(df.wild.cg, df.cg)
df.wild.cg$altitude.type <- substr(df.wild.cg$altitude.type, 0,3)

# prep t.tests
library(rstatix)
library(dplyr)
stat.test <-  df.wild.cg %>% dplyr::group_by(altitude.type, species,sex) %>% t_test(aspect.ratio~type.reared) %>% add_significance()
stat.test <-  stat.test %>% add_xy_position(x="altitude.type")

test <- df.wild.cg %>% dplyr::group_by(altitude.type, species,sex) 
# plot 
reared.wild.sex.sp.boxplot <- ggplot(aes(x=altitude.type, y=aspect.ratio, colour=type.reared), data = subset(df.wild.cg, aspect.ratio!=""))+
  geom_boxplot(alpha=.7,  aes(fill=type.reared))+ 
  geom_point(aes(group=type.reared, fill=type.reared), position = position_dodge(width = 0.75))+
  stat_pvalue_manual(stat.test,label="{p.signif}")+
  ylab("Aspect ratio")+ xlab("Altitude")+labs(colour="",fill="")+
  facet_wrap(~species+sex)+
  scale_x_discrete(labels=c("Highland", "Lowland"))+
  scale_color_manual(values = c("black", "grey"))+scale_fill_manual(values = c("black", "grey"))+
  theme_bw()+ theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
                    panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
                    panel.grid.major = element_blank(), axis.text=element_text(size=12),
                    panel.grid.minor = element_blank(), axis.title=element_text(size=12,face="bold",vjust = 1),
                    legend.text.align = 0, legend.position = "right",strip.text = element_text(size=12,face="bold", hjust = 0, vjust = -.5),
                    strip.background =element_rect(fill="transparent", color="transparent"));reared.wild.sex.sp.boxplot 

############################## 2. analyses common garden broods #######
############### 2.1 erato lmm #####
## checks, prep ##
levels(df$species) <- c("erato", "melpomene")

## scale explanatory variables, only
era.df.more3 <- subset(subset(df, aspect.ratio!=""& mother.id!=""&species=="erato"&type.reared=="reared"), mother.id %in% summ[summ$n>2,]$mother.id)

nrow(era.df.more3); length(unique(era.df.more3$mother.id))
hist(era.df.more3$altitude)
era.df.more3$altitude.scaled <- scale(era.df.more3$altitude ); hist(era.df.more3$altitude.scaled)
era.df.more3$aspect.ratio.scaled <- scale(era.df.more3$aspect.ratio ); hist(era.df.more3$aspect.ratio.scaled)
era.df.more3$area.scaled <- scale(era.df.more3$area ); hist(era.df.more3$area.scaled)
era.df.more3$dev.time.scaled <- scale(era.df.more3$dev.time ); hist(era.df.more3$dev.time)

# colinearity
library(olsrr)
m1.lm <- lm(aspect.ratio ~ altitude.type +area.scaled+sex+dev.time.scaled , data=subset(era.df.more3, species=="erato") ) ; summary(m1.lm)
ols_vif_tol(m1.lm)
ols_eigen_cindex(m1.lm)

## 1. lmerTest erato model selection ##

m1.lmer <- lmerTest::lmer(aspect.ratio ~ altitude.type +area.scaled+sex+dev.time.scaled  + (1|rearing.batch) + (1|altitude.type:rearing.batch) + (1|altitude.type:rearing.batch:mother.id), 
                          data=subset(era.df.more3, species=="erato") ) ; plot(m1.lmer);summary(m1.lmer)

# equivalent , because for every mother.id there is only ever ONE altitude.type and ONE rearing batch
m1.lmer1 <- lmerTest::lmer(aspect.ratio ~ altitude.type +area.scaled+sex+dev.time.scaled  + (1|rearing.batch/mother.id) + (1|mother.id), 
                           data=subset(era.df.more3, species=="erato") ) ; plot(m1.lmer1);summary(m1.lmer1)

# 2. for model random selection use REML fit
lmerTest::ranova(m1.lmer1) # drop nesting
m1.lmer1 <- lmerTest::lmer(aspect.ratio ~ altitude.type +area.scaled+sex+dev.time.scaled  + (1|rearing.batch) + (1|mother.id), 
                           data=subset(era.df.more3, species=="erato") ) ; plot(m1.lmer1);summary(m1.lmer1)
lmerTest::ranova(m1.lmer1) # drop rearing batch
m2.lmer <- lmerTest::lmer(aspect.ratio ~ altitude.type +area.scaled+sex+dev.time.scaled + (1|mother.id), 
                          data=subset(era.df.more3, species=="erato") ) ; plot(m2.lmer);summary(m2.lmer)

# 3. for fixed effect selection use ML fit
drop1(update(m2.lmer, REML = F), test = "Chisq") # drop nothing

m3.lmer <- lmerTest::lmer(aspect.ratio ~ altitude.type +area.scaled+sex+dev.time.scaled + (1|mother.id), 
                          data=subset(era.df.more3, species=="erato") ) ; plot(m3.lmer);summary(m3.lmer)

# checks
library(ggResidpanel)
resid_panel(m3.lmer)

# results to report
summary(m3.lmer)
lmerTest::ranova(m3.lmer) 

## R2 mumin ##
# marginal:  proportion of variance explained by the fixed factor(s)
# conditional: proportion of variance explained by both the fixed and random factors. info: https://jonlefcheck.net/2013/03/13/r2-for-linear-mixed-effects-models/comment-page-1/
MuMIn::r.squaredGLMM((m3.lmer))

# similar values when estimated with other packages
performance::r2((m3.lmer))
piecewiseSEM::rsquared((m3.lmer))

############### 2.2 melpomene lmm #################
## checks, prep ##
## scale explanatory variables, only
mel.df.more3 <- subset(subset(df, aspect.ratio!=""& mother.id!=""&species=="melpomene"&dev.time!=""), mother.id %in% summ[summ$n>2,]$mother.id)
# development time does not get selected so we can drop it and use more datapoints
mel.df.more3 <- subset(subset(df, aspect.ratio!=""& mother.id!=""&species=="melpomene"&(type.reared=="reared")), mother.id %in% summ[summ$n>2,]$mother.id)
nrow(mel.df.more3); length(unique(mel.df.more3$mother.id))
mel.df.more3$altitude.scaled <- scale(mel.df.more3$altitude, center = T ); hist(mel.df.more3$altitude.scaled)
mel.df.more3$aspect.ratio.scaled <- scale(mel.df.more3$aspect.ratio, center = T  ); hist(mel.df.more3$aspect.ratio.scaled)
mel.df.more3$area.scaled <- scale(mel.df.more3$area, center = T  ); hist(mel.df.more3$area.scaled)
mel.df.more3$dev.time.scaled <- scale(mel.df.more3$dev.time, center = T  ); hist(mel.df.more3$dev.time)


## 1. lmerTest melpomene model selection ##
m1.lmer <- lmerTest::lmer(aspect.ratio ~ altitude.type +area.scaled+sex+dev.time.scaled  + (1|mother.id), 
                          data=subset(mel.df.more3, species=="melpomene") ) ; plot(m1.lmer);summary(m1.lmer)

# 2. for model random selection use REML fit
lmerTest::ranova(m1.lmer) # dont drop
m2.lmer <- lmerTest::lmer(aspect.ratio ~ altitude.type +area.scaled+sex+dev.time.scaled + (1|mother.id), 
                          data=subset(mel.df.more3, species=="melpomene") ) ; plot(m2.lmer);summary(m2.lmer)

# 3. for fixed effect selection use ML fit
drop1(update(m2.lmer, REML = F), test = "Chisq") # drop dev time

m3.lmer <- lmerTest::lmer(aspect.ratio ~ altitude.type +area.scaled+sex + (1|mother.id), 
                          data=subset(mel.df.more3, species=="melpomene") ) ; plot(m3.lmer);summary(m3.lmer)
# checks
library(ggResidpanel)
resid_panel(m3.lmer)

# results to report
summary(m3.lmer)
lmerTest::ranova(m3.lmer) 

## R2 mumin ##
# marginal:  proportion of variance explained by the fixed factor(s)
# conditional: proportion of variance explained by both the fixed and random factors. info: https://jonlefcheck.net/2013/03/13/r2-for-linear-mixed-effects-models/comment-page-1/
MuMIn::r.squaredGLMM((m3.lmer))

# similar values when estimated with other packages
performance::r2((m3.lmer))
piecewiseSEM::rsquared((m3.lmer))



############### 2.3 anova / repeatability   #####
# size
size.aov<-aov(area ~ mother.id, data = subset(df.more3, species=="H. erato") ); summary(size.aov)
rep.size<- rpt(area  ~ area+sex+(1|mother.id), grname=c("mother.id"), datatype="Gaussian", 
                nboot = 100, data=subset(era.df.more3, species=="erato") ); print(rep.size) ; plot(rep.size, cex.main = 1)

size.aov<-aov(area ~ mother.id, data = subset(df.more3, species=="H. melpomene") ); summary(size.aov)
rep.size<- rpt(area  ~ area+sex+(1|mother.id), grname=c("mother.id"), datatype="Gaussian", 
                nboot = 100, data=subset(df.more3, mother.id!=""&species=="H. melpomene"&type.reared=="reared") ); print(rep.size) ; plot(rep.size, cex.main = 1)

# shape
shape.aov<-aov(aspect.ratio ~ mother.id, data = subset(df.more3, species=="H. erato") ); summary(shape.aov)
rep.shape<- rpt(aspect.ratio  ~ area.scaled +altitude.type+ dev.time.scaled+ sex+(1|mother.id), grname=c("mother.id"), datatype="Gaussian", 
                nboot = 100, data=subset(era.df.more3, species=="erato") ); print(rep.shape) ; plot(rep.shape, cex.main = 1)

shape.aov<-aov(aspect.ratio ~ mother.id, data = subset(df.more3, species=="H. melpomene") ); summary(shape.aov)
rep.shape<- rpt(aspect.ratio  ~ area.scaled +altitude.type+  sex+(1|mother.id), grname=c("mother.id"), datatype="Gaussian", 
                nboot = 100, data=subset(mel.df.more3, species=="melpomene")  ); print(rep.shape) ; plot(rep.shape, cex.main = 1)


############### ###############
