#### packages ####
rm(list=ls())
dev.off()

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(devtools)
library(BEDASSLE)
library(ggplot2)
library(vcfR)
library(fossil)
library(viridis)
library(patchwork)
library(sp)
library(raster)
library(ecodist)
library(radiator)
library(wesanderson)
library(gridExtra)
#devtools::install_github("dkahle/ggmap")
library(ggmap)
register_google(key="AIzaSyAz8Urlhyb4VqYuiy_dBIv-ietj7eY4YBo", write = TRUE)
options(scipen = 999)
library(RgoogleMaps)
library(topoDistance)
library(elevatr)
library(googleway)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
require(mapdata)
library(ggspatial)
library(ggsn)
library(cowplot)

# must register for google APIs (free for education) 
# https://developers.google.com/maps/documentation/javascript/get-api-key
# https://lucidmanager.org/data-science/geocoding-with-ggmap/
# enable other apis https://stackoverflow.com/questions/60061173/error-in-ggmap-must-be-an-array-and-http-400-bad-request

##### functions ##### 
#### 0. prep data ####
setwd("/Users/gabrielamontejokovacevich/Dropbox (Cambridge University)/git/shape-genomics-2020/")
coll <- read.csv("1.data/4.mothers.haplotagging.formap.csv")
coll <- subset(coll, !(is.na(latitude)))

#### 1. Fig. 1A sampling map #####
mean(coll$latitude);mean(coll$longitude)
coll.map <- get_googlemap(center=c(lon=-77.8, lat=-1.24), zoom = 9, color="bw", maptype = "terrain",
                          style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")

coll.x <- c(-78.4, -77.4); coll.y <- c(-1.8, -.5)
gc(); coll.map.gg <- ggmap(coll.map) +
  geom_jitter(data = subset(coll, experiment=="haplotagging"),
  aes(x = longitude, y = latitude, shape = experiment, fill=pattern, group=pattern),
            size=2, alpha=.8, width = .02, height = .02, stroke=.5)+
  scale_shape_manual(values = c(24,21))+
  scale_fill_manual(values=c("black", "#009E73", "#E69F00"))+
  scale_x_continuous(limits = c(-78.4, -77.4, expand = c(0, 0))) +
  scale_y_continuous(limits = c(-1.8, -.5), expand = c(0, 0)) +
  # add common garden mother triangles
  geom_jitter(data = subset(coll, experiment=="cg"),
  aes(x = longitude, y = latitude, shape = experiment, fill=pattern, group=pattern),
  size=4, alpha=.85, width = .02, height = .02, stroke=.5, fill="#E69F00")+
  # add cross for ikiam
  geom_point(data = subset(coll, substr(locality,0,3)=="Iki"),
              aes(x = longitude, y = latitude, shape = experiment, fill=pattern, group=pattern),
              size=4, alpha=.8, width = .02, height = .02, stroke=3, fill="black", shape=4)+
  # this will add a scalebar of 50km
  scalebar(dist = 25, dist_unit = "km",location="bottomleft",
           st.bottom = TRUE, st.color = "black", height=0.01, st.dist = 0.03,
           x.min = coll.x[1], x.max = coll.x[2], size=1,box.fill=c("black", "black"),
           y.min = coll.y[1]+0.07, y.max = coll.y[2], 
           transform = TRUE, model = "WGS84") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title = element_blank()) +
  guides(size = guide_legend(order=2)); coll.map.gg

# ggsave("fig.1.a.png", scale = 1, dpi = 300)


## Fig. 1A south america insert with box matching ## 
extend(coll.map.gg); GetMap.bbox(coll.map.gg)
map <- get_googlemap(center=c(lon=-74.54, lat=-0.2602158), zoom = 5, color="bw",  maptype = "terrain",
                          style = "feature:road|visibility:off&style=element:labels|visibility:off&style=visibility|feature:off&style=administrative.country:on")
gc(); p1 <- ggmap(map) +
  scale_x_continuous(limits = c(-88, -65), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-9, 13), expand = c(0, 0)) +
  labs(size = "Elevation (m)") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none", 
        axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title = element_blank()) +
  geom_rect(mapping = aes(xmin=coll.x[1], xmax=coll.x[2], ymin=coll.y[1], ymax=coll.y[2]),fill=NA,color="black", alpha=.5)+
  guides(size = guide_legend(order=2)); p1

#ggsave("fig.1.a.SA.png")

#### 2. Fig. 1C elevation topographic cline #####
# prep data, mean lat/lon per population, and pop names
# make bins according to transect length, these will be used for getting the path (no need for pattern differentiation)
coll$location.bin <- cut(coll$transect.position,10, labels=c(1:10))
summ.haplo.bin <- dplyr::summarise(dplyr::group_by(subset(coll, experiment=="haplotagging"&transect.position!=""), location.bin),
                            n=dplyr::n(),
                            longitude=mean(longitude),
                            latitude=mean(latitude),
                            altitude=mean(altitude),
                            transect.position=mean(transect.position)); summ.haplo.bin
# settings for elevatr
prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# get points into right format
summ.haplo.bin.sp <- SpatialPoints(summ.haplo.bin[,3:4], proj4string = CRS(prj_dd))
# obtain dem(digital elevation model) 
summ.haplo.bin.dem <-  get_elev_raster(summ.haplo.bin.sp, z = 9, prj = prj_dd)

# LEFT: plot elevation map, check points along transect
plot(summ.haplo.bin.dem )
points(summ.haplo.bin[,3:4])

# use summ.haplo.bin.dem obtained with elavatr for the map in the preamble
# get topographic distance, paths=TRUE to then plot paths
# order by location bin, only estimate path from point 1 to 10 of transect
setorder(summ.haplo.bin, location.bin)
tdist <- topoDist(summ.haplo.bin.dem$layer, summ.haplo.bin[,3:4][c(1,10),], paths = TRUE); tdist

# plot map with path, checks
topoPathMap(summ.haplo.bin.dem$layer,summ.haplo.bin[,3:4], topoPaths = tdist, type = "hillshade",
            pathWidth = 4, cex = 2, bg = "blue")
# plot topo surface, can choose type = "plotly", checks
topoProfile(summ.haplo.bin.dem$layer, topoPaths = tdist, pts = 1000, 
            type = "base", singlePlot = TRUE )

### extract data from tdist for ggplot plotting
topoPaths <- tdist[[2]]
pathLength <- SpatialLinesLengths(topoPaths, longlat = TRUE)
samplePts <- spsample(topoPaths[1], n = 1000, type = "regular")
pathCoordinates <- samplePts@coords
elevations<- extract(summ.haplo.bin.dem$layer, samplePts)
pathDists <- seq(0, pathLength[1], by = pathLength[1]/(1000 -  1))

#make dataframe with topo data
path.df <- data_frame(elevations=elevations, pathDists=pathDists, longitude=pathCoordinates[,1], latitude=pathCoordinates[,2])

# find closest longitude in path to longitude in sites of summary haplo, but with colour pattern info
summ.haplo.bin.pattern <- summarise(group_by(subset(coll, experiment=="haplotagging"&transect.position!=""), location.bin,pattern),
                                    n=n(),
                                    longitude=mean(longitude),
                                    latitude=mean(latitude),
                                    altitude=mean(altitude),
                                    transect.position=mean(transect.position)); summ.haplo.bin.pattern
# closest longitude
closest<-function(xv,sv){xv[which(abs(xv-sv)==min(abs(xv-sv)))] }
closest.lon <- list()
for (i in 1:nrow(summ.haplo.bin.pattern)){
  closest.lon[i] <- closest(path.df$longitude, summ.haplo.bin.pattern$longitude[i])
}; head(closest.lon)
summ.haplo.bin.pattern$longitude.match.path <- unlist(closest.lon)
# add  closest path position for plotting to each binned datapoint
summ.haplo.bin.pattern$path.dist <- path.df$pathDists[match(summ.haplo.bin.pattern$longitude.match.path, path.df$longitude)]


ggplot(aes(x=path.dist, y=altitude,  size=n, colour=pattern, fill=pattern), data=subset(summ.haplo.bin.pattern, transect.position!=""))+
  annotate("rect", xmin=-4,xmax=0.00, ymin=0, ymax=1600, colour='transparent',fill="tan4", alpha=0.4)+
  geom_area(inherit.aes = F, aes(x=pathDists, y=elevations), data=path.df , colour='transparent',fill="tan4", alpha=0.4)+
  coord_cartesian(ylim=c(330,1450), xlim=c(-4,79))+
  #ylim(150,1450)+xlim(-15,75)+
  geom_jitter(shape=21, colour="black", width = 2.1, height = 0.1,alpha=.8)+
  #geom_point(inherit.aes = F, aes(x=transect.position, y=altitude), data=summ.haplo.bin[c(1,9),], size=5, colour='tan4',fill="tan4")+
  theme_classic()+
  scale_size(range = c(4, 20))+ labs(size = "Sample Size")+
  #scale_shape_manual(values = c(22))+
  scale_x_continuous(expand = c(0, 0))+scale_y_continuous(expand = c(0, 0))+ 
  #coord_fixed(10/100)+
  scale_fill_manual(values=c("black", "#009E73", "#E69F00"))+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = c("right"), axis.text = element_text(size=12, face="bold", colour = "black"),
        #axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.line = element_line( colour = "black"),
        legend.background = element_rect(fill = "transparent")) # get rid of legend bg

#ggsave("fig.1.c.hybrid.cline.png",  bg = "transparent", scale=1, width = 8, height = 4)



#### NOT IN USE 3. barlpot aspect ratio #### 
library(ggpubr)
setwd("~/Dropbox/PhD/29.shape.paper/")
df <- read.csv("~/Dropbox/PhD/29.shape.paper/data/wild.previous/era.mel.races.wild.haplo.csv", na.strings = c("", " ", NA))
head(df)
df$experiment <- coll$experiment[match(df$cam.id,coll$cam.id)]

gghistogram(subset(df, !(sex==""&area.mm2!="")&unit.type=="wild"&
                     (subsp=="lativitta"|subsp=="hybrid"|subsp=="notabilis")&
                     (group=="haplotagging"|group=="ps.wild")), x = "aspect.ratio",add = "mean", rug = TRUE,add.params=list(size=1.5, linetype="solid"),
            color = "transparent", fill="subsp",alpha=.5, bins=10, facet.by = c("sex"), scales="free_y")

gghistogram(subset(df, !(sex==""&area.mm2!="")&unit.type=="wild"&
                     (subsp=="lativitta"|subsp=="hybrid"|subsp=="notabilis")&
                     (group=="haplotagging")), x = "aspect.ratio",add = "mean", rug = TRUE,add.params=list(size=1.5, linetype="solid"),
            color = "transparent", fill="subsp",alpha=.5, bins=10, facet.by = c("sex"), scales="free_y")+theme_bw()+
  #scale_x_discrete(limits=c("notabilis", "hybrid", "lativitta"))+
  scale_fill_manual(values=c("forestgreen", "darkorange","black"))+
  scale_colour_manual(values=c("forestgreen", "darkorange","black"))+
  ylab("Aspect ratio")+xlab("")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text.x=element_text(size=12),axis.text.y=element_text(size=14),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold"),axis.title.x = element_text(), axis.title.y = element_text(size=14),
        legend.text.align = 0, legend.position = "none",strip.text = element_text(size=14,face="bold"),
        strip.background =element_rect(fill="transparent", color="transparent")  )

gghistogram(subset(df, !(sex==""&area.mm2!="")&unit.type=="wild"&
                     (subsp=="malleti"|subsp=="hybrid"|subsp=="plesseni")&
                     (group=="haplotagging"|group=="ps.wild")), x = "aspect.ratio",add = "mean", rug = TRUE,add.params=list(size=1.5, linetype="solid"),
            color = "transparent", fill="subsp",alpha=.5, bins=10, facet.by = c("sex"), scales="free_y")


levels(summ$species) <- c("H. erato", "H. melpomene")

my_comparisons1 <- list( c("notabilis", "lativitta"), c("notabilis", "hybrid"), c( "hybrid", "lativitta") )
era <- ggplot(aes(x=subsp, y=aspect.ratio, fill=subsp), data= subset(df, !(sex==""&area.mm2!="")&unit.type=="wild"&
                                                                       (subsp=="lativitta"|subsp=="hybrid"|subsp=="notabilis")&(group=="haplotagging")))+
  geom_violin()+
  stat_compare_means(comparisons = my_comparisons1, method = "t.test", label = "p.signif", bracket.size = 0.25, step.increase=.05, tip.length = 0, vjust = 0.2, size=4)+theme_bw()+
  scale_x_discrete(limits=c("notabilis", "hybrid", "lativitta"))+
  scale_fill_manual(values=c("forestgreen", "darkorange","black"))+
  scale_colour_manual(values=c("forestgreen", "darkorange","black"))+
  ylab("Aspect ratio")+xlab("")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text.x=element_text(size=12,face="bold.italic"),axis.text.y=element_text(size=14),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold"),axis.title.x = element_text(), axis.title.y = element_text(size=14),
        legend.text.align = 0, legend.position = "none",strip.text = element_text(size=14,face="bold.italic"),
        strip.background =element_rect(fill="transparent", color="transparent")  ); era 

my_comparisons2 <- list( c("plesseni", "malleti"), c("plesseni", "hybrid"), c("hybrid", "malleti") )
mel <- ggplot(aes(x=subsp, y=aspect.ratio, fill=subsp), data= subset(df, !(sex==""&area.mm2!="")&unit.type=="wild"&species=="melpomene"&
                                                                       (subsp=="plesseni"|subsp=="hybrid"|subsp=="malleti")&(group=="haplotagging"|group=="ps.wild")))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons2, method = "t.test", label = "p.signif", bracket.size = 0.25, step.increase=.05, tip.length = 0, vjust = 0.2, size=4)+theme_bw()+
  scale_x_discrete(limits=c("plesseni",  "hybrid","malleti"))+
  scale_fill_manual(values=c("forestgreen", "darkorange","black"))+
  scale_colour_manual(values=c("forestgreen", "darkorange","black"))+
  ylab("Aspect ratio")+xlab("")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.8),
        panel.background = element_blank(), legend.text = element_text(size=12), legend.title =  element_text(size=12),
        panel.grid.major = element_blank(), axis.text.x=element_text(size=12,face="bold.italic"),axis.text.y=element_text(size=14),
        panel.grid.minor = element_blank(), axis.title=element_text(size=14,face="bold"),axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.text.align = 0, legend.position = "none",strip.text = element_text(size=14,face="bold.italic"),
        strip.background =element_rect(fill="transparent", color="transparent")  ); mel

plot_grid(era, mel, labels = c("A", "B"), align = "hv")
ggsave("figs/supplementary/haplo.races.AR.png", width = 6, height = 4)
