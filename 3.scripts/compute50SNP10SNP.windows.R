##!Rscript --vanilla

# Load packages required
require("data.table")
require("devtools")
require("windowscanr")

# Read in the vcf file name
args<-commandArgs(TRUE)
file<-args[1]
ass.dat <- fread(file, header=T, sep = '\t', stringsAsFactors=F)

# modify headers
# grab chromosome no., for erato 7,8 ; for mel 6,7
ass.dat$CHR <- as.integer(substr(ass.dat$Chromosome,7,8))

# Rename LRTscore to LRT if present (only for CC)
names(ass.dat)[names(ass.dat)=="LRTscore"]<-"LRT"

# remove mitochondrial genome
ass.dat <- subset(ass.dat, ass.dat$CHR!="_m")

### replace -999 (failed) with NA ###
ass.dat[ass.dat==-999] <- NA
ass.dat <- subset(ass.dat, LRT!="")

# convert LRT to pvalues (to be able to average)
ass.dat$LRT.pval <- pchisq(ass.dat$LRT,1, lower.tail = F)

#### Generate window averages #####

# windows, position fixed
ass.dat.windows.pos <- winScan(x = ass.dat,  
                               position = NULL,
                               values = c("LRT.pval", "Position"), 
                               groups = "Chromosome",
                               win_size = 50,
                               win_step = 10,
                               cores=20,
                               funs = c("min","max","median"))                   

# Write out results of window analyses:
write.table(ass.dat.windows.pos,file=paste0(file,".50snp.10snp.pos"), sep = '\t',quote=FALSE)

