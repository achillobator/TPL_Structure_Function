#Script 1: General script to import yeast cytometry data and graph time-course series
source("http://bioconductor.org/biocLite.R")
biocLite("flowCore")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("flowTime")

library(ggthemes)
library(ggplot2)
library(flowCore)
library(flowTime)

#Read in the data using flowtime and write the annotation CSV
read.plateSet <- function(path = getwd(), pattern = "\\."){ 
  files <- grep(pattern = pattern, x = list.files(path), value = TRUE)
  for(file in files) {
    plate <- read.flowSet(path = paste(path, file, sep = "/"), alter.names = TRUE)
    plate_num <- which(files == file) 
    sampleNames(plate) <- paste0(plate_num, sampleNames(plate))
    pData(plate)$name <- sampleNames(plate)
    pData(plate)$folder <- file
    if('flow_set' %in% ls()) flow_set <- rbind2(flow_set, plate) else flow_set <- plate
  }
  return(flow_set)
}
setwd("~/Documents/folder_with_multiple_fcs_files")
flow_set <- read.plateSet(path = "date_TC/", pattern = "date")

annotation <- read.csv("date_annotation.csv")
annotation$name <- gsub(x = annotation$name, pattern = "_", replacement = "")

flow_set <- annotateFlowSet(flow_set, annotation)

write.flowSet(flow_set, outdir = "date_name")

#write the annotation file
flow_set <- read.flowSet(path = "date_name", phenoData = "annotation.txt")

#dat_sum combines the annotation with the actual data
dat_sum <- summarizeFlow(flow_set, channel = "FL2.A", ploidy = "diploid", only = "yeast")

#nice data table :)
dat_sum

median_cl_boot <- function(x, conf = 0.95) {
  lconf <- (1 - conf)/2
  uconf <- 1 - lconf
  require(boot)
  bmedian <- function(x, ind) median(x[ind])
  bt <- boot(x, bmedian, 1000)
  bb <- boot.ci(bt, type = "perc")
  data.frame(y = median(x), ymin = quantile(bt$t, lconf), ymax = quantile(bt$t, uconf))
}

#make a boxplot of data to check data imported properly:
plot1<- ggplot(eventsfix, aes(as.factor(strain), y = FL2.Amean, color=as.factor(strain))) + 
  geom_boxplot(outlier.size = 0.1) + 
  stat_summary(fun.data = median_cl_boot, geom = "errorbar", width=.1) +
  stat_summary(fun.y=median, geom="point", shape=23, size=4) +
  theme_classic(base_family = 'Arial Bold', base_size = 10)  + 
  ylim(c(0,40000)) + 
  ylab('Venus Fluorescence (AU)') + 
  theme(legend.position="right")
plot1+theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Plot all to see changes in trends over time
plotall <- ggplot(data = eventsfix, mapping = aes(x = time, y = FL2.Amean, fill = as.factor(strain), color = TPL_type)) +
  geom_point(size=2) + 
  geom_line() + 
  ylab("Fluorescence (a.u.)") + 
  xlab("Time (min)") + 
  labs(color = "TPL type", subtitle = bquote("[Auxin]"), title = "Experiment") + 
  theme_classic(base_family = 'Arial Bold', base_size = 12) +
  geom_errorbar(aes(ymin=FL2.Amean-(FL2.Asd/sqrt(events)), ymax=FL2.Amean+(FL2.Asd/sqrt(events)))) + 
  scale_fill_discrete(guide = FALSE) + #this deletes the fill legend from the initial ggplot function (we only want the color legend)
  #scale_color_manual(values = A_colors) + #this allows us to manually put in a vector of colors that we want
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), strip.background  = element_blank()) #hjust = 0.5 centers the title and the subtitle instead of being left justified. strip.background = element_blank takes the box away from around the auxin concentration numbers
plotall
plotall+ facet_wrap(~"variable", ncol = 3)

#Colors- to use a manual scale color pallete
A_colors <- c( "#8512C0","#E43025","#E43025","#E43025","#E43025","#E43025","#E43025","#E43025","#E43025","#0A5526","#000000")

##############

#Script 2: representative code for plotting data as histograms in thise case for figure 4, where the inducible MED21 was utilized in SLARC strains.

#To generate histograms:
dat.SS_yeast_clip <- steadyState(flowset = flow_set, ploidy = 'diploid', only = 'yeast')
#subset data to only have TPL_type with TPL length N188 as an example
A_N188 <- c("N188")
N188_dat <- subset(dat.SS_yeast_clip, TPL_type %in% A_N188)

## Basic histogram from the vector "rating". Each bin is 500 wide.
hist1<-ggplot(N188_dat, aes(x=FSC.H)) + 
  geom_histogram(binwidth=500) + 
  theme_classic(base_family = 'Arial Bold', base_size = 15)
hist1 
hist1 + facet_wrap(~ Rapa, nrow =1)

# Density curve
hist2<-ggplot(N188_dat, aes(x=FSC.H)) + 
  geom_density() + 
  theme_classic(base_family = 'Arial Bold', base_size = 15)
hist2

#Histogram and density plots with multiple groups
hist3<-ggplot(N188_dat, aes(x=FL2.A, color=as.factor(iMED21))) + 
  geom_density() + 
  theme_classic(base_family = 'Arial Bold', base_size = 12) 
hist3+ scale_x_log10() 
hist3+ scale_x_log10() + facet_wrap(~ Rapa, nrow =1)
hist3+ scale_x_log10() + facet_grid(rows = vars(Rapa), cols = vars(Best))
hist3

hist3b<-ggplot(N188_dat, aes(x=FL2.A, color=as.factor(iMED21))) + geom_density() + theme_classic(base_family = 'Arial Bold', base_size = 12) 
hist3b+ scale_x_log10() 
hist3b+ scale_x_log10() + facet_wrap(~ Rapa, nrow =2) + theme(strip.background = element_blank(), strip.text.x = element_blank())
hist3b+ scale_x_log10() + facet_grid(rows = vars(Rapa), cols = vars(Best))
hist3b

A_colors <- c("#000000", "#33FF00", "#2570e8", "#FF0000", "#FFCC00","#FF9900","#CC9900","#0033CC","#0099FF","#006600","#000000")

hist3c<-ggplot(N188_dat, aes(x=FL2.A, color=as.factor(iMED21), fill=as.factor(iMED21))) + 
  geom_density(size=0.5, alpha=0.2) + 
  theme_classic(base_family = 'Arial Bold', base_size = 12) + 
  scale_color_manual(values = A_colors) + 
  scale_fill_manual(values = A_colors)
hist3c+ scale_x_log10() + theme(  strip.background = element_blank(), strip.text.x = element_blank())

# facet wrap by rapa treatment
hist3c+ scale_x_log10() + facet_wrap(~ Rapa, nrow =2) + theme(  strip.background = element_blank(), strip.text.x = element_blank()) 

hist5<-ggplot(dat.SS_yeast_clip, aes(x=FL2.A, color=as.factor(iMED21))) + 
  geom_density() + 
  theme_classic(base_family = 'Arial Bold', base_size = 8)
hist5+ scale_x_log10() + facet_grid(rows = vars(Rapa), cols = vars(Best))

A_colors <- c( "#000000", "#33FF00", "#2570e8", "#FFCC00","#FF9900","#CC9900","#0033CC","#0099FF","#006600","#000000")

#Alternative plots as bargraphs as per SG request
plot3c<- ggplot(N188_dat, aes(as.factor(AAtag), y = FL2.A, color=as.factor(iMED21))) + 
  geom_boxplot(outlier.size = 0.1) + 
  theme_classic(base_family = 'Arial Bold', base_size = 10)  + 
  facet_wrap(~ Rapa, nrow =1)+
  coord_cartesian(ylim = c(0, 60000)) +
  ylab('Venus Fluorescence (AU)') + 
  theme(legend.position="right")+
  scale_color_manual(values = A_colors)
plot3c+theme(axis.text.x = element_text(angle = 90, hjust = 1))

#######
#Script 3. Representative code for analyzing statistical significance between genotypes of yeast in cytometry, in thise case for figure 4, where the inducible MED21 was utilized in SLARC strains.


library(multcompView)
library(rcompanion)
library(lsmeans)
library(multcomp)

#Run ANNOVA
attach(dat)
data(dat)
str(dat)
tapply(iMED21, FL2.A, mean) 
tapply(iMED21, FL2.A, var)
tapply(iMED21, FL2.A, length)
boxplot(FL2.A ~ iMED21, data=dat)

lm.out = with(dat, lm(FL2.A ~ iMED21))
aov.out = aov(FL2.A ~ iMED21, data=dat)
oneway.test(FL2.A ~ iMED21, data=dat)
is.factor(FL2.A)
is.factor(iMED21)

aov.out
summary(aov.out)
TukeyHSD(aov.out)

summary.lm(aov.out)

model=lm( dat$FL2.A ~ dat$iMED21 )
ANOVA=aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'dat$genotype', conf.level=0.99)
plot(TUKEY , las=1 , col="brown")

marginal = lsmeans(model, ~ iMED21)
pairs(marginal, adjust="tukey")
CLD = cld(marginal, alpha = 0.01, Letters = letters, adjust  = "tukey")
CLD

# untreated Anova
model2=lm( dat$FL2.A ~ dat$iMED21 )
ANOVA2=aov(model2)
TUKEY2 <- TukeyHSD(x=ANOVA2, 'dat$genotype', conf.level=0.99)
plot(TUKEY2 , las=1 , col="brown")

marginal2 = lsmeans(model2, ~ iMED21)
pairs(marginal2, adjust="tukey")
CLD2 = cld(marginal2, alpha = 0.01, Letters = letters, adjust  = "tukey")
CLD2

########################
#Script 4: Plotting and statistical tests for Lateral root density:

library(ggthemes)
library(ggplot2)
library("viridis") 

dat4t1 <- read.csv("Supp_Datafile_LRD.csv")

B_line <- c("312","295","298","299","Col-0","slr","5057","252")
dat_B <- subset(dat8, genotype %in% B_line)
dat_B$genotype <- factor(dat_B$genotype, levels=c("Col-0","5057","298","299","252","295","312","slr"),ordered=TRUE)

plot1<- ggplot(dat4t1, aes(as.factor(genotype), y = LRD, color=as.factor(TPL_type), show.legend = FALSE)) + geom_boxplot(outlier.size = NULL) +
  geom_jitter(size =.2, width = 0.2) +
  #scale_y_continuous(limits=c(0,0.6)) +
  ylab('LRD') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), strip.background  = element_blank()) +
  theme_classic(base_family = 'Arial Bold', base_size = 10)+
  scale_color_viridis(discrete = TRUE, option = "D")
plot1+theme(axis.text.x = element_text(angle = 90, hjust = 1))+facet_wrap(~dpg, nrow=1)+ theme(legend.position="none")

plot2<- ggplot(dat_B, aes(as.factor(genotype), y = LRD, color =as.factor(genotype))) +   geom_boxplot(outlier.shape = NA, lwd=0.6) +
  geom_jitter(size=.5, width=.2, color="black", alpha=0.6) +
  ylab('LRD') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), strip.background  = element_blank()) +
  theme_classic(base_family = 'Arial Bold', base_size = 10)+
  scale_color_brewer(palette = "Paired")
plot2

#Run ANNOVA On dat_B

dat_B
attach(dat_B)
data(dat_B)
str(dat_B)
tapply(genotype, LRD, mean) 
tapply(genotype, LRD, var)
tapply(genotype, LRD, length)
boxplot(LRD ~ genotype, data=dat_B)

lm.out = with(dat_B, lm(LRD ~ genotype))
aov.out = aov(LRD ~ genotype, data=dat_B)
oneway.test(LRD ~ genotype, data=dat_B)
is.factor(LRD)
is.factor(genotype)

aov.out
summary(aov.out)
TukeyHSD(aov.out)
summary.lm(aov.out)

model=lm( dat_B$LRD ~ dat_B$genotype )
ANOVA=aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'dat_B$genotype', conf.level=0.95)
plot(TUKEY , las=1 , col="brown")
library(multcompView)
marginal = lsmeans(model, ~ genotype)
pairs(marginal, adjust="tukey")
CLD = cld(marginal, alpha = 0.05, Letters = letters, adjust  = "tukey")
CLD

###################
#Script 5. Sorting image files for measuring fluorescence from Anchor Away & SLARC plating experiments. Images have to be saved as text image, and the identical region of plate is saved in both the red and green channel with the format "strain_treatment_color.txt" (ie 3758_aux_green.txt where aux indicated the addition of auxin in the plate). Additionally an  annotation file is created to provide identities for the strains named "annotation.csv".
#Strain	Filter	aux	AAtag	TPL
#3758-1	green	0	None	iaa14
#3758-1	red	0	None	iaa14
#3758-1	green	10	None	iaa14
#3758-1	red	10	None	iaa14
#3759-1	green	0	None	H15
#3759-1	red	0	None	H15
#3759-1	green	10	None	H15
#3759-1	red	10	None	H15
#3760-1	green	0	None	n188
#3760-1	red	0	None	n188
#3760-1	green	10	None	n188
#3760-1	red	10	None	n188
#then

library(dplyr)
# combining image data

# set to local directory where files are located
directory <- "~/Documents/raw_image_data/"

# create list of files ending with "red.csv" or "green.csv"
files <- list.files(directory, pattern="(red|green)\\.txt$") 

# construnct the header rows of the file
properties <- sub("\\.txt$", "", files)
properties <- stringr::str_split(properties, "_", n=3, simplify=TRUE)
properties[, 2] <- sub("noaux", "0", properties[, 2]) # replace "no" with 0
properties[, 2] <- sub("aux", "10", properties[, 2]) 
properties <- t(properties[, c(1,3,2)])

# linearize and combine the data
df <- data.frame()
for (i in 1:length(files)) {
  file <- files[i]
  data <- read.table(file.path(directory, file), sep="\t")
  column <- c(t(data))
  df[1:length(column),i] <- column
}

# load annotations
annotationFile <- ("annotation.csv")
anno <- read.csv(annotationFile, sep=",", header=TRUE, stringsAsFactors=FALSE)

# # OLD FORMAT
# # attach the header to the top of the data
# combinedData <- rbind(properties, df, stringsAsFactors=FALSE)
# row.names(combinedData) <- as.numeric(row.names(combinedData)) - 3
# row.names(combinedData)[1:3] <- c()
# # save a tab separated file
# write.table(combinedData, file.path(directory, "combined_data.tsv"), col.names=FALSE, sep="\t")

# new format:
# combine properties and data
output <- cbind(as.data.frame(t(properties), stringsAsFactors=FALSE), 
                as.data.frame(t(df), stringsAsFactors=FALSE))
# melt data so each value has it's own row
output <- reshape2::melt(output, id.vars=c("V1","V2","V3"), value.name="FL", 
                         variable.name="pixel", na.rm=TRUE)
# re-order so rows with the same Strain, Filter, and Rapa are together
output <- output[with(output, order(V1, V2, V3)), ]
colnames(output)[1:3] <- c("Strain", "Filter", "aux")
# convert strain and aux to integers
output$aux <- as.integer(output$aux)
#output$Strain <- as.integer(output$Rapa)
# merge annotation to add cel, aa and slac columns
output2 <- left_join(output, anno)


write.table(output2, file.path(directory, "new_format_combined_data2.txt"), 
            col.names=TRUE, sep="\t", quote=FALSE, row.names=FALSE)

####################
#Script 6. Graphing plate data


library(ggthemes)
library(ggplot2)
library(flowCore)
library(flowTime)
library(Hmisc)
library(dplyr)
library(lsr)
setwd("~/Documents/wd")
#from above, load in the data file generated form the raw image files. 
Datafile_HG <- read.csv("new_format_combined_data2.txt", sep = "\t")
#the next step divides the green by the red values for the associated red and green pixel at the same location in the image
Datafile_HG = Datafile_HG[,c(2:8)]
Test = (longToWide(Datafile_HG, formula = FL~Filter))
Test$Norm_FL= (Test$FL_green/Test$FL_red)
test2<- Test
test2$TPL <- factor(test2$TPL, levels=c("IAA14","H1-5", "N188"),ordered=TRUE)
test2$AAtag <- factor(test2$AAtag, levels=c("None","med21", "del3med21","del5med21","del7med21"),ordered=TRUE)

plot1<- ggplot(test2, aes(as.factor(AAtag), y = Norm_FL, color=as.factor(aux))) + 
  geom_boxplot(outlier.size = 0.1) + 
  theme_classic(base_family = 'Arial Bold', base_size = 10)  + 
  ylab('Normalized Venus Fluorescence (8-bit)') + 
  theme(legend.position="right") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = B_colors)
plot1
plot1 + facet_wrap(~TPL, nrow=1)
plot1 + facet_grid(rows = vars(aux), cols = vars(TPL))
B_colors <- c("#112bf0", "#ff5233") #blue and red for normalized colors
