
# __author__ = "Aivett Bilbao, aivett.bilbao@pnnl.gov"
# This script was implemented to perform analyses accompanying the manuscript:
# "Online Hydrophilic Interaction Chromatography (HILIC) Enhanced Top-Down Mass Spectrometry Characterization of SARS-Cov-2 Spike Receptor Binding Domain".
# The code reads a list of deconvoluted peaks from intact protein MS analyses (matrix of mass, abundance, and elution time slice using Protein Metrics Intact Mass software version 4.2)
# and performs comparisons, filtering, statistical calculations and generates figures and csv reports.

# install.packages("ggplot2")
# install.packages("ggVennDiagram")

dataBasePath = "E:/glyco/a_Rscripts/Intact-peak-analysis/data"

#dataFolder = "Ray_N331Q"
#inputReconstruction = "N331Q_glycopeptide_reconstruction.csv"

#dataFolder = "Ray_WT"
#inputReconstruction = "RayWT_Nglycan_glycopeptide_reconstruction.csv"

#dataFolder = "Sino_N501Y"
#inputReconstruction = "SinoN501_site4_glycopeptide_reconstruction.csv"

dataFolder = "Sino_WT"
inputReconstruction = "SinoWT_site4_glycopeptide_reconstruction.csv"

# Regex to extract each LC method name from the csv file name
regexMethod = ".*__(.+)__.*" # e.g.

setwd(file.path(dataBasePath, dataFolder))

masstol = 2
minReplicates = 2

outputFile = paste(dataFolder, "_massTol-", masstol, "_minReps-", minReplicates, "_figures.pdf", sep = "")
outputPeaks = paste(dataFolder,"_massTol-", masstol, "_minReps-", minReplicates, "_outputPeaks.csv", sep = "")

allfiles = list.files(pattern = ".csv")
if(length(grep("output", allfiles)) > 0)
  allfiles = allfiles[-(grep("output", allfiles))]

if(length(grep("reconstruction", allfiles)) > 0)
  allfiles = allfiles[-(grep("reconstruction", allfiles))]

if(inputReconstruction != "" && grepl("\\.csv", inputReconstruction))
  recons = read.csv(file = inputReconstruction, sep = ',', stringsAsFactors = FALSE)

# -----------------------------------------------------------
# Create a data frame in long shape (LC methods and technical replicates as multiple rows)
dat = NULL
for(k in 1:length(allfiles))
{
  temp = read.csv(file = allfiles[k], sep = ',', stringsAsFactors = FALSE)
  temp = temp[-(1:2),] # remove first 2 rows
  
  print(allfiles[k])
  
  # Parse the replicates and RT slices from the abundance columns
  tb = NULL
  rtx = 1
  replicate_previous = 1
  for(c in 3:length(temp[1,]))
  {
    tbx = temp[,c(1,c)]
    tbx = tbx[!is.na(tbx[,2]),]
    colnames(tbx) = c("Mass", "Abundance")

    strColname = colnames(temp)[c]
    strColname = sub("\\.[0-9]+$", "", strColname)
    strColname = sub("_[a-z]+$", "", strColname)
    replicatenum = as.numeric(sub(".+_([0-9]+)$", "\\1", strColname))
    tbx$Replicate = replicatenum
    
    if(replicatenum > replicate_previous)
      rtx = 1
    tbx$RT = rtx
    rtx = rtx + 1
    replicate_previous = replicatenum
    tb = rbind(tb, tbx)
  }
  
  tb$Method = sub(regexMethod, "\\1", allfiles[k])
  dat = rbind(dat, tb)
}
tb = NULL
tbx = NULL
temp = NULL

# -----------------------------------------------------------
# Cluster and count replicates per feature, within each method:
dat = dat[with(dat, order(Method, -Abundance)), ]
dat$featID = 1:length(dat[,1])
dat$Count = 0
dat$Cluster = 0
dat$Mass = as.numeric(dat$Mass)
dat = dat[dat$Abundance > 0, ]
for(k in 1:length(dat[,1]))
{
  indexes = which(dat$Count == 0 &
                  dat$Method == dat$Method[k] &
                  abs(dat$Mass - dat$Mass[k]) <= masstol)
  if(length(indexes) > 0)
  {
    dat$Count[indexes] = length(unique(dat$Replicate[indexes]))
    dat$Cluster[indexes] = k
  }
}

# Remove repetitions, keep most abundant feature per replicate, calculate mean abundance and mass:
dat = dat[with(dat, order(Method, Cluster, Replicate, -Abundance)), ]
dat = dat[!duplicated(dat[,c("Method", "Cluster", "Replicate")]), ]
dat$AvgAbundance = ave(dat$Abundance, dat$Method, dat$Cluster, FUN=mean)
dat$Abundance = dat$AvgAbundance
dat$MedianMass = ave(dat$Mass, dat$Method, dat$Cluster, FUN=median)
dat$Mass = floor(dat$MedianMass)
dat$AvgAbundance = NULL
dat$MedianMass = NULL
dat = dat[!duplicated(dat[,c("Method", "Cluster")]), ]

datunfiltered = dat
dat = dat[which(dat$Count >= minReplicates), ] # keep features found in n replicates
dat$Replicate = NULL
dat$Cluster = NULL
dat$Count = NULL

# -----------------------------------------------------------
# Cluster masses across methods to make Venn diagram, but do not remove repetitions:
dat = dat[with(dat, order(-Abundance)), ]
dat$featID = 1:length(dat[,1])
dat$CountMethods = 0
dat$Cluster = 0
for(k in 1:length(dat[,1]))
{
  indexes = which(dat$CountMethods == 0 &
                    abs(dat$Mass - dat$Mass[k]) <= masstol)
  if(length(indexes) > 0)
  {
    dat$CountMethods[indexes] = length(unique(dat$Method[indexes]))
    dat$Cluster[indexes] = k
  }
}
dat = dat[with(dat, order(Method, Cluster)), ]
dat$Mass_original = dat$Mass
dat$MedianMass = ave(dat$Mass, dat$Cluster, FUN=median)
dat$Mass = floor(dat$MedianMass)
dat$featID = NULL
dat$Cluster = NULL
# Remove repetitions after taking floor of Mass:
dat = dat[with(dat, order(Method, Mass, -Abundance)), ]
dat = dat[!duplicated(dat[,c("Method", "Mass")]), ]

# -----------------------------------------------------------
# Create a matrix to format results in wide shape (LC methods as multiple columns): 
widedat = NULL
for(xmethod in unique(dat$Method))
{
  x = dat[which(dat$Method == xmethod), ]
  x = x[,c("Mass", "Abundance", "RT", "Mass_original")]
  colnames(x)[2:4] = paste(colnames(x)[2:4], xmethod, sep = "_")
  if(is.null(widedat))
  { 
    widedat = x
  }else{ 
    widedat = merge(widedat, x, by="Mass", all = TRUE)
  }
}
widedat = cbind(data.frame(Mass=widedat$Mass), widedat[, sort(colnames(widedat)[-1])])

# -----------------------------------------------------------
# Add annotations from reconstructed spectrum:
recons$mass = floor(recons$mass) # binning masses for consistency with integer deconvoluted masses
recons = recons[with(recons, order(mass, -probability)), ]
recons = recons[!duplicated(recons[,c("mass")]), ] # keep the row with highest probability (most abundant)

annotdat = data.frame(Mass=widedat$Mass, ReconstMass = 0, ReconstAbund = 0, ReconstAnnotation="", Methods = "")
if(inputReconstruction != "")
  for(k in 1:length(annotdat[,1]))
  {
    indexes = which(abs(recons$mass - annotdat$Mass[k]) <= masstol)
    if(length(indexes) > 0)
    {
      # Take the most abundant reconstructed peak within mass tolerance:
      maxIndex = which.max(recons$probability[indexes])
      annotdat$ReconstMass[k] = recons$mass[indexes][maxIndex]
      annotdat$ReconstAbund[k] = recons$probability[indexes][maxIndex]
      annotdat$ReconstAnnotation[k] = recons$Glycan[indexes][maxIndex]
    }
  }

# -----------------------------------------------------------
# Reformat method name by concatenating the names if found multiple ones:
colindexes = grep("Abundance", colnames(widedat))
for(k in 1:length(annotdat[,1]))
{
  methods = paste(sub("Abundance_", "", colnames(widedat)[colindexes][!is.na(widedat[k,colindexes])]), collapse = " & ")
  annotdat$Methods[k] = methods
}

widedat = merge(annotdat, widedat, by = "Mass")
write.csv(widedat, file = outputPeaks, row.names = FALSE)

dat$MassString = as.character(dat$Mass)

# Plot figures:
library(ggplot2)

# Save plots as pdf:
pdf(outputFile, paper="a4r", useDingbats=FALSE) # create .pdf file

p = ggplot(dat, aes(x=Mass, y=Abundance, colour=Method, fill=Method))
p = p + theme_bw()
p = p + theme(text=element_text(size = 12), axis.title = element_text(size = 12))
p = p + facet_grid(Method ~ ., scales = "free_y")
plot(p + geom_bar(stat = "identity", size = 0.5))

# Plot peaks that were filtered out for QC:
p = ggplot(datunfiltered, aes(x=factor(Count), y=log10(Abundance), fill=Method))
p = p + theme_bw()
p = p + xlab("Number of replicates")
p = p + theme(text=element_text(size = 12), axis.title = element_text(size = 12))
p = p + facet_grid(. ~ Method)
plot(p + geom_violin(colour="black", draw_quantiles = c(0.25, 0.5, 0.75)))

p = ggplot(datunfiltered, aes(x=factor(Count), y=log10(Abundance), fill=Method))
p = p + theme_bw()
p = p + xlab("Number of replicates")
p = p + theme(text=element_text(size = 12), axis.title = element_text(size = 12))
p = p + facet_grid(. ~ Method)
plot(p + geom_boxplot(colour="black"))

p = ggplot(dat, aes(x=RT, y=Mass, colour=Method))
p = p + theme_bw()
p = p + theme(text=element_text(size = 12), axis.title = element_text(size = 12))
p = p + facet_grid(Method ~ .)
plot(p + geom_point(size = 1.5), alpha = 0.7)

# Factor dividing abundance by range in 5 intervals:
dat$Log10Abundance = cut(log10(dat$Abundance), breaks=5)
p = ggplot(dat, aes(x=RT, y=Mass, colour=Method))
p = p + theme_bw()
p = p + theme(text=element_text(size = 12), axis.title = element_text(size = 12))
p = p + facet_grid(Method ~ .)
plot(p + geom_point(aes(size = Log10Abundance), alpha = 0.7))

p = ggplot(dat, aes(x=RT, y=Mass, colour=Method))
p = p + theme_bw()
p = p + theme(text=element_text(size = 12), axis.title = element_text(size = 12))
p = p + facet_grid(Method ~ .)
plot(p + geom_point(aes(alpha = Log10Abundance)))

p = ggplot(dat, aes(x=RT, y=Mass, colour=Log10Abundance))
p = p + theme_bw()
p = p + theme(text=element_text(size = 12), axis.title = element_text(size = 12))
p = p + facet_grid(Method ~ .)
p = p + scale_color_brewer()
plot(p + geom_point(alpha = 0.7))

p = ggplot(dat, aes(x=log10(Abundance), fill=Method))
p = p + theme_bw()
p = p + theme(text=element_text(size = 12), axis.title = element_text(size = 12))
p = p + facet_grid(Method ~ ., scales = "free_y")
plot(p + geom_histogram(bins = 30, colour = "black"))

p = ggplot(dat, aes(x=Mass, fill=Method))
p = p + theme_bw()
p = p + theme(text=element_text(size = 12), axis.title = element_text(size = 12))
p = p + facet_grid(Method ~ ., scales = "free_y")
plot(p + geom_histogram(bins = 30, colour = "black"))

# -----------------------------------------------------------
# Calculate Pearson correlation:
pearsoncor = NULL
colindexes = grep("Abundance", colnames(widedat))
for(k in colindexes)
{
  rowindexes = which(widedat[,k] > 0)
  x = widedat[rowindexes,k]
  y = widedat$ReconstAbund[rowindexes]
  corresult = cor(x, y, use = "everything",
      method = "pearson")
  pearsoncor = rbind(pearsoncor, data.frame(Method = sub("Abundance_", "", colnames(widedat)[k]), 
                                            PearsonCorrelation = corresult,
                                            CountAnnotPeaks = length(which(widedat$ReconstAbund[rowindexes] > 0))))
}
pearsoncor$textlabel = format(pearsoncor$PearsonCorrelation, nsmall = 0, digits=5, scientific = FALSE)
pearsoncor$textlabel = paste(pearsoncor$textlabel, "\n(", pearsoncor$CountAnnotPeaks, " annot peaks)", sep = '')
p = ggplot(pearsoncor, aes(x=Method, y=PearsonCorrelation, fill=Method))
p = p + theme_bw()
p = p + theme(text=element_text(size = 12), axis.title = element_text(size = 12))
plot(p + geom_bar(stat = "identity", colour="black")
     + geom_text(aes(x=Method, y=PearsonCorrelation+0.04, label = textlabel), 
                 color="blue", position=position_dodge(.2), hjust=.5))

library(ggVennDiagram)

# List of items
myMethods = unique(dat$Method)

x = replicate(length(myMethods), NULL)
for(k in 1:length(myMethods))
{
  x[[k]] = (dat$MassString[which(dat$Method == myMethods[k])])
}
x = setNames(x, as.list(myMethods))

# 4D Venn diagram
ggVennDiagram(x)


dev.off() # close .pdf file

