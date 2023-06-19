# This script is evaluating and plotting the Transient versus Asymptotic characteristics of 
# tropical and subtropical coral populations from Japan and Australia.
# This script utilises the outputs generated using the script 'Constructing IPMs'.
# To speed up processing times, the outputs from this previous script where generated iteratively and so need to be condense back together. This process is handled in STEP 1 of this current script.
# whilst STEP 2 carries out the subsequent analysis of this data.

# Author: James Cant (jic2@st-andrews.ac.uk)
# Date last modified: May 2021

# ----------------------------------------------------------------------------------

# clear workspace
rm(list=ls(all=TRUE))

# load required packages
library(popdemo)
library(gmodels)
library(reshape2)
library(caret)
library(ggplot2)
library(plotrix)
library(Rmisc)
library(plyr)
library(dplyr)
library(plsdepot)
library(visdat)
library(fuzzySim)


# set working directory for data extraction
setwd("Data_file_location")

############################################
# STEP 1: Extract and condense the transient and asymptotic outputs for analyses
############################################

# define an indexing vector for condensing all the data below.
index_start <- seq(1,1000, by = 20)
index_end <- index_start + 19
# row indexing for separating extracted matrices
row_index_start <- seq(1,4000, by = 200); row_index_end <- row_index_start + 199

# 1. Extract already calculated measures
# Generate storage outputs
C.outputs <- ST.outputs <- W.outputs <- matrix(NA, nrow = 1000, ncol = 36)
colnames(W.outputs) <- colnames(ST.outputs) <- colnames(C.outputs) <- c("AS Lambda", "AT Lambda", "JS Lambda", "JT Lambda",
                                                                        "AS U.Kreiss", "AT U.Kreiss", "JS U.Kreiss", "JT U.Kreiss",
                                                                        "AS L.Kreiss", "AT L.Kreiss", "JS L.Kreiss", "JT L.Kreiss",
                                                                        "AS TE", "AT TE", "JS TE", "JT TE",
                                                                        "R.AS", "R.AT", "R.JS", "R.JT",
                                                                        "AS reac", "AT reac", "JS reac", "JT reac",
                                                                        "AS att", "AT att", "JS att", "JT att",
                                                                        "AS maxamp", "AT maxamp", "JS maxamp", "JT maxamp",
                                                                        "AS maxatt", "AT maxatt", "JS maxatt", "JT maxatt")

# Extract data
data_list <- data_list1 <- data_list2 <- list()
data_list[[1]] <- read.csv("C outputs_1-20.csv"); data_list1[[1]] <- read.csv("ST outputs_1-20.csv"); data_list2[[1]] <- read.csv("W outputs_1-20.csv")
data_list[[2]] <- read.csv("C outputs_21-40.csv"); data_list1[[2]] <- read.csv("ST outputs_21-40.csv"); data_list2[[2]] <- read.csv("W outputs_21-40.csv")
data_list[[3]] <- read.csv("C outputs_41-60.csv"); data_list1[[3]] <- read.csv("ST outputs_41-60.csv"); data_list2[[3]] <- read.csv("W outputs_41-60.csv")
data_list[[4]] <- read.csv("C outputs_61-80.csv"); data_list1[[4]] <- read.csv("ST outputs_61-80.csv"); data_list2[[4]] <- read.csv("W outputs_61-80.csv")
data_list[[5]] <- read.csv("C outputs_81-100.csv"); data_list1[[5]] <- read.csv("ST outputs_81-100.csv"); data_list2[[5]] <- read.csv("W outputs_81-100.csv")
data_list[[6]] <- read.csv("C outputs_101-120.csv"); data_list1[[6]] <- read.csv("ST outputs_101-120.csv"); data_list2[[6]] <- read.csv("W outputs_101-120.csv")
data_list[[7]] <- read.csv("C outputs_121-140.csv"); data_list1[[7]] <- read.csv("ST outputs_121-140.csv"); data_list2[[7]] <- read.csv("W outputs_121-140.csv")
data_list[[8]] <- read.csv("C outputs_141-160.csv"); data_list1[[8]] <- read.csv("ST outputs_141-160.csv"); data_list2[[8]] <- read.csv("W outputs_141-160.csv")
data_list[[9]] <- read.csv("C outputs_161-180.csv"); data_list1[[9]] <- read.csv("ST outputs_161-180.csv"); data_list2[[9]] <- read.csv("W outputs_161-180.csv")
data_list[[10]] <- read.csv("C outputs_181-200.csv"); data_list1[[10]] <- read.csv("ST outputs_181-200.csv"); data_list2[[10]] <- read.csv("W outputs_181-200.csv")
data_list[[11]] <- read.csv("C outputs_201-220.csv"); data_list1[[11]] <- read.csv("ST outputs_201-220.csv"); data_list2[[11]] <- read.csv("W outputs_201-220.csv")
data_list[[12]] <- read.csv("C outputs_221-240.csv"); data_list1[[12]] <- read.csv("ST outputs_221-240.csv"); data_list2[[12]] <- read.csv("W outputs_221-240.csv")
data_list[[13]] <- read.csv("C outputs_241-260.csv"); data_list1[[13]] <- read.csv("ST outputs_241-260.csv"); data_list2[[13]] <- read.csv("W outputs_241-260.csv")
data_list[[14]] <- read.csv("C outputs_261-280.csv"); data_list1[[14]] <- read.csv("ST outputs_261-280.csv"); data_list2[[14]] <- read.csv("W outputs_261-280.csv")
data_list[[15]] <- read.csv("C outputs_281-300.csv"); data_list1[[15]] <- read.csv("ST outputs_281-300.csv"); data_list2[[15]] <- read.csv("W outputs_281-300.csv")
data_list[[16]] <- read.csv("C outputs_301-320.csv"); data_list1[[16]] <- read.csv("ST outputs_301-320.csv"); data_list2[[16]] <- read.csv("W outputs_301-320.csv")
data_list[[17]] <- read.csv("C outputs_321-340.csv"); data_list1[[17]] <- read.csv("ST outputs_321-340.csv"); data_list2[[17]] <- read.csv("W outputs_321-340.csv")
data_list[[18]] <- read.csv("C outputs_341-360.csv"); data_list1[[18]] <- read.csv("ST outputs_341-360.csv"); data_list2[[18]] <- read.csv("W outputs_341-360.csv")
data_list[[19]] <- read.csv("C outputs_361-380.csv"); data_list1[[19]] <- read.csv("ST outputs_361-380.csv"); data_list2[[19]] <- read.csv("W outputs_361-380.csv")
data_list[[20]] <- read.csv("C outputs_381-400.csv"); data_list1[[20]] <- read.csv("ST outputs_381-400.csv"); data_list2[[20]] <- read.csv("W outputs_381-400.csv")
data_list[[21]] <- read.csv("C outputs_401-420.csv"); data_list1[[21]] <- read.csv("ST outputs_401-420.csv"); data_list2[[21]] <- read.csv("W outputs_401-420.csv")
data_list[[22]] <- read.csv("C outputs_421-440.csv"); data_list1[[22]] <- read.csv("ST outputs_421-440.csv"); data_list2[[22]] <- read.csv("W outputs_421-440.csv")
data_list[[23]] <- read.csv("C outputs_441-460.csv"); data_list1[[23]] <- read.csv("ST outputs_441-460.csv"); data_list2[[23]] <- read.csv("W outputs_441-460.csv")
data_list[[24]] <- read.csv("C outputs_461-480.csv"); data_list1[[24]] <- read.csv("ST outputs_461-480.csv"); data_list2[[24]] <- read.csv("W outputs_461-480.csv")
data_list[[25]] <- read.csv("C outputs_481-500.csv"); data_list1[[25]] <- read.csv("ST outputs_481-500.csv"); data_list2[[25]] <- read.csv("W outputs_481-500.csv")
data_list[[26]] <- read.csv("C outputs_501-520.csv"); data_list1[[26]] <- read.csv("ST outputs_501-520.csv"); data_list2[[26]] <- read.csv("W outputs_501-520.csv")
data_list[[27]] <- read.csv("C outputs_521-540.csv"); data_list1[[27]] <- read.csv("ST outputs_521-540.csv"); data_list2[[27]] <- read.csv("W outputs_521-540.csv")
data_list[[28]] <- read.csv("C outputs_541-560.csv"); data_list1[[28]] <- read.csv("ST outputs_541-560.csv"); data_list2[[28]] <- read.csv("W outputs_541-560.csv")
data_list[[29]] <- read.csv("C outputs_561-580.csv"); data_list1[[29]] <- read.csv("ST outputs_561-580.csv"); data_list2[[29]] <- read.csv("W outputs_561-580.csv")
data_list[[30]] <- read.csv("C outputs_581-600.csv"); data_list1[[30]] <- read.csv("ST outputs_581-600.csv"); data_list2[[30]] <- read.csv("W outputs_581-600.csv")
data_list[[31]] <- read.csv("C outputs_601-620.csv"); data_list1[[31]] <- read.csv("ST outputs_601-620.csv"); data_list2[[31]] <- read.csv("W outputs_601-620.csv")
data_list[[32]] <- read.csv("C outputs_621-640.csv"); data_list1[[32]] <- read.csv("ST outputs_621-640.csv"); data_list2[[32]] <- read.csv("W outputs_621-640.csv")
data_list[[33]] <- read.csv("C outputs_641-660.csv"); data_list1[[33]] <- read.csv("ST outputs_641-660.csv"); data_list2[[33]] <- read.csv("W outputs_641-660.csv")
data_list[[34]] <- read.csv("C outputs_661-680.csv"); data_list1[[34]] <- read.csv("ST outputs_661-680.csv"); data_list2[[34]] <- read.csv("W outputs_661-680.csv")
data_list[[35]] <- read.csv("C outputs_681-700.csv"); data_list1[[35]] <- read.csv("ST outputs_681-700.csv"); data_list2[[35]] <- read.csv("W outputs_681-700.csv")
data_list[[36]] <- read.csv("C outputs_701-720.csv"); data_list1[[36]] <- read.csv("ST outputs_701-720.csv"); data_list2[[36]] <- read.csv("W outputs_701-720.csv")
data_list[[37]] <- read.csv("C outputs_721-740.csv"); data_list1[[37]] <- read.csv("ST outputs_721-740.csv"); data_list2[[37]] <- read.csv("W outputs_721-740.csv")
data_list[[38]] <- read.csv("C outputs_741-760.csv"); data_list1[[38]] <- read.csv("ST outputs_741-760.csv"); data_list2[[38]] <- read.csv("W outputs_741-760.csv")
data_list[[39]] <- read.csv("C outputs_761-780.csv"); data_list1[[39]] <- read.csv("ST outputs_761-780.csv"); data_list2[[39]] <- read.csv("W outputs_761-780.csv")
data_list[[40]] <- read.csv("C outputs_781-800.csv"); data_list1[[40]] <- read.csv("ST outputs_781-800.csv"); data_list2[[40]] <- read.csv("W outputs_781-800.csv")
data_list[[41]] <- read.csv("C outputs_801-820.csv"); data_list1[[41]] <- read.csv("ST outputs_801-820.csv"); data_list2[[41]] <- read.csv("W outputs_801-820.csv")
data_list[[42]] <- read.csv("C outputs_821-840.csv"); data_list1[[42]] <- read.csv("ST outputs_821-840.csv"); data_list2[[42]] <- read.csv("W outputs_821-840.csv")
data_list[[43]] <- read.csv("C outputs_841-860.csv"); data_list1[[43]] <- read.csv("ST outputs_841-860.csv"); data_list2[[43]] <- read.csv("W outputs_841-860.csv")
data_list[[44]] <- read.csv("C outputs_861-880.csv"); data_list1[[44]] <- read.csv("ST outputs_861-880.csv"); data_list2[[44]] <- read.csv("W outputs_861-880.csv")
data_list[[45]] <- read.csv("C outputs_881-900.csv"); data_list1[[45]] <- read.csv("ST outputs_881-900.csv"); data_list2[[45]] <- read.csv("W outputs_881-900.csv")
data_list[[46]] <- read.csv("C outputs_901-920.csv"); data_list1[[46]] <- read.csv("ST outputs_901-920.csv"); data_list2[[46]] <- read.csv("W outputs_901-920.csv")
data_list[[47]] <- read.csv("C outputs_921-940.csv"); data_list1[[47]] <- read.csv("ST outputs_921-940.csv"); data_list2[[47]] <- read.csv("W outputs_921-940.csv")
data_list[[48]] <- read.csv("C outputs_941-960.csv"); data_list1[[48]] <- read.csv("ST outputs_941-960.csv"); data_list2[[48]] <- read.csv("W outputs_941-960.csv")
data_list[[49]] <- read.csv("C outputs_961-980.csv"); data_list1[[49]] <- read.csv("ST outputs_961-980.csv"); data_list2[[49]] <- read.csv("W outputs_961-980.csv")
data_list[[50]] <- read.csv("C outputs_981-1000.csv"); data_list1[[50]] <- read.csv("ST outputs_981-1000.csv"); data_list2[[50]] <- read.csv("W outputs_981-1000.csv")

# Condense the data together
for(ii in 1:length(data_list)){
  start <- index_start[ii]
  end <- index_end[ii]
  index <- seq(start,end,by = 1)
  # extract the required data
  for(x in index) {
    try(C.outputs[x, ] <- unlist(data_list[[ii]][x, ]))
    try(ST.outputs[x, ] <- unlist(data_list1[[ii]][x, ]))
    try(W.outputs[x, ] <- unlist(data_list2[[ii]][x, ]))
  }
}

# 2. Now evaluate damping ratio estimates for each model.
damping.ratios <- matrix(NA, nrow = 1000, ncol = 12)
colnames(damping.ratios) <- c("CAS","CAT","CJS","CJT",
                              "STAS","STAT","STJS","STJT",
                              "WAS","WAT","WJS","WJT")

# Loop through all K kernels to estimate damping ratios
for (ii in 1:50) {

  # read in kernel
  try(CAS.Kernel <- read.csv(paste0("CAS K kernel_", index_start[ii], "-", index_end[ii], ".csv", sep = "")))
  try(CAT.Kernel <- read.csv(paste0("CAT K kernel_", index_start[ii], "-", index_end[ii], ".csv", sep = "")))
  try(CJS.Kernel <- read.csv(paste0("CJS K kernel_", index_start[ii], "-", index_end[ii], ".csv", sep = "")))
  try(CJT.Kernel <- read.csv(paste0("CJT K kernel_", index_start[ii], "-", index_end[ii], ".csv", sep = "")))
  try(STAS.Kernel <- read.csv(paste0("STAS K kernel_", index_start[ii], "-", index_end[ii], ".csv", sep = "")))
  try(STAT.Kernel <- read.csv(paste0("STAT K kernel_", index_start[ii], "-", index_end[ii], ".csv", sep = "")))
  try(STJS.Kernel <- read.csv(paste0("STJS K kernel_", index_start[ii], "-", index_end[ii], ".csv", sep = "")))
  try(STJT.Kernel <- read.csv(paste0("STJT K kernel_", index_start[ii], "-", index_end[ii], ".csv", sep = "")))
  try(WAS.Kernel <- read.csv(paste0("wAS K kernel_", index_start[ii], "-", index_end[ii], ".csv", sep = "")))
  try(WAT.Kernel <- read.csv(paste0("wAT K kernel_", index_start[ii], "-", index_end[ii], ".csv", sep = "")))
  try(WJS.Kernel <- read.csv(paste0("wJS K kernel_", index_start[ii], "-", index_end[ii], ".csv", sep = "")))
  try(WJT.Kernel <- read.csv(paste0("wJT K kernel_", index_start[ii], "-", index_end[ii], ".csv", sep = "")))
  
  # determine the rows being dealt with 
  index <- seq(index_start[ii], index_end[ii], by = 1)
    
  # extract the required data
  for(xx in index) {
    # progress read out 
    cat(ii,".",xx)
    
    # find the desired matrix
    row.index <- which(index == xx)
    row.start <- row_index_start[row.index]
    row.end <- row_index_end[row.index]
    
    # estimate and store the damping ratio
    try(damping.ratios[xx, "CAS"] <- dr(as.matrix(CAS.Kernel[row.start:row.end,])))
    try(damping.ratios[xx, "CAT"] <- dr(as.matrix(CAT.Kernel[row.start:row.end,])))
    try(damping.ratios[xx, "CJS"] <- dr(as.matrix(CJS.Kernel[row.start:row.end,])))
    try(damping.ratios[xx, "CJT"] <- dr(as.matrix(CJT.Kernel[row.start:row.end,])))
    try(damping.ratios[xx, "STAS"] <- dr(as.matrix(STAS.Kernel[row.start:row.end,])))
    try(damping.ratios[xx, "STAT"] <- dr(as.matrix(STAT.Kernel[row.start:row.end,])))
    try(damping.ratios[xx, "STJS"] <- dr(as.matrix(STJS.Kernel[row.start:row.end,])))
    try(damping.ratios[xx, "STJT"] <- dr(as.matrix(STJT.Kernel[row.start:row.end,])))
    try(damping.ratios[xx, "WAS"] <- dr(as.matrix(WAS.Kernel[row.start:row.end,])))
    try(damping.ratios[xx, "WAT"] <- dr(as.matrix(WAT.Kernel[row.start:row.end,])))
    try(damping.ratios[xx, "WJS"] <- dr(as.matrix(WJS.Kernel[row.start:row.end,])))
    try(damping.ratios[xx, "WJT"] <- dr(as.matrix(WJT.Kernel[row.start:row.end,])))
    # Also for some reason the STJS kernels were missed in estimates of maximal amplification
    # so that can be addressed here.
    #try(ST.outputs[xx, "JS maxamp"] <- maxamp(as.matrix(STJS.Kernel[row.start:row.end,])))
  }
}

#### and remove un-needed objects to save memory space
rm(data_list, data_list1, data_list2, x, start, end, ii, index, index_end, index_start,
   row_index_end, row_index_start, row.index, row.end, row.start, xx, WAS.Kernel, WAT.Kernel, WJS.Kernel, WJT.Kernel,
   CAS.Kernel, CAT.Kernel, CJS.Kernel, CJT.Kernel, STAS.Kernel, STAT.Kernel, STJS.Kernel, STJT.Kernel)

#### Save outputs as a checkpoint before moving onto STEP 2

# switch working directory
setwd("save_file_location")

# and save output files
write.csv(C.outputs, file = "Competitive outputs.csv", row.names = FALSE)
write.csv(ST.outputs, file = "Stress Tolerant outputs.csv", row.names = FALSE) # checkpoint
write.csv(W.outputs, file = "Weedy outputs.csv", row.names = FALSE)
write.csv(damping.ratios, file = "Damping ratios.csv", row.names = FALSE)

# Now add the maximal amplification estimates that have been calculated separately for the Stress Tolerant population
# first define the indexing used
index_start <- seq(1,49, by = 2)
index_end <- seq(2,50, by = 2)
row_start <- seq(1,961, by = 40)
row_end <- seq(40,1000, by = 40)

# now run a loop to do the extracting and merging.
for (ii in 1:25) {
  
  # read in the data
  try(STJS_maxamp <- read.csv(paste0("Stress Tolerant outputs_", index_start[ii], "-", index_end[ii], ".csv", sep = "")))
  
  # determine the rows that have data 
  index <- seq(row_start[ii], row_end[ii], by = 1)
  
  # extract the required data
  for(xx in index) {
    # progress read out 
    cat(ii,".",xx)
    
    # extract and store the maximal amplification estimate
    try(ST.outputs[xx, "JS maxamp"] <- STJS_maxamp[xx,])
  }
}

#### and again remove un-needed objects to save memory space
rm(ii, index, index_end, index_start, row_end, row_start, xx, STJS_maxamp)

# I also know from previous run throughs that one estimate hasn't worked and is clearly an outlier in the distribution
hist(ST.outputs[, "JS maxamp"])
hist(ST.outputs[which(ST.outputs[, "JS maxamp"] <= 1.0e+13), "JS maxamp"])
ST.outputs[which(ST.outputs[, "JS maxamp"] >= 1.0e+13), "JS maxamp"] <- NA

# and finally restore the now complete datafile
write.csv(ST.outputs, file = "Stress Tolerant outputs.csv", row.names = FALSE)

##########################################
# STEP 2: Evaluate spatial patterns in transient and asymptotic characteristics
##########################################

###### If no new data needed start here
setwd("save_file_location")

# Reload output files if needed
C.outputs <- read.csv("Competitive outputs.csv")
ST.outputs <- read.csv("Stress Tolerant outputs.csv")
W.outputs <- read.csv("Weedy outputs.csv")
damping.ratios <- read.csv("Damping ratios.csv")

# load GPS information
GPS_coords <- read.csv("Site SST info.csv", row.names = 1) # Generated using the 'SST variance extraction and plotting' script

# How does lambda compare across the groups
CI(na.omit(C.outputs[,"AS.Lambda"]))
CI(na.omit(C.outputs[,"AT.Lambda"]))
CI(na.omit(C.outputs[,"JS.Lambda"]))
CI(na.omit(C.outputs[,"JT.Lambda"]))
CI(na.omit(ST.outputs[,"AS.Lambda"]))
CI(na.omit(ST.outputs[,"AT.Lambda"]))
CI(na.omit(ST.outputs[,"JS.Lambda"]))
CI(na.omit(ST.outputs[,"JT.Lambda"]))
CI(na.omit(W.outputs[,"AS.Lambda"]))
CI(na.omit(W.outputs[,"AT.Lambda"]))
CI(na.omit(W.outputs[,"JS.Lambda"]))
CI(na.omit(W.outputs[,"JT.Lambda"]))

# And how does maximal amplification compare across the groups
CI(na.omit(C.outputs[,"AS.maxamp"]))
CI(na.omit(C.outputs[,"AT.maxamp"]))
CI(na.omit(C.outputs[,"JS.maxamp"]))
CI(na.omit(C.outputs[,"JT.maxamp"]))
CI(na.omit(ST.outputs[,"AS.maxamp"]))
CI(na.omit(ST.outputs[,"AT.maxamp"]))
CI(na.omit(ST.outputs[,"JS.maxamp"]))
CI(na.omit(ST.outputs[,"JT.maxamp"]))
CI(na.omit(W.outputs[,"AS.maxamp"]))
CI(na.omit(W.outputs[,"AT.maxamp"]))
CI(na.omit(W.outputs[,"JS.maxamp"]))
CI(na.omit(W.outputs[,"JT.maxamp"]))

# And how does transient envelope compare across the groups
CI(na.omit(C.outputs[,"AS.TE"]))
CI(na.omit(C.outputs[,"AT.TE"]))
CI(na.omit(C.outputs[,"JS.TE"]))
CI(na.omit(C.outputs[,"JT.TE"]))
CI(na.omit(ST.outputs[,"AS.TE"]))
CI(na.omit(ST.outputs[,"AT.TE"]))
CI(na.omit(ST.outputs[,"JS.TE"]))
CI(na.omit(ST.outputs[,"JT.TE"]))
CI(na.omit(W.outputs[,"AS.TE"]))
CI(na.omit(W.outputs[,"AT.TE"]))
CI(na.omit(W.outputs[,"JS.TE"]))
CI(na.omit(W.outputs[,"JT.TE"]))

## So the overall aim is first test for significant differences in the amplification characteristics of tropical and subtropical populations
## before the exploring trends in the lambda and transient envelope and damping ratios of the different coral groups using a pls regression approach to see how the recovery resistance and performance of coral populations relates to temperature variability.
# Before this I need to reformat the necessary data.

#### 1a. Reformat data -----------------------------------------------------------------
# Variables of interest: Country, Region, Lambda, Damping ratio, Max Amplification, Max Attenuation, 
#                        Attenuation, Reactivity, Transient Envelope, Mean SST & SST CV.

# Define vector for country and region variables
Country <- rep(c("Australia","Japan"), each = 2000)
Region <- rep(rep(c("Subtropical", "Tropical"), each = 1000), length.out = 4000) 
# and for abiotic variables
MeanSST <- rep(c(GPS_coords["AS","Mean"],GPS_coords["AT","Mean"],GPS_coords["JS","Mean"],GPS_coords["JT","Mean"]), each = 1000)
cvSST <- rep(c(GPS_coords["AS","CV"],GPS_coords["AT","CV"],GPS_coords["JS","CV"],GPS_coords["JT","CV"]), each = 1000)
SSTFreq <- rep(c(GPS_coords["AS","Freq"],GPS_coords["AT","Freq"],GPS_coords["JS","Freq"],GPS_coords["JT","Freq"]), each = 1000)
SSTAuto <- rep(c(GPS_coords["AS","Auto"],GPS_coords["AT","Auto"],GPS_coords["JS","Auto"],GPS_coords["JT","Auto"]), each = 1000)

# Competitive
Competitive_lambda <- melt(C.outputs[,c("AS.Lambda", "AT.Lambda", "JS.Lambda", "JT.Lambda")])
Competitive_DR <- melt(damping.ratios[, c("CAS", "CAT", "CJS", "CJT")])
Competitive_TE <- melt(C.outputs[, c("AS.TE", "AT.TE", "JS.TE", "JT.TE")])
Competitive_R <- melt(C.outputs[, c("AS.reac", "AT.reac", "JS.reac", "JT.reac")])
Competitive_Att <- melt(C.outputs[, c("AS.att", "AT.att", "JS.att", "JT.att")])
Competitive_MAmp <- melt(C.outputs[, c("AS.maxamp", "AT.maxamp", "JS.maxamp", "JT.maxamp")])
Competitive_MAtt <- melt(C.outputs[, c("AS.maxatt", "AT.maxatt", "JS.maxatt", "JT.maxatt")])
# condense back together
Competitive_stats <- data.frame(Country, Region, Competitive_lambda[,2], Competitive_DR[,2], Competitive_TE[,2],
                                Competitive_R[,2], Competitive_Att[,2], Competitive_MAmp[,2], Competitive_MAtt[,2],
                                MeanSST, cvSST, SSTFreq, SSTAuto)
# reassign column names
colnames(Competitive_stats) <- c("Country", "Region", "Lambda", "DR", "TE", "Reac", "Att", "MaxAmp", "MaxAtt", "SSTMean", "SSTCV", "SSTFreq", "SSTAuto")
Competitive_stats$LHS <- "C"

# Stress Tolerant
Tolerant_lambda <- melt(ST.outputs[,c("AS.Lambda", "AT.Lambda", "JS.Lambda", "JT.Lambda")])
Tolerant_DR <- melt(damping.ratios[, c("STAS", "STAT", "STJS", "STJT")])
Tolerant_TE <- melt(ST.outputs[, c("AS.TE", "AT.TE", "JS.TE", "JT.TE")])
Tolerant_R <- melt(ST.outputs[, c("AS.reac", "AT.reac", "JS.reac", "JT.reac")])
Tolerant_Att <- melt(ST.outputs[, c("AS.att", "AT.att", "JS.att", "JT.att")])
Tolerant_MAmp <- melt(ST.outputs[, c("AS.maxamp", "AT.maxamp", "JS.maxamp", "JT.maxamp")], id.vars = NULL)
Tolerant_MAtt <- melt(ST.outputs[, c("AS.maxatt", "AT.maxatt", "JS.maxatt", "JT.maxatt")])
# condense back together
Tolerant_stats <- data.frame(Country, Region, Tolerant_lambda[,2], Tolerant_DR[,2],
                             Tolerant_TE[,2], Tolerant_R[,2], Tolerant_Att[,2], Tolerant_MAmp[,2],
                             Tolerant_MAtt[,2], MeanSST, cvSST, SSTFreq, SSTAuto)
# reassign column names
colnames(Tolerant_stats) <- c("Country", "Region", "Lambda", "DR", "TE", "Reac", "Att", "MaxAmp", "MaxAtt", "SSTMean", "SSTCV", "SSTFreq", "SSTAuto")
Tolerant_stats$LHS <- "ST"

# Weedy
Weedy_lambda <- melt(W.outputs[,c("AS.Lambda", "AT.Lambda", "JS.Lambda", "JT.Lambda")])
Weedy_DR <- melt(damping.ratios[, c("WAS", "WAT", "WJS", "WJT")])
Weedy_TE <- melt(W.outputs[, c("AS.TE", "AT.TE", "JS.TE", "JT.TE")])
Weedy_R <- melt(W.outputs[, c("AS.reac", "AT.reac", "JS.reac", "JT.reac")])
Weedy_Att <- melt(W.outputs[, c("AS.att", "AT.att", "JS.att", "JT.att")])
Weedy_MAmp <- melt(W.outputs[, c("AS.maxamp", "AT.maxamp", "JS.maxamp", "JT.maxamp")], id.vars = NULL)
Weedy_MAtt <- melt(W.outputs[, c("AS.maxatt", "AT.maxatt", "JS.maxatt", "JT.maxatt")])
# condense back together
Weedy_stats <- data.frame(Country, Region, Weedy_lambda[,2], Weedy_DR[,2], Weedy_TE[,2],
                          Weedy_R[,2], Weedy_Att[,2], Weedy_MAmp[,2], Weedy_MAtt[,2], 
                          MeanSST, cvSST, SSTFreq, SSTAuto)
# reassign column names
colnames(Weedy_stats) <- c("Country", "Region", "Lambda", "DR", "TE", "Reac", "Att", "MaxAmp", "MaxAtt", "SSTMean", "SSTCV", "SSTFreq", "SSTAuto")
# add LHS variable
Weedy_stats$LHS <- "W"

# stitch together the three dataframes
coral_stats <- do.call(rbind, list(Competitive_stats, Tolerant_stats, Weedy_stats))

# and tidy up
rm(Weedy_lambda, Weedy_DR, Weedy_TE, Weedy_R, Weedy_Att, Weedy_MAmp, Weedy_MAtt, 
   MeanSST, cvSST, SSTFreq, SSTAuto,
   Tolerant_lambda, Tolerant_DR, Tolerant_TE, Tolerant_R, Tolerant_Att, Tolerant_MAmp, Tolerant_MAtt,
   Competitive_lambda, Competitive_DR, Competitive_TE, Competitive_R, Competitive_Att, Competitive_MAmp, Competitive_MAtt, Country, Region)

#### 1b. Transform if necessary -----------------------------------------------------------------
# The main variables being used at this stage are the maximal amplification, lambda and the transient envelope

# a. Lambda 
hist(coral_stats$Lambda) # there are some values here outside the limits expected for lambda - thus these Jackknifed model variants need to be omitted 
                         # from further analyses
dim(coral_stats[which(coral_stats$Lambda > 2 | coral_stats$Lambda < 0),])[1] # this will remove 26 entries
coral_stats[which(coral_stats$Lambda > 2 | coral_stats$Lambda < 0), c("Lambda", "DR", "TE", "Reac", "Att", "MaxAmp", "MaxAtt")] <- NA
hist(coral_stats$Lambda) # much better

# b. Maximal Amplification
hist(coral_stats$MaxAmp) # <- needs work
BoxCoxTrans(coral_stats$MaxAmp, na.rm = TRUE)
hist(log(coral_stats$MaxAmp))
hist(1/coral_stats$MaxAmp^0.5) # much better
# make the change
coral_stats$tMaxAmp <- 1/coral_stats$MaxAmp^0.5

# c. Transient envelope
hist(coral_stats$TE) # <- needs work
BoxCoxTrans(coral_stats$TE, na.rm = TRUE)
hist(1/coral_stats$TE^0.1) # much better
hist(log(coral_stats$TE)) # also works
# make the change
coral_stats$tTE <- 1/(coral_stats$TE^0.1)

# D. Damping ratio
hist(coral_stats$DR) # <- Needs work
BoxCoxTrans(coral_stats$DR, na.rm = TRUE)
hist(log(coral_stats$DR))
hist(1/(coral_stats$DR^2))
#make the change
coral_stats$tDR <- 1/(coral_stats$DR^2)

# convert character variables into factors
coral_stats$Country <- as.factor(coral_stats$Country)
coral_stats$Region <- as.factor(coral_stats$Region)
coral_stats$LHS <- as.factor(coral_stats$LHS)

# 2. Three-way ANOVAs ----------------------------------------------------------------------------------------------- 

# Firstly run a three-way ANOVA comparing population growth rates across the three coral groups, countries and regions
Lambda.mod <- aov(Lambda ~ LHS * Region * Country, data = coral_stats)
summary(Lambda.mod)
TukeyHSD(Lambda.mod)

# and extract the mean and variance estimates for each populations growth characteristics
lambda.stats <- summarySE(coral_stats, measurevar = "Lambda", groupvars = c("Country","Region","LHS"), na.rm = TRUE)

# Now Run a three-way ANOVA comparing Maximal Amplification across the three coral groups, countries and regions
maxamp.mod <- aov(tMaxAmp ~ LHS * Region * Country, data = coral_stats)
summary(maxamp.mod)
TukeyHSD(maxamp.mod)

# Extract the necessary data for plotting the interaction
interaction.plot.data <- coral_stats %>% 
                            group_by(LHS, Country, Region) %>% 
                             dplyr::summarise(tMAmp = mean(tMaxAmp, na.rm = TRUE),
                                              tconf.low = CI(na.omit(tMaxAmp))[3],
                                              tconf.high = CI(na.omit(tMaxAmp))[1],
                                              tSE = std.error(tMaxAmp, na.rm = TRUE),
                                              tSD = sd(tMaxAmp, na.rm = TRUE),
                                              MAmp = mean(MaxAmp, na.rm = TRUE),
                                              conf.low = CI(na.omit(MaxAmp))[3],
                                              conf.high = CI(na.omit(MaxAmp))[1],
                                              SE = std.error(MaxAmp, na.rm = TRUE),
                                              SD = sd(MaxAmp, na.rm = TRUE),)

# Plot the data
interaction.plot <- ggplot(interaction.plot.data, aes(x = Region, y = tMAmp)) +
  # 95% intervals
  geom_errorbar(aes(ymin = tconf.low, ymax = tconf.high), size = 0.9, width=0.15,
                position=position_dodge(0.65)) +
  geom_point(aes(col = LHS, shape = Region), size = 6) +
  geom_line(aes(group = LHS, col = LHS), size = 1.5, linetype = "dashed") +
  facet_wrap(~Country) +
  # Formatting
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0,1)) +
  scale_color_manual(name = "", values = c("blue2", "orange2", "red2"), 
                     labels = c("Competitive    ", "Stress-Tolerant","Weedy          "),
                     guide = guide_legend(override.aes = list(shape = 15, size = 10))) +
  scale_shape_manual(name = "", values = c(16, 17), 
                     labels = c("Subtropical", "Tropical"),
                     guide = NULL) + # tropical = triangle, and subtropical = circle
  xlab("") + 
  ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(size = 20, colour = "black"), axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0.9),axis.line.y.right = element_line(colour = "black", size = 1),
        legend.background = element_blank(), legend.key = element_blank(),legend.position = "bottom",
        legend.box = "horizontal")

# 3. PLS Regression ---------------------------------------------------------------------------------------------------

# Firstly I need to ensure there is not colinearity across my variables
# firstly the demographic variables
multicol(vars = coral_stats[,c("tTE","Lambda","tDR")]) # all Variance inflation factors estimates are <10 and the variable tolerances are > 0.2
# abiotic variables (Obviously the different temporal time series are all colinear so they need to be explored separately)
multicol(vars = coral_stats[,c("SSTMean", "SSTCV", "SSTFreq", "SSTAuto")]) # these values cannot be used together
multicol(vars = coral_stats[,c("SSTMean", "SSTCV", "SSTFreq")]) # these are fine to use
# so it appears the autocorrelation variable is to colinear with the other variables so will need to be dropped.

####### Run regression
# How much data is missing
vis_miss(coral_stats[,c("tTE","Lambda","tDR","SSTMean","SSTCV","SSTFreq")]) # a small amount fo data is missing from the entries set to NA due to Lambda inaccuracies
# remove these 
coral_stats <- coral_stats[complete.cases(coral_stats[,c("tTE","Lambda","tDR","SSTMean","SSTCV","SSTFreq")]),]
# and run the model
coral.pls <- plsreg2(predictors = coral_stats[,c("SSTMean","SSTCV","SSTFreq")], responses = coral_stats[,c("tTE","Lambda","tDR")]) 
# That worked but what relationship has been found (if any)?
plot(coral.pls)

# Extract the necessary data for plotting the interaction (combine with previous dataframe)
coral_stats$xy <- coral.pls$y.scores[,1] # how the y variables map onto the 1st component
coral_stats$yy <- coral.pls$y.scores[,2] # and how the y variables map onto the 2nd components

# Extract loadings
X.loadings <- as.data.frame(coral.pls[["x.loads"]])
Y.loadings <- as.data.frame(coral.pls[["y.loads"]])

# generate plot
PLS.plot <- ggplot(coral_stats, aes(x = xy, y = yy)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_point(aes(col = LHS), size = 2, alpha = 0.6) +
  # add loading arrows
  # X variables
  geom_segment(data = X.loadings, aes(x = 0, y = 0, xend = (p1*8), yend = (p2*8)), arrow = arrow(length = unit(0.9, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (X.loadings[,1]*9), y = (X.loadings[,2]*9),
           #label = c("  Mean SST", "SST CV", "SST Freq"), color ="black", 
           # size = 7, fontface = "bold") +
  # Y variables
  geom_segment(data = Y.loadings[2,], aes(x = 0, y = 0, xend = (c1*22), yend = (c2*22)), arrow = arrow(length = unit(0.9, "picas"), type = "closed"),
               color = "gray30", size = 1) + #Lambda wasn't inversely transformed
  #annotate("text", x = (Y.loadings[2,1]*24), y = (Y.loadings[2,2]*24),
            #label = c("Lambda"), color ="gray30", 
           #size = 7, fontface = "italic")+
  geom_segment(data = Y.loadings[c(1,3),], aes(x = 0, y = 0, xend = (c1*-22), yend = (c2*-22)), arrow = arrow(length = unit(0.9, "picas"), type = "closed"),
               color = "gray30", size = 1) + # Transient envelope and damping ratio were
  #annotate("text", x = (Y.loadings[c(1,3),1]*-22), y = (Y.loadings[c(1,3),2]*-22),
           #label = c("TE", "DR"), color ="gray30", 
           #size = 7, fontface = "italic") +
  # Formatting
  scale_x_continuous(expand = expansion(mult = c(0.05,0.05)), limits = c(-8,13)) +
  scale_y_continuous(expand = expansion(mult = c(0.05,0.05)), limits = c(-12,7)) +
  scale_color_manual(name = "", values = c("blue2", "orange2", "red2"), 
                     labels = c("Competitive    ", "Stress-Tolerant","Weedy          "),
                     guide = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 10))) +
  xlab("") + 
  ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0.9),axis.line.y.right = element_line(colour = "black", size = 1),
        legend.background = element_blank(), legend.key = element_blank(),legend.position = "bottom",
        legend.box = "horizontal")

# And display the model fit
coral.pls$expvar

# 4. Repeat PLS Regression ---------------------------------------------------------------------------------------------------
# On this repeat I just one one value for each measure for each population (the mean obtained through Jackknifing)

# Firstly I need to obtain summaries for each population
# Mean estimates
coral_stats_mean <- aggregate(coral_stats, list(coral_stats$Country, coral_stats$Region, coral_stats$LHS), mean)[,c(1:3,6:16,18:20)]

# Now double check that there is no colinearity across my variables
# firstly the demographic variables
multicol(vars = coral_stats_mean[,c("TE","Lambda","DR")]) # all Variance inflation factors estimates are <10 and the variable tolerances are > 0.2
# abiotic variables (Obviously the different temporal time series are all colinear so they need to be explored separately)
multicol(vars = coral_stats_mean[,c("SSTMean", "SSTCV", "SSTFreq", "SSTAuto")]) # these values cannot be used together
multicol(vars = coral_stats_mean[,c("SSTMean", "SSTCV", "SSTFreq")]) # these are fine to use
# so it appears the autocorrelation variable is to colinear with the other variables so will need to be dropped.

####### Run regression
# How much data is missing
vis_miss(coral_stats_mean[,c("TE","Lambda","DR","SSTMean","SSTCV","SSTFreq")]) # Nothing is missing
# and run the model
coral.pls2 <- plsreg2(predictors = coral_stats_mean[,c("SSTMean","SSTCV","SSTFreq")], responses = coral_stats_mean[,c("TE","Lambda","DR")]) 
# That worked but what relationship has been found (if any)?
plot(coral.pls2)

# Extract the necessary data for plotting the interaction (combine with previous dataframe)
coral_stats_mean$xy <- coral.pls2$y.scores[,1] # how the y variables map onto the 1st component
coral_stats_mean$yy <- coral.pls2$y.scores[,2] # and how the y variables map onto the 2nd components

# Extract loadings
X.loadings2 <- as.data.frame(coral.pls2[["x.loads"]])
Y.loadings2 <- as.data.frame(coral.pls2[["y.loads"]])

# generate plot
#PLS.plot2 <- 
ggplot(coral_stats_mean, aes(x = xy, y = yy)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  stat_ellipse(aes(lty = Group.2, fill = Group.2, color = Group.2), geom = "polygon", size = 1.2, type = "t", alpha = 0.15) +
  geom_point(aes(col = Group.3, shape = Group.2), size = 6, alpha = 0.6) +
  # add loading arrows
  # X variables
  geom_segment(data = X.loadings2, aes(x = 0, y = 0, xend = (p1*6), yend = (p2*6)), arrow = arrow(length = unit(0.8, "picas"), type = "closed"),
               color = "black", size = 1) +
  #annotate("text", x = (X.loadings2[,1]*5), y = (X.loadings2[,2]*5),
  #label = c("  Mean SST", "SST CV", "SST Freq"), color ="black", 
  #size = 7, fontface = "bold") +
  # Y variables
  geom_segment(data = Y.loadings2, aes(x = 0, y = 0, xend = (c1*8), yend = (c2*8)), arrow = arrow(length = unit(0.8, "picas"), type = "closed"),
               color = "gray30", size = 1) + #Lambda wasn't inversely transformed
  #annotate("text", x = (Y.loadings2[,1]*6), y = (Y.loadings2[,2]*6),
  #label = c("TE", "Lambda", "DR"), color ="gray30", 
  #size = 7, fontface = "italic") +
  # Formatting
  scale_x_continuous(expand = expansion(mult = c(0.05,0.05)), limits = c(-4,6)) +
  scale_y_continuous(expand = expansion(mult = c(0.05,0.05)), limits = c(-6,6)) +
  scale_color_manual(name = "", values = c("blue2", "orange2", "gray87", "gray75", "red2"), 
                     labels = c("Competitive","Stress-Tolerant","Subtropical","Tropical","Weedy"),
                     guide = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 10))) +
  scale_shape_manual(name = "", values = c(16, 17), 
                     labels = c("Subtropical", "Tropical"),
                     guide = NULL) + # tropical = triangle, and subtropical = circle
  scale_linetype_manual(name = "", values = c("dashed", "solid"), 
                      labels = c("Subtropical", "Tropical"),
                      guide = NULL) + # tropical = solid, and subtropical = dotted
  scale_fill_manual(name = "", values = c("gray87", "gray75"), 
                    labels = c("Subtropical", "Tropical"),
                    guide = NULL) + # tropical = darker, and subtropical = lighter
  xlab("") + 
  ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(size = 20, colour = "black"), axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0.9),axis.line.y.right = element_line(colour = "black", size = 1),
        legend.background = element_blank(), legend.key = element_blank(),legend.position = "bottom",
        legend.box = "horizontal")

# And display the model fit
coral.pls2$expvar

###################################################################################################################################