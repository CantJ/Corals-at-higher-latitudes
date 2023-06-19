# This script is for evaluating and visualizing the vital rate sensitivities of long-term performance and transient variability seen 
# across the tropical and subtropical coral populations from Japan and Australia.
# This script relies on outputs generated using the 'Computing sensitivities' script

# Author: James Cant (jic2@st-andrews.ac.uk)
# Date last modified: May 2021

# ----------------------------------------------------------------------------------

# clear workspace
rm(list=ls(all=TRUE))

# load required packages
library(popdemo)
library(fields)
library(ggplot2)
library(dplyr)
library(Rmisc)

# Load IPM construction functions
source("file_pathway/IPM construction functions.R")

############################################
# STEP 1: Extract and create sensitivity and elasticity matrices for maximal amplification
############################################

# source relevant ARC output files
setwd("data_file_location")

# define matrix dimensions
ndim <- 200

# create output matrices
# Sensitivity matrices
CAS.MA.sens <- CAT.MA.sens <- CJS.MA.sens <- CJT.MA.sens <- 
  STAS.MA.sens <- STAT.MA.sens <- STJS.MA.sens <- STJT.MA.sens <- 
  WAS.MA.sens <- WAT.MA.sens <- WJS.MA.sens <- WJT.MA.sens <- matrix(NA, ncol = ndim, nrow = ndim)
# Elasticity matrices
CAS.MA.elas <- CAT.MA.elas <- CJS.MA.elas <- CJT.MA.elas <- 
  STAS.MA.elas <- STAT.MA.elas <- STJS.MA.elas <- STJT.MA.elas <- 
  WAS.MA.elas <- WAT.MA.elas <- WJS.MA.elas <- WJT.MA.elas <- matrix(NA, ndim, ndim)
# Mean sensitivity matrices
CA.MA.sens <- CJ.MA.sens <- STA.MA.sens <- STJ.MA.sens <- WA.MA.sens <- WJ.MA.sens <- matrix(NA, ndim, ndim)
# Mean Elasticity matrices
CA.MA.elas <- CJ.MA.elas <- STA.MA.elas <- STJ.MA.elas <- WA.MA.elas <- WJ.MA.elas <- matrix(NA, ndim, ndim)

# define loop indexing
index.start <- seq(1, 2400, by = ndim)
index.end <- seq(ndim, 2400, by = ndim)

# run loop to extract data
for (ii in 1:ndim) {
  # progress countdown
  print((ndim+1)-ii)

  # read in desired column
  coldata <- read.csv(paste0("Sensitivity matrices_col-",ii,".csv", sep = ""))
  
  # extract data
  # Base sensitivities and elasticities
  CAS.MA.sens[,ii] <- coldata[index.start[1]:index.end[1], 1]; CAS.MA.elas[,ii] <- coldata[index.start[1]:index.end[1], 2]
  CAT.MA.sens[,ii] <- coldata[index.start[2]:index.end[2], 1]; CAT.MA.elas[,ii] <- coldata[index.start[2]:index.end[2], 2]
  CJS.MA.sens[,ii] <- coldata[index.start[3]:index.end[3], 1]; CJS.MA.elas[,ii] <- coldata[index.start[3]:index.end[3], 2]
  CJT.MA.sens[,ii] <- coldata[index.start[4]:index.end[4], 1]; CJT.MA.elas[,ii] <- coldata[index.start[4]:index.end[4], 2]
  STAS.MA.sens[,ii] <- coldata[index.start[5]:index.end[5], 1]; STAS.MA.elas[,ii] <- coldata[index.start[5]:index.end[5], 2]
  STAT.MA.sens[,ii] <- coldata[index.start[6]:index.end[6], 1]; STAT.MA.elas[,ii] <- coldata[index.start[6]:index.end[6], 2]
  STJS.MA.sens[,ii] <- coldata[index.start[7]:index.end[7], 1]; STJS.MA.elas[,ii] <- coldata[index.start[7]:index.end[7], 2]
  STJT.MA.sens[,ii] <- coldata[index.start[8]:index.end[8], 1]; STJT.MA.elas[,ii] <- coldata[index.start[8]:index.end[8], 2]
  WAS.MA.sens[,ii] <- coldata[index.start[9]:index.end[9], 1]; WAS.MA.elas[,ii] <- coldata[index.start[9]:index.end[9], 2]
  WAT.MA.sens[,ii] <- coldata[index.start[10]:index.end[10], 1]; WAT.MA.elas[,ii] <- coldata[index.start[10]:index.end[10], 2]
  WJS.MA.sens[,ii] <- coldata[index.start[11]:index.end[11], 1]; WJS.MA.elas[,ii] <- coldata[index.start[11]:index.end[11], 2]
  WJT.MA.sens[,ii] <- coldata[index.start[12]:index.end[12], 1]; WJT.MA.elas[,ii] <- coldata[index.start[12]:index.end[12], 2]
  # Mean sensitivities and elasticities
  CA.MA.sens[,ii] <- coldata[index.start[1]:index.end[1], 3]; CA.MA.elas[,ii] <- coldata[index.start[1]:index.end[1], 4]
  CJ.MA.sens[,ii] <- coldata[index.start[3]:index.end[3], 3]; CJ.MA.elas[,ii] <- coldata[index.start[3]:index.end[3], 4]
  STA.MA.sens[,ii] <- coldata[index.start[5]:index.end[5], 3]; STA.MA.elas[,ii] <- coldata[index.start[5]:index.end[5], 4]
  STJ.MA.sens[,ii] <- coldata[index.start[7]:index.end[7], 3]; STJ.MA.elas[,ii] <- coldata[index.start[7]:index.end[7], 4]
  WA.MA.sens[,ii] <- coldata[index.start[9]:index.end[9], 3]; WA.MA.elas[,ii] <- coldata[index.start[9]:index.end[9], 4]
  WJ.MA.sens[,ii] <- coldata[index.start[11]:index.end[11], 3]; WJ.MA.elas[,ii] <- coldata[index.start[11]:index.end[11], 4]
  
  # tidy up
  rm(coldata)
}


############################################
# STEP 2: Extract and create sensitivity and elasticity matrices for population growth
############################################

# source relevant ARC output files
setwd("F:/James/Loop Outputs")

# Read in K kernels 
try(CAS.K <- read.csv("CAS K kernel_1-20.csv"))
try(CAT.K <- read.csv("CAT K kernel_1-20.csv"))
try(CJS.K <- read.csv("CJS K kernel_1-20.csv"))
try(CJT.K <- read.csv("CJT K kernel_1-20.csv"))
try(STAS.K <- read.csv("STAS K kernel_1-20.csv"))
try(STAT.K <- read.csv("STAT K kernel_1-20.csv"))
try(STJS.K <- read.csv("STJS K kernel_1-20.csv"))
try(STJT.K <- read.csv("STJT K kernel_1-20.csv"))
try(WAS.K <- read.csv("WAS K kernel_1-20.csv"))
try(WAT.K <- read.csv("WAT K kernel_1-20.csv"))
try(WJS.K <- read.csv("WJS K kernel_1-20.csv"))
try(WJT.K <- read.csv("WJT K kernel_1-20.csv"))

# And extract the required kernels (I am only here interested in the first model)
try(CAS.K <- as.matrix(CAS.K[1:200,]))
try(CAT.K <- as.matrix(CAT.K[1:200,]))
try(CJS.K <- as.matrix(CJS.K[1:200,]))
try(CJT.K <- as.matrix(CJT.K[1:200,]))
try(STAS.K <- as.matrix(STAS.K[1:200,]))
try(STAT.K <- as.matrix(STAT.K[1:200,]))
try(STJS.K <- as.matrix(STJS.K[1:200,]))
try(STJT.K <- as.matrix(STJT.K[1:200,]))
try(WAS.K <- as.matrix(WAS.K[1:200,]))
try(WAT.K <- as.matrix(WAT.K[1:200,]))
try(WJS.K <- as.matrix(WJS.K[1:200,]))
try(WJT.K <- as.matrix(WJT.K[1:200,]))
# and arrange kernels as a list
ListK <- list(CAS.K, CAT.K, CJS.K, CJT.K, STAS.K, STAT.K, STJS.K, STJT.K, WAS.K, WAT.K, WJS.K, WJT.K)
# arrange a second list for corresponding tropical kernels (for the loop to estimate the mean matrices)
ListKtrop <- list(CAT.K, CAT.K, CJT.K, CJT.K, STAT.K, STAT.K, STJT.K, STJT.K, WAT.K, WAT.K, WJT.K, WJT.K)
# save space
rm(CAS.K, CAT.K, CJS.K, CJT.K, STAS.K, STAT.K, STJS.K, STJT.K, WAS.K, WAT.K, WJS.K, WJT.K)

# Estimate the sensitivities and elasticities of lambda for the 12 base models
CAS.L.sens <- sens(ListK[[1]], eval = "max"); CAS.L.elas <- elas(ListK[[1]], eval = "max")
CAT.L.sens <- sens(ListK[[2]], eval = "max"); CAT.L.elas <- elas(ListK[[2]], eval = "max")
CJS.L.sens <- sens(ListK[[3]], eval = "max"); CJS.L.elas <- elas(ListK[[3]], eval = "max")
CJT.L.sens <- sens(ListK[[4]], eval = "max"); CJT.L.elas <- elas(ListK[[4]], eval = "max")
STAS.L.sens <- sens(ListK[[5]], eval = "max"); STAS.L.elas <- elas(ListK[[5]], eval = "max")
STAT.L.sens <- sens(ListK[[6]], eval = "max"); STAT.L.elas <- elas(ListK[[6]], eval = "max")
STJS.L.sens <- sens(ListK[[7]], eval = "max"); STJS.L.elas <- elas(ListK[[7]], eval = "max")
STJT.L.sens <- sens(ListK[[8]], eval = "max"); STJT.L.elas <- elas(ListK[[8]], eval = "max")
WAS.L.sens <- sens(ListK[[9]], eval = "max"); WAS.L.elas <- elas(ListK[[9]], eval = "max")
WAT.L.sens <- sens(ListK[[10]], eval = "max"); WAT.L.elas <- elas(ListK[[10]], eval = "max")
WJS.L.sens <- sens(ListK[[11]], eval = "max"); WJS.L.elas <- elas(ListK[[11]], eval = "max")
WJT.L.sens <- sens(ListK[[12]], eval = "max"); WJT.L.elas <- elas(ListK[[12]], eval = "max")

# and repeat for the mean models
CA.L.sens <- sens((ListK[[1]]+ListK[[2]])/2, eval = "max"); CA.L.elas <- elas((ListK[[1]]+ListK[[2]])/2, eval = "max") 
CJ.L.sens <- sens((ListK[[3]]+ListK[[4]])/2, eval = "max"); CJ.L.elas <- elas((ListK[[3]]+ListK[[4]])/2, eval = "max") 
STA.L.sens <- sens((ListK[[5]]+ListK[[6]])/2, eval = "max"); STA.L.elas <- elas((ListK[[5]]+ListK[[6]])/2, eval = "max") 
STJ.L.sens <- sens((ListK[[7]]+ListK[[8]])/2, eval = "max"); STJ.L.elas <- elas((ListK[[7]]+ListK[[8]])/2, eval = "max") 
WA.L.sens <- sens((ListK[[9]]+ListK[[10]])/2, eval = "max"); WA.L.elas <- elas((ListK[[9]]+ListK[[10]])/2, eval = "max") 
WJ.L.sens <- sens((ListK[[11]]+ListK[[12]])/2, eval = "max"); WJ.L.elas <- elas((ListK[[11]]+ListK[[12]])/2, eval = "max") 


############################################
# STEP 3: Prepare variables and storage files for LTRE decomposition
############################################
# Condense all sensitivity and elasticity matrices together
# Base models
L.SensList <- list(CAS.L.sens, CAT.L.sens, CJS.L.sens, CJT.L.sens, 
                   STAS.L.sens, STAT.L.sens, STJS.L.sens, STJT.L.sens, 
                   WAS.L.sens, WAT.L.sens, WJS.L.sens, WJT.L.sens)
L.ElasList <- list(CAS.L.elas, CAT.L.elas, CJS.L.elas, CJT.L.elas, 
                   STAS.L.elas, STAT.L.elas, STJS.L.elas, STJT.L.elas, 
                   WAS.L.elas, WAT.L.elas, WJS.L.elas, WJT.L.elas)
MA.SensList <- list(CAS.MA.sens, CAT.MA.sens, CJS.MA.sens, CJT.MA.sens, 
                   STAS.MA.sens, STAT.MA.sens, STJS.MA.sens, STJT.MA.sens, 
                   WAS.MA.sens, WAT.MA.sens, WJS.MA.sens, WJT.MA.sens)
MA.ElasList <- list(CAS.MA.elas, CAT.MA.elas, CJS.MA.elas, CJT.MA.elas, 
                   STAS.MA.elas, STAT.MA.elas, STJS.MA.elas, STJT.MA.elas, 
                   WAS.MA.elas, WAT.MA.elas, WJS.MA.elas, WJT.MA.elas)
# Mean matrices
L.meanSensList <- list(CA.L.sens, CJ.L.sens, 
                       STA.L.sens, STJ.L.sens, 
                       WA.L.sens, WJ.L.sens)
L.meanElasList <- list(CA.L.elas, CJ.L.elas, 
                       STA.L.elas, STJ.L.elas, 
                       WA.L.elas, WJ.L.elas)
MA.meanSensList <- list(CA.MA.sens, CJ.MA.sens, 
                       STA.MA.sens, STJ.MA.sens, 
                       WA.MA.sens, WJ.MA.sens)
MA.meanElasList <- list(CA.MA.elas, CJ.MA.elas, 
                       STA.MA.elas, STJ.MA.elas, 
                       WA.MA.elas, WJ.MA.elas)

# and tidy up memory
rm(CAS.L.sens, CAT.L.sens, CJS.L.sens, CJT.L.sens, STAS.L.sens, STAT.L.sens, STJS.L.sens, STJT.L.sens, 
   WAS.L.sens, WAT.L.sens, WJS.L.sens, WJT.L.sens, CAS.L.elas, CAT.L.elas, CJS.L.elas, CJT.L.elas, 
   STAS.L.elas, STAT.L.elas, STJS.L.elas, STJT.L.elas, WAS.L.elas, WAT.L.elas, WJS.L.elas, WJT.L.elas,
   CAS.MA.sens, CAT.MA.sens, CJS.MA.sens, CJT.MA.sens, STAS.MA.sens, STAT.MA.sens, STJS.MA.sens, STJT.MA.sens, 
   WAS.MA.sens, WAT.MA.sens, WJS.MA.sens, WJT.MA.sens, CAS.MA.elas, CAT.MA.elas, CJS.MA.elas, CJT.MA.elas, 
   STAS.MA.elas, STAT.MA.elas, STJS.MA.elas, STJT.MA.elas, WAS.MA.elas, WAT.MA.elas, WJS.MA.elas, WJT.MA.elas,
   CA.L.sens, CJ.L.sens, STA.L.sens, STJ.L.sens, WA.L.sens, WJ.L.sens, 
   CA.L.elas, CJ.L.elas, STA.L.elas, STJ.L.elas, WA.L.elas, WJ.L.elas,
   CA.MA.sens, CJ.MA.sens, STA.MA.sens, STJ.MA.sens, WA.MA.sens, WJ.MA.sens,
   CA.MA.elas, CJ.MA.elas, STA.MA.elas, STJ.MA.elas, WA.MA.elas, WJ.MA.elas)

# read in lambda and maximal amplification estimates previous obtained for each model variant
C.outputs <- read.csv("C:/Users/rhian/OneDrive/Desktop/James/Chapter 4/Data/Competitive outputs.csv")
ST.outputs <- read.csv("C:/Users/rhian/OneDrive/Desktop/James/Chapter 4/Data/Stress Tolerant outputs.csv")
W.outputs <- read.csv("C:/Users/rhian/OneDrive/Desktop/James/Chapter 4/Data/Weedy outputs.csv")
# arrange together for the loop below
L_file <- do.call(cbind, list(C.outputs[,1:4], ST.outputs[,1:4], W.outputs[,1:4]))
MA_file <- do.call(cbind, list(C.outputs[,29:32], ST.outputs[,29:32], W.outputs[,29:32]))
colnames(L_file) <- colnames(MA_file) <- c("CAS", "CAT", "CJS", "CJT", 
                                           "STAS", "STAT", "STJS", "STJT", 
                                           "WAS", "WAT", "WJS", "WJT")
# and tidy up
rm(C.outputs, ST.outputs, W.outputs)

# and define max and min model boundaries (for this I need to load and format the original datafile)
# read in the full demographic dataset
demographic.data <- read.csv("C:/Users/rhian/OneDrive/Desktop/James/Chapter 4/Data/Demographic data.csv", stringsAsFactors = TRUE)
# subset demographic datafile to focus on only the relevant populations.
demographic.data <- demographic.data[which(demographic.data$Group == "Scleractinia" & 
                                             demographic.data$LHS != "UA"),]
# remove data from April 2016, all other surveys are for annual intervals.
demographic.data <- demographic.data[which(demographic.data$time.t %in% c("2016", "2017", "2018")),]
# Format variables
# remove entries corresponding with lost corals. These can't be used in estimating demographic parameters.
demographic.data <- subset(demographic.data, Size.t != "Lost" | is.na(Size.t)) #<- this removes all rows containing "Lost" in the size t section, but importantly retains any NA entries - these are really important!
demographic.data <- subset(demographic.data, Size.t1 != "Lost" | is.na(Size.t1)) #<- this does the same for Size at t+1
# convert variables from the dataframe into their required formats
demographic.data$Size.t <- as.numeric(paste(demographic.data$Size.t))
demographic.data$Size.t1 <- as.numeric(paste(demographic.data$Size.t1))
demographic.data$Tag.year <- as.factor(demographic.data$Tag.year)
demographic.data$No.frag.t1 <- as.numeric(paste(demographic.data$No.frag.t1))
demographic.data$LHS <- droplevels(demographic.data$LHS, exclude = "UA") # just to remove the Un-assigned category which has been dropped from the data.
# covert size data into cm^2 and transform 
demographic.data$Size.t <- demographic.data$Size.t/100; demographic.data$Size.t1 <- demographic.data$Size.t1/100 
demographic.data$Size.t <- log(demographic.data$Size.t)
demographic.data$Size.t1 <- log(demographic.data$Size.t1) 
# remove generalist corals 
demographic.data <- demographic.data[which(demographic.data$LHS %in% c("C", "ST", "W")),]
# The data is now ready for analyses

# Competitive
C.L <- 1.1 * min(demographic.data[which(demographic.data$LHS == "C"), c("Size.t","Size.t1")], na.rm = T) #minimum size - on the log scale small sizes are negative and so dividing would result in lower limits inside the actual size range.
C.U = 1.1 * max(demographic.data[which(demographic.data$LHS == "C"), c("Size.t","Size.t1")], na.rm = T) #maximum size
# because we have a couple of census periods the maximum and minimum sizes observed at any point during those censuses will be used for the limits of every projection
# (because since they have been observed they are feasible in any projection of this population)
# Stress-Tolerant
ST.L = 1.1 * min(demographic.data[which(demographic.data$LHS == "ST"), c("Size.t","Size.t1")], na.rm = T) 
ST.U = 1.1 * max(demographic.data[which(demographic.data$LHS == "ST"), c("Size.t","Size.t1")], na.rm = T) 
# Weedy
W.L = 1.1 * min(demographic.data[which(demographic.data$LHS == "W"), c("Size.t","Size.t1")], na.rm = T) 
W.U = 1.1 * max(demographic.data[which(demographic.data$LHS == "W"), c("Size.t","Size.t1")], na.rm = T) 
# And store the parameters together
L.store <- rep(c(C.L, ST.L, W.L), each = 2)
U.store <- rep(c(C.U, ST.U, W.U), each = 2)

# and remove datafile
rm(demographic.data)

# generate storage outputs
sens_outputs <- list()

# define indexing variables
LHS_index <- rep(c("C","ST","W"), each = 2)
Country_index <- rep(c("A","J"), length.out = 6)
Region_index <- rep(c("S","T"), each = 6)
file_select_start <- seq(1, 981, by = 20)
file_select_end <- seq(20, 1000, by = 20)
row_index_start <- seq(1, 4000, by = 200)
row_index_end <- seq(200, 4000, by = 200)

############################################
# STEP 4: LTRE decomposition
############################################
# start progress count
count <- 0

for (ii in 1:6) {
  row_count <- 0
  for (xx in 1:50) {
    # start by reading in the model variant parameters
    try(params_subtrop <- read.csv(paste0(LHS_index[ii], Country_index[ii], Region_index[ii]," par_", file_select_start[xx],"-", file_select_end[xx],".csv", sep = "")))
    try(params_trop <- read.csv(paste0(LHS_index[ii], Country_index[ii], Region_index[ii+6]," par_", file_select_start[xx],"-", file_select_end[xx],".csv", sep = "")))
    # and model variant F kernels (so that reproduction remains consistent across variant indexing)
    try(F_subtrop <- read.csv(paste0(LHS_index[ii], Country_index[ii], Region_index[ii]," F kernel_", file_select_start[xx],"-", file_select_end[xx],".csv", sep = "")))
    try(F_trop <- read.csv(paste0(LHS_index[ii], Country_index[ii], Region_index[ii+6]," F kernel_", file_select_start[xx],"-", file_select_end[xx],".csv", sep = "")))
    
    # repeat loop for each available parameter set
    for (jj in 1:20) {
      # advance indexing place holders
      count <- count+1
      row_count <- row_count+1
      # progress bar
      print(count)
      
      # Extract parameters
      try(subtrop_par <- as.numeric(params_subtrop[jj,]))
      try(trop_par <- as.numeric(params_trop[jj,]))
      try(names(trop_par) <- names(subtrop_par) <- names(params_trop))
      # Extract associated F matrix
      try(matF_trop <- as.matrix(F_trop[row_index_start[jj]:row_index_end[jj],]))
      try(matF_subtrop <- as.matrix(F_subtrop[row_index_start[jj]:row_index_end[jj],]))
      
      # Construct the models (Here I am only interested in recreating the survival, growth and fragmentation functions)
      try(trop_mod <- mk_P(m = ndim, m.par = trop_par, L = L.store[ii], U = U.store[ii]))
      try(subtrop_mod <- mk_P(m = ndim, m.par = subtrop_par, L = L.store[ii], U = U.store[ii]))
  
      # Determine vital rate level differences between the models
      
      ###### 1. Reproduction
      # Reproduction is easiest here as it is self contained independently within separate matrices
      # the subtropical model is being taken away from the tropical control (thus positive contributions suggest increases in the subtropics)
      try(ReprDiff <- matF_trop - matF_subtrop)
      
      ###### 2. Survival
      # Firstly Create a stage-based survival vector for each model
      try(sigma_trop <- colSums(trop_mod$P))
      try(sigma_subtrop <- colSums(subtrop_mod$P))

      # Now calculate a survival independent kernel for each model
      Surv_trop <- Surv_subtrop <- matrix(NA, ndim, ndim)
      # extract survival component from each element
      for(i in 1:ndim){
        try(Surv_trop[,i] <- trop_mod$P[,i]/sigma_trop[i])
        try(Surv_subtrop[,i] <- subtrop_mod$P[,i]/sigma_subtrop[i])
      }
      # and determine the element level differences
      try(SurvDiff <- Surv_trop - Surv_subtrop)
    
      ###### 3. Fragmentation
      # Fragmentation requires decomposition of the H kernel to remove survival
      Frag_trop <- Frag_subtrop <- matrix(NA, ndim, ndim)
      # extract survival component from each element
      for(i in 1:ndim){
        try(Frag_trop[,i] <- trop_mod$H[,i]/sigma_trop[i])
        try(Frag_subtrop[,i] <- subtrop_mod$H[,i]/sigma_subtrop[i])
      }
      # and determine the element level differences
      try(FragDiff <- Frag_trop - Frag_subtrop)
      
      # I also need to extract a stage-based fragmentation vector for each model to use in determining growth
      try(Frag_trop_vec <- colSums(Frag_trop))
      try(Frag_subtrop_vec <- colSums(Frag_subtrop))
      
      ###### 4. Growth
      # Now I can process growth (or more appropriately how do size transitions of colonies differ between the populations)
      # again remembering to remove survival.
      # First create a blank matrix for each model ready to store the elements following the removal of survival and fragmentation
      Growth_trop <- Growth_subtrop <- matrix(NA, ndim, ndim)
      # generate the growth independent models
      for(i in 1:ndim){
        try(Growth_trop[,i] <- trop_mod$G[,i]/(sigma_trop[i]*(1-Frag_trop_vec[i])))
        try(Growth_subtrop[,i] <- subtrop_mod$G[,i]/(sigma_subtrop[i]*(1-Frag_subtrop_vec[i])))
        # dividing each element of the original G kernels by the associated stage survival and fragmentation returns the element size transition value
      }
      # and determine the element level differences
      try(GrowDiff <- Growth_trop - Growth_subtrop)
      
      ###### and calculate LTREs by multiplying differences by mean matrix sensitivities
      # condense outputs into a storage list
      try(output <- list(LHS = LHS_index[ii],
                     Country = Country_index[ii],
                     Lambda_trop = L_file[row_count, paste0(LHS_index[ii], Country_index[ii], Region_index[ii+6])],
                     Lambda_subtrop = L_file[row_count, paste0(LHS_index[ii], Country_index[ii], Region_index[ii])],
                     MaxAmp_trop = MA_file[row_count, paste0(LHS_index[ii], Country_index[ii], Region_index[ii+6])],
                     MaxAmp_trop = MA_file[row_count, paste0(LHS_index[ii], Country_index[ii], Region_index[ii])],
                    
                     # Maximal Amplification
                     MA.sens.Surv = colSums(matrix(SurvDiff * MA.meanSensList[[ii]], ndim, ndim)),
                     MA.sens.Frag = colSums(matrix(FragDiff * MA.meanSensList[[ii]], ndim, ndim)),
                     MA.sens.Grow = colSums(matrix(GrowDiff * MA.meanSensList[[ii]], ndim, ndim)),
                     MA.sens.Repr = colSums(matrix(ReprDiff * MA.meanSensList[[ii]], ndim, ndim)),
                     MA.elas.Surv = colSums(matrix(SurvDiff * MA.meanElasList[[ii]], ndim, ndim)),
                     MA.elas.Frag = colSums(matrix(FragDiff * MA.meanElasList[[ii]], ndim, ndim)),
                     MA.elas.Grow = colSums(matrix(GrowDiff * MA.meanElasList[[ii]], ndim, ndim)),
                     MA.elas.Repr = colSums(matrix(ReprDiff * MA.meanElasList[[ii]], ndim, ndim)),
                     # Lambda
                     L.sens.Surv = colSums(matrix(SurvDiff * L.meanSensList[[ii]], ndim, ndim)),
                     L.sens.Frag = colSums(matrix(FragDiff * L.meanSensList[[ii]], ndim, ndim)),
                     L.sens.Grow = colSums(matrix(GrowDiff * L.meanSensList[[ii]], ndim, ndim)),
                     L.sens.Repr = colSums(matrix(ReprDiff * L.meanSensList[[ii]], ndim, ndim)),
                     L.elas.Surv = colSums(matrix(SurvDiff * L.meanElasList[[ii]], ndim, ndim)),
                     L.elas.Frag = colSums(matrix(FragDiff * L.meanElasList[[ii]], ndim, ndim)),
                     L.elas.Grow = colSums(matrix(GrowDiff * L.meanElasList[[ii]], ndim, ndim)),
                     L.elas.Repr = colSums(matrix(ReprDiff * L.meanElasList[[ii]], ndim, ndim))))
      
      # and store for extraction
      try(sens_outputs[[count]] <- output)
      
      # and tidy up
      try(rm(subtrop_par, trop_par, matF_subtrop, matF_trop, trop_mod, subtrop_mod, ReprDiff, SurvDiff, GrowDiff, FragDiff, output, 
         Growth_subtrop, Growth_trop, Frag_subtrop, Frag_trop, Surv_subtrop, Surv_trop, Frag_subtrop_vec, Frag_trop_vec, sigma_subtrop, sigma_trop))
    }
    
    # further tidying up
    try(rm(params_subtrop, params_trop, F_subtrop, F_trop))
  }
}

############################################
# STEP 5: Extract elements from the list for plotting
############################################

# First extract vectors of LHS and country to attach to extracted dataframes
LHS_vec <- unlist(lapply(sens_outputs, '[[', 1))
Country_vec <- unlist(lapply(sens_outputs, '[[', 2))

# Now extract the vital rate contributions (for now I'm only interested in Maximal Amplification)
SurvSensLTRE <- matrix(unlist(lapply(sens_outputs, '[[', 7)), ncol = 200, byrow = TRUE)
FragSensLTRE <- matrix(unlist(lapply(sens_outputs, '[[', 8)), ncol = 200, byrow = TRUE)
GrowSensLTRE <- matrix(unlist(lapply(sens_outputs, '[[', 9)), ncol = 200, byrow = TRUE)
ReprSensLTRE <- matrix(unlist(lapply(sens_outputs, '[[', 10)), ncol = 200, byrow = TRUE)
SurvElasLTRE <- matrix(unlist(lapply(sens_outputs, '[[', 11)), ncol = 200, byrow = TRUE)
FragElasLTRE <- matrix(unlist(lapply(sens_outputs, '[[', 12)), ncol = 200, byrow = TRUE)
GrowElasLTRE <- matrix(unlist(lapply(sens_outputs, '[[', 13)), ncol = 200, byrow = TRUE)
ReprElasLTRE <- matrix(unlist(lapply(sens_outputs, '[[', 14)), ncol = 200, byrow = TRUE)
# determine row sums for each run through
SurvSensSums <- rowSums(SurvSensLTRE)
FragSensSums <- rowSums(FragSensLTRE)
GrowSensSums <- rowSums(GrowSensLTRE)
ReprSensSums <- rowSums(ReprSensLTRE)
SurvElasSums <- rowSums(SurvElasLTRE)
FragElasSums <- rowSums(FragElasLTRE)
GrowElasSums <- rowSums(GrowElasLTRE)
ReprElasSums <- rowSums(ReprElasLTRE)

# And attach the LHS and country vectors and rowSums
LTRE.contributions <- data.frame(LHS = LHS_vec,
                                 Country = Country_vec,
                                 SurvSens = SurvSensSums,
                                 FragSens = FragSensSums,
                                 GrowSens = GrowSensSums,
                                 ReprSens = ReprSensSums,
                                 SurvElas = SurvElasSums,
                                 FragElas = FragElasSums,
                                 GrowElas = GrowElasSums,
                                 ReprElas = ReprElasSums)
# and tidy up data variation
LTRE.contributions[sapply(LTRE.contributions, is.nan)] <- NA
LTRE.contributions[sapply(LTRE.contributions, is.infinite)] <- NA

# estimate the proportional contributions
LTRE.contributions$SensSum <- NA
LTRE.contributions$ElasSum <- NA
#calculate the total contributions calculated using sensitivities and elasticities
for(xx in 1:dim(LTRE.contributions)[1]){
  LTRE.contributions[xx, "SensSum"] <- sum(sapply(LTRE.contributions[xx, 3:6], abs))
  LTRE.contributions[xx, "ElasSum"] <- sum(sapply(LTRE.contributions[xx, 7:10],abs))
}
# and now estimate the proportion contributed by each vital rate
LTRE.contributions$SurvSprop <- NA
LTRE.contributions$FragSprop <- NA
LTRE.contributions$GrowSprop <- NA
LTRE.contributions$ReprSprop <- NA
LTRE.contributions$SurvEprop <- NA
LTRE.contributions$FragEprop <- NA
LTRE.contributions$GrowEprop <- NA
LTRE.contributions$ReprEprop <- NA
# run a loop to calculate the proportions
for(xx in 1:dim(LTRE.contributions)[1]){
  LTRE.contributions[xx, "SurvSprop"] <- LTRE.contributions[xx, "SurvSens"]/LTRE.contributions[xx, "SensSum"]
  LTRE.contributions[xx, "FragSprop"] <- LTRE.contributions[xx, "FragSens"]/LTRE.contributions[xx, "SensSum"]
  LTRE.contributions[xx, "GrowSprop"] <- LTRE.contributions[xx, "GrowSens"]/LTRE.contributions[xx, "SensSum"]
  LTRE.contributions[xx, "ReprSprop"] <- LTRE.contributions[xx, "ReprSens"]/LTRE.contributions[xx, "SensSum"]
  LTRE.contributions[xx, "SurvEprop"] <- LTRE.contributions[xx, "SurvElas"]/LTRE.contributions[xx, "ElasSum"]
  LTRE.contributions[xx, "FragEprop"] <- LTRE.contributions[xx, "FragElas"]/LTRE.contributions[xx, "ElasSum"]
  LTRE.contributions[xx, "GrowEprop"] <- LTRE.contributions[xx, "GrowElas"]/LTRE.contributions[xx, "ElasSum"]
  LTRE.contributions[xx, "ReprEprop"] <- LTRE.contributions[xx, "ReprElas"]/LTRE.contributions[xx, "ElasSum"]
}

# Now sort through each coral group and country combination to determine the mean and variance for each vital rate variable
# create output
LTRE.output <- data.frame(LHS = rep(NA, length.out = 8*6),
                          Country = rep(NA, length.out = 8*6),
                          Vital_rate = rep(NA, length.out = 8*6), 
                          CI_Lower.est = rep(NA, length.out = 8*6),
                          Mean.est = rep(NA, length.out = 8*6),
                          CI_Upper.est = rep(NA, length.out = 8*6),
                          CI_Lower.prop = rep(NA, length.out = 8*6),
                          Mean.prop = rep(NA, length.out = 8*6),
                          CI_Upper.prop = rep(NA, length.out = 8*6))

# define output indexing
output.index.start <- seq(1,48,by = 8)
output.index.end <- seq(8,48,by = 8)

# run loop
for(ii in 1:6) {
  # subset the data according to the desired LHS and country combination
  sample.select <- LTRE.contributions[which(LTRE.contributions$LHS == LHS_index[ii] & LTRE.contributions$Country == Country_index[ii]),]
  
  # define output storage location
  index <- c(output.index.start[ii]:output.index.end[ii])
  
  # and calculate mean and CI for the LTRE contributions
  LTRE.output[index, "LHS"] <- LHS_index[ii]
  LTRE.output[index, "Country"] <- Country_index[ii]
  LTRE.output[index, "Vital_rate"] <- c("SurvS", "FragS", "GrowS", "ReprS","SurvE", "FragE", "GrowE", "ReprE")
  # Actual estimates
  LTRE.output[index[1],"CI_Lower.est"] <- CI(na.omit(sample.select$SurvSens))[3]; LTRE.output[index[1],"Mean.est"] <- CI(na.omit(sample.select$SurvSens))[2]; LTRE.output[index[1],"CI_Upper.est"] <- CI(na.omit(sample.select$SurvSens))[1]
  LTRE.output[index[2],"CI_Lower.est"] <- CI(na.omit(sample.select$FragSens))[3]; LTRE.output[index[2],"Mean.est"] <- CI(na.omit(sample.select$FragSens))[2]; LTRE.output[index[2],"CI_Upper.est"] <- CI(na.omit(sample.select$FragSens))[1]
  LTRE.output[index[3],"CI_Lower.est"] <- CI(na.omit(sample.select$GrowSens))[3]; LTRE.output[index[3],"Mean.est"] <- CI(na.omit(sample.select$GrowSens))[2]; LTRE.output[index[3],"CI_Upper.est"] <- CI(na.omit(sample.select$GrowSens))[1]
  LTRE.output[index[4],"CI_Lower.est"] <- CI(na.omit(sample.select$ReprSens))[3]; LTRE.output[index[4],"Mean.est"] <- CI(na.omit(sample.select$ReprSens))[2]; LTRE.output[index[4],"CI_Upper.est"] <- CI(na.omit(sample.select$ReprSens))[1]
  LTRE.output[index[5],"CI_Lower.est"] <- CI(na.omit(sample.select$SurvElas))[3]; LTRE.output[index[5],"Mean.est"] <- CI(na.omit(sample.select$SurvElas))[2]; LTRE.output[index[5],"CI_Upper.est"] <- CI(na.omit(sample.select$SurvElas))[1]
  LTRE.output[index[6],"CI_Lower.est"] <- CI(na.omit(sample.select$FragElas))[3]; LTRE.output[index[6],"Mean.est"] <- CI(na.omit(sample.select$FragElas))[2]; LTRE.output[index[6],"CI_Upper.est"] <- CI(na.omit(sample.select$FragElas))[1]
  LTRE.output[index[7],"CI_Lower.est"] <- CI(na.omit(sample.select$GrowElas))[3]; LTRE.output[index[7],"Mean.est"] <- CI(na.omit(sample.select$GrowElas))[2]; LTRE.output[index[7],"CI_Upper.est"] <- CI(na.omit(sample.select$GrowElas))[1]
  LTRE.output[index[8],"CI_Lower.est"] <- CI(na.omit(sample.select$ReprElas))[3]; LTRE.output[index[8],"Mean.est"] <- CI(na.omit(sample.select$ReprElas))[2]; LTRE.output[index[8],"CI_Upper.est"] <- CI(na.omit(sample.select$ReprElas))[1]
  # proportional estimates
  LTRE.output[index[1],"CI_Lower.prop"] <- CI(na.omit(sample.select$SurvSprop))[3]; LTRE.output[index[1],"Mean.prop"] <- CI(na.omit(sample.select$SurvSprop))[2]; LTRE.output[index[1],"CI_Upper.prop"] <- CI(na.omit(sample.select$SurvSprop))[1]
  LTRE.output[index[2],"CI_Lower.prop"] <- CI(na.omit(sample.select$FragSprop))[3]; LTRE.output[index[2],"Mean.prop"] <- CI(na.omit(sample.select$FragSprop))[2]; LTRE.output[index[2],"CI_Upper.prop"] <- CI(na.omit(sample.select$FragSprop))[1]
  LTRE.output[index[3],"CI_Lower.prop"] <- CI(na.omit(sample.select$GrowSprop))[3]; LTRE.output[index[3],"Mean.prop"] <- CI(na.omit(sample.select$GrowSprop))[2]; LTRE.output[index[3],"CI_Upper.prop"] <- CI(na.omit(sample.select$GrowSprop))[1]
  LTRE.output[index[4],"CI_Lower.prop"] <- CI(na.omit(sample.select$ReprSprop))[3]; LTRE.output[index[4],"Mean.prop"] <- CI(na.omit(sample.select$ReprSprop))[2]; LTRE.output[index[4],"CI_Upper.prop"] <- CI(na.omit(sample.select$ReprSprop))[1]
  LTRE.output[index[5],"CI_Lower.prop"] <- CI(na.omit(sample.select$SurvEprop))[3]; LTRE.output[index[5],"Mean.prop"] <- CI(na.omit(sample.select$SurvEprop))[2]; LTRE.output[index[5],"CI_Upper.prop"] <- CI(na.omit(sample.select$SurvEprop))[1]
  LTRE.output[index[6],"CI_Lower.prop"] <- CI(na.omit(sample.select$FragEprop))[3]; LTRE.output[index[6],"Mean.prop"] <- CI(na.omit(sample.select$FragEprop))[2]; LTRE.output[index[6],"CI_Upper.prop"] <- CI(na.omit(sample.select$FragEprop))[1]
  LTRE.output[index[7],"CI_Lower.prop"] <- CI(na.omit(sample.select$GrowEprop))[3]; LTRE.output[index[7],"Mean.prop"] <- CI(na.omit(sample.select$GrowEprop))[2]; LTRE.output[index[7],"CI_Upper.prop"] <- CI(na.omit(sample.select$GrowEprop))[1]
  LTRE.output[index[8],"CI_Lower.prop"] <- CI(na.omit(sample.select$ReprEprop))[3]; LTRE.output[index[8],"Mean.prop"] <- CI(na.omit(sample.select$ReprEprop))[2]; LTRE.output[index[8],"CI_Upper.prop"] <- CI(na.omit(sample.select$ReprEprop))[1]
  
  # and tidy up
  rm(sample.select)
}

# Now split the data for easy plotting
LTRE.outputSens <- LTRE.output[which(LTRE.output$Vital_rate %in% c("SurvS", "FragS", "GrowS", "ReprS")),]
LTRE.outputElas <- LTRE.output[which(LTRE.output$Vital_rate %in% c("SurvE", "FragE", "GrowE", "ReprE")),]

# and plot the data 
rate_order <- c('SurvE', 'GrowE', 'FragE', 'ReprE')
# create plot
LTRE.contributions.plots <- ggplot(LTRE.outputElas, aes(x = Vital_rate)) +
  geom_bar(aes(y = Mean.prop, fill = LHS, colour = LHS), stat = "identity", alpha = 0.9) +
  geom_errorbar(aes(ymin = CI_Lower.prop, ymax = CI_Upper.prop), width = .2, size = 0.7) +
  facet_grid(LHS ~ Country) +
  # Formatting
  scale_x_discrete(limits = rate_order, 
                   labels=c("??","??","??","??")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(-0.6,0.6),
                     breaks = c(-0.5,0,0.5), labels = c(-0.5,0,0.5)) +
  scale_fill_manual(name = "", values = c("blue2","orange3","red2"), 
                    labels = c("Competitive    ","Stress-Tolerant","Weedy          "), 
                    guide = FALSE) +
  scale_color_manual(name = "", values = c("blue2","orange3","red2"), 
                     labels = c("Competitive    ","Stress-Tolerant","Weedy          "), 
                     guide = FALSE) +
  xlab("") + 
  ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(size = 15, colour = "black", face = "italic"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0.9),axis.line.y.right = element_line(colour = "black", size = 1),
        legend.background = element_blank(), legend.key = element_blank(),legend.position = "bottom",
        legend.box = "horizontal", panel.spacing.y = unit(1.5, "lines"))

################################################################### End of Code #######################
