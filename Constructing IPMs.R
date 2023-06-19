# This script is for defining the demographic functions that will be used to construct integral projection subkernels 
# representing the survival, growth, fragmentation and recruitment dynamics of tropical and subtropical coral populations in both Australia and Japan.
# These functions will then be used to construct IPMs for estimating values of lambda and Kriess Bounds for each target coral assemblage and their variance (using Jack-knife resampling).
# This script requires the vital rate regression coefficients extracted using the 'Quantifying vital rate patterns' script. 

# Author: James Cant (jic2@st-andrews,ac.uk)
# Date last modified: April 2020
# -----------------------------------------------------------------------------------------

# set working directory
setwd("Data_file_location")

# Load required packages
library(fields)
library(popdemo)
library(plyr)
library(stats)
library(tidyr)
library(lme4)

# Load required data
m.par.store <- read.csv("vital rate parameters.csv", row.names = 1)
demographic.data <- read.csv("Formatted data.csv", stringsAsFactors = TRUE)
fecundity.data <- read.csv("Fecundity data.csv", stringsAsFactors = TRUE)

# set random number seed 
set.seed(43400)
# define indexing for the loop below
a = 2
b = 1000

#######################################
# STEP 1: Define IPM functions
#######################################

# 1. Define function for discretising and combining the survival + growth (G) fragmentation (H) (G + H = P) and recruitment (F) kernels
mk_K <- function(m, m.par, L, U) {
  # define the mesh points 
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  # determine the number of new larvae settling 
  recruit.no <- recruits(m.par)
  # Define the survival and growth kernel
  G <- h * (outer(meshpts, meshpts, G_z1z, m.par))
  # Define the fragmentation kernel
  H <- h * (outer(meshpts, meshpts, H_z1z, m.par))
  # Define the recruitment kernel
  F <- h * (outer(meshpts, meshpts, F_z1z, m.par, recruit.no))
  # define adult only kernel (P and H)
  P <- G + H 
  # combine all sub kernels together to build the full model.
  K <- P + F 
  return(list(K = K, P = P, meshpts = meshpts, G = G, H = H, F = F, R = recruit.no))
} # The model is also returning the number of recruits added so that I can 
# store them and have them associated with corresponding lambda/Kreiss values.

# 2. Define the survival/growth kernel 
# For this kernel I am only interested in the survival and growth of colonies that DONT fragment.
G_z1z <- function (Size.t1, Size.t, m.par) {
  return( (1 - h_z(Size.t, m.par)) * s_z(Size.t, m.par) * g_z1z(Size.t1, Size.t, m.par) ) 
}        # Colonies don't fragment      # survive            # and grow/shrink

# 3. Define the dynamics of colonies that DO fragment (a measure of clonality/asexual reproduction)
H_z1z <- function (Size.t1, Size.t, m.par) {
  return( s_z(Size.t, m.par) * h_z(Size.t, m.par) * pf_z(Size.t, m.par) * f_z1(Size.t1, Size.t, m.par) )
}        # Colonies survive     # and fragment      # producing x number of fragments of a given size

# 4. Define the production of new recruits into the populations.
F_z1z <- function (Size.t1, Size.t, m.par, recruit.no) {
  R <- settleR(Size.t, m.par, recruit.no) # estimate the proportional fecundity of colonies
  return( F <- R * c0_z1(Size.t1, m.par) )
}     # Number of recruits # Size of successful recruits


#######################################
# STEP 2: Define vital rate functions
#######################################

# 1. Colony survival
s_z = function(Size.t, m.par) {
  linear_p <- m.par["surv.int"] + (m.par["surv.slope"] * Size.t) 
  p <- 1/(1+exp(-linear_p)) # this accounts for the logistic nature of the survival regression model                          
  return(p)
}

# 2. Colony growth
g_z1z = function(Size.t1, Size.t, m.par) {
  mean <- m.par["grow.int"] + (m.par["grow.slope"] * Size.t) + (m.par["grow.slope2"] * (Size.t^2))
  sd <- sqrt(pi/2) * exp((m.par["grow.sd.int"] + m.par["grow.sd.slope"] * Size.t)) # this allows for the variance in new size to change with initial size.
  p_den_grow <- dnorm(Size.t1, mean = mean, sd = sd)
  return(p_den_grow)
} #linear growth

# 3. Colony Fragmentation
h_z <- function(Size.t, m.par){
  linear_h <- m.par["frag.int"] + (m.par["frag.slope"] * Size.t) + (m.par["frag.slope2"] * (Size.t^2))
  h <- 1/(1+exp(-linear_h)) # logistic regression
  return(h)
}

# 4. Number of fragments produced.
pf_z <- function(Size.t, m.par){
  p <- exp(m.par["frag.no.int"] + m.par["frag.no.slope"] * Size.t)
  return(p)
} 

# 5. Fragment size
f_z1 <- function(Size.t1, Size.t, m.par){
  mean <- m.par["frag.size.int"] + (m.par["frag.size.slope"] * Size.t) + (m.par["frag.size.slope2"] * (Size.t^2))
  sd <- sqrt(pi/2) * exp((m.par["frag.sd.int"] + m.par["frag.sd.slope"] * Size.t)) # this allows for the variance in fragment size to change with adult colony size.
  p_den_fsz <- dnorm(Size.t1, mean = mean, sd = sd) # calculate the probability of your size given any starting size. The Dnorm function creates a normal size distribution from this data - because remember IPMs create population density distributions.
  return(p_den_fsz)
}

# 6. A function that randomly selects a number of settling larvae (from within observed limits)
recruits <- function(m.par){
  R1 <- m.par["Rmax"]
  R2 <- m.par["Rmin"]
  rawR <- runif(1, min = R2, max = R1) # generate a random estimate from with the estimated range
  roundR <- round(rawR) # ensures we are dealing with a whole number of recruits.
  return(roundR)
}

# 7.  This function converts size-specific fecundity into a proportional recruit contribution (this function completes the IPM loop between 
# the dynamics of existing colonies and observed recruitment patterns. 
settleR <- function(Size.t, m.par, recruit.no){
  # firstly estimate colony fecundity
  rawf <- m.par["fec.a"] * exp(m.par["fec.b"] * Size.t)
  # Adjust this larval output according to the randomized estimate of larval counts
  maxf <- sum(rawf)
  ratiof <- recruit.no/maxf # generates a ratio translating larval volume output (colony fecundity) into the observed number of settling larvae
  truef <- rawf * ratiof # these are now proportional recruit contributions given each colonies size.
  return(truef)
} 

# 8. Recruit Size
# the pdf of size z1 recruits next summer - 
c0_z1 <- function(Size.t1, m.par){
  mean <- m.par["rec.size"]
  sd <- m.par["rec.size.sd"]
  p_den_rcsz <- dnorm(Size.t1, mean = mean, sd = sd) 
  return(p_den_rcsz)
}

#######################################
# STEP 3: Estimate parameters defining the dimension of model kernels
#######################################

# Firstly how different are the size ranges of each population
# Competitive
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "C"),"Size.t1"])
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "C"),"Size.t1"])
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "C"),"Size.t1"])
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "C"),"Size.t1"])
# Stress -Tolerant
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "ST"),"Size.t1"])
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "ST"),"Size.t1"])
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "ST"),"Size.t1"])
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "ST"),"Size.t1"])
# Weedy
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "W"),"Size.t1"])
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "W"),"Size.t1"])
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "W"),"Size.t1"])
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "W"),"Size.t1"])
# Within LHS groups the size ranges are similar across regions

# Now to define the boundary sizes and the number of bins to be used in discretising the model kernels for each coral group.
# I will keep dimensions constant across countries and ecoregions but will allow them to differ as needed across LHS groups.
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
L.store <- c(C.L, ST.L, W.L)
U.store <- c(C.U, ST.U, W.U)

# Estimate the number of bins (m). Bin number will remain consistent within but not between LHS groups.
# This is done using a histogram of the relevant populations to tell us the appropriate (minimum) bin sizes.
# Competitive
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "C"),"Size.t1"], breaks = 162)
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "C"),"Size.t1"], breaks = 150)
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "C"),"Size.t1"], breaks = 186)
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "C"),"Size.t1"], breaks = 138)
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "C"),"Size.t"], breaks = 162)
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "C"),"Size.t"], breaks = 142)
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "C"),"Size.t"], breaks = 164)
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "C"),"Size.t"], breaks = 128)
# Stress -Tolerant
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "ST"),"Size.t1"], breaks = 154)
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "ST"),"Size.t1"], breaks = 175)
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "ST"),"Size.t1"], breaks = 161)
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "ST"),"Size.t1"], breaks = 168)
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "ST"),"Size.t1"], breaks = 154)
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "ST"),"Size.t1"], breaks = 175)
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "ST"),"Size.t1"], breaks = 161)
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "ST"),"Size.t1"], breaks = 168)
# Weedy
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "W"),"Size.t1"], breaks = 136)
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "W"),"Size.t1"], breaks = 122)
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "W"),"Size.t1"], breaks = 132)
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "W"),"Size.t1"], breaks = 122)#
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "W"),"Size.t1"], breaks = 136)
hist(demographic.data[which(demographic.data$Country == "Australia" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "W"),"Size.t1"], breaks = 122)
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Subtropical" & demographic.data$LHS == "W"),"Size.t1"], breaks = 132)
hist(demographic.data[which(demographic.data$Country == "Japan" & demographic.data$Ecoregion == "Tropical" & demographic.data$LHS == "W"),"Size.t1"], breaks = 122)

# Now define the number of bins 
# The number of bins used needs to be slightly higher than optimum to account for the extended size range for integration.
# And will be standardized across all models  
m <- 200 

#######################################
# STEP 4: Implement IPM models and estimate Lambda and Kriess bounds
#######################################
# Generate output storage.
C.outputs <- matrix(NA, nrow = 1000, ncol = 36)
ST.outputs <- matrix(NA, nrow = 1000, ncol = 36)
W.outputs <- matrix(NA, nrow = 1000, ncol = 36)
colnames(C.outputs) <- colnames(ST.outputs) <- colnames(W.outputs) <- c("AS Lambda", "AT Lambda", "JS Lambda", "JT Lambda",
                                                                        "AS U.Kreiss", "AT U.Kreiss", "JS U.Kreiss", "JT U.Kreiss",
                                                                        "AS L.Kreiss", "AT L.Kreiss", "JS L.Kreiss", "JT L.Kreiss",
                                                                        "AS TE", "AT TE", "JS TE", "JT TE",
                                                                        "R.AS", "R.AT", "R.JS", "R.JT",
                                                                        "AS reac", "AT reac", "JS reac", "JT reac",
                                                                        "AS att", "AT att", "JS att", "JT att",
                                                                        "AS maxamp", "AT maxamp", "JS maxamp", "JT maxamp",
                                                                        "AS maxatt", "AT maxatt", "JS maxatt", "JT maxatt") 
# Moment of truth

### 1. Construct each IPM and display lambda to check the models have worked
# Competitive
CAS.IPM <- mk_K(m = m, m.par = unlist(m.par.store[1,]), L = L.store[1], U = U.store[1])
CAS.lambdaK <- Re(eigen(CAS.IPM$K)$values)[1]; CAS.lambdaK; C.outputs[1,1] <- CAS.lambdaK

CAT.IPM <- mk_K(m = m, m.par = unlist(m.par.store[2,]), L = L.store[1], U = U.store[1])
CAT.lambdaK <- Re(eigen(CAT.IPM$K)$values)[1]; CAT.lambdaK; C.outputs[1,2] <- CAT.lambdaK

CJS.IPM <- mk_K(m = m, m.par = unlist(m.par.store[3,]), L = L.store[1], U = U.store[1])
CJS.lambdaK <- Re(eigen(CJS.IPM$K)$values)[1]; CJS.lambdaK; C.outputs[1,3] <- CJS.lambdaK

CJT.IPM <- mk_K(m = m, m.par = unlist(m.par.store[4,]), L = L.store[1], U = U.store[1])
CJT.lambdaK <- Re(eigen(CJT.IPM$K)$values)[1]; CJT.lambdaK; C.outputs[1,4] <- CJT.lambdaK

# Stress Tolerant
STAS.IPM <- mk_K(m = m, m.par = unlist(m.par.store[5,]), L = L.store[2], U = U.store[2])
STAS.lambdaK <- Re(eigen(STAS.IPM$K)$values)[1]; STAS.lambdaK; ST.outputs[1,1] <- STAS.lambdaK

STAT.IPM <- mk_K(m = m, m.par = unlist(m.par.store[6,]), L = L.store[2], U = U.store[2])
STAT.lambdaK <- Re(eigen(STAT.IPM$K)$values)[1]; STAT.lambdaK; ST.outputs[1,2] <- STAT.lambdaK

STJS.IPM <- mk_K(m = m, m.par = unlist(m.par.store[7,]), L = L.store[2], U = U.store[2])
STJS.lambdaK <- Re(eigen(STJS.IPM$K)$values)[1]; STJS.lambdaK; ST.outputs[1,3] <- STJS.lambdaK

STJT.IPM <- mk_K(m = m, m.par = unlist(m.par.store[8,]), L = L.store[2], U = U.store[2])
STJT.lambdaK <- Re(eigen(STJT.IPM$K)$values)[1]; STJT.lambdaK; ST.outputs[1,4] <- STJT.lambdaK

# Weedy
WAS.IPM <- mk_K(m = m, m.par = unlist(m.par.store[9,]), L = L.store[3], U = U.store[3])
WAS.lambdaK <- Re(eigen(WAS.IPM$K)$values)[1]; WAS.lambdaK; W.outputs[1,1] <- WAS.lambdaK

WAT.IPM <- mk_K(m = m, m.par = unlist(m.par.store[10,]), L = L.store[3], U = U.store[3])
WAT.lambdaK <- Re(eigen(WAT.IPM$K)$values)[1]; WAT.lambdaK; W.outputs[1,2] <- WAT.lambdaK

WJS.IPM <- mk_K(m = m, m.par = unlist(m.par.store[11,]), L = L.store[3], U = U.store[3])
WJS.lambdaK <- Re(eigen(WJS.IPM$K)$values)[1]; WJS.lambdaK; W.outputs[1,3] <- WJS.lambdaK

WJT.IPM <- mk_K(m = m, m.par = unlist(m.par.store[12,]), L = L.store[3], U = U.store[3])
WJT.lambdaK <- Re(eigen(WJT.IPM$K)$values)[1]; WJT.lambdaK; W.outputs[1,4] <- WJT.lambdaK

# Estimate the Kreiss bounds for each IPM.
# Competitive
C.outputs[1,5] <- Kreiss(CAS.IPM$K, bound = "upper"); C.outputs[1,9] <- Kreiss(CAS.IPM$K, bound = "lower")
C.outputs[1,6] <- Kreiss(CAT.IPM$K, bound = "upper"); C.outputs[1,10] <- Kreiss(CAT.IPM$K, bound = "lower")
C.outputs[1,7] <- Kreiss(CJS.IPM$K, bound = "upper"); C.outputs[1,11] <- Kreiss(CJS.IPM$K, bound = "lower")
C.outputs[1,8] <- Kreiss(CJT.IPM$K, bound = "upper"); C.outputs[1,12] <- Kreiss(CJT.IPM$K, bound = "lower")
# Stress-Tolerant
ST.outputs[1,5] <- Kreiss(STAS.IPM$K, bound = "upper"); ST.outputs[1,9] <- Kreiss(STAS.IPM$K, bound = "lower")
ST.outputs[1,6] <- Kreiss(STAT.IPM$K, bound = "upper"); ST.outputs[1,10] <- Kreiss(STAT.IPM$K, bound = "lower")
ST.outputs[1,7] <- Kreiss(STJS.IPM$K, bound = "upper"); ST.outputs[1,11] <- Kreiss(STJS.IPM$K, bound = "lower")
ST.outputs[1,8] <- Kreiss(STJT.IPM$K, bound = "upper"); ST.outputs[1,12] <- Kreiss(STJT.IPM$K, bound = "lower")
# Weedy
W.outputs[1,5] <- Kreiss(WAS.IPM$K, bound = "upper"); W.outputs[1,9] <- Kreiss(WAS.IPM$K, bound = "lower")
W.outputs[1,6] <- Kreiss(WAT.IPM$K, bound = "upper"); W.outputs[1,10] <- Kreiss(WAT.IPM$K, bound = "lower")
W.outputs[1,7] <- Kreiss(WJS.IPM$K, bound = "upper"); W.outputs[1,11] <- Kreiss(WJS.IPM$K, bound = "lower")
W.outputs[1,8] <- Kreiss(WJT.IPM$K, bound = "upper"); W.outputs[1,12] <- Kreiss(WJT.IPM$K, bound = "lower")

# determine the Transient envelopes of each population 
# Competitive
C.outputs[1,13] <- C.outputs[1,5] - C.outputs[1,9]
C.outputs[1,14] <- C.outputs[1,6] - C.outputs[1,10]
C.outputs[1,15] <- C.outputs[1,7] - C.outputs[1,11]
C.outputs[1,16] <- C.outputs[1,8] - C.outputs[1,12]
# Stress-Tolerant
ST.outputs[1,13] <- ST.outputs[1,5] - ST.outputs[1,9]
ST.outputs[1,14] <- ST.outputs[1,6] - ST.outputs[1,10]
ST.outputs[1,15] <- ST.outputs[1,7] - ST.outputs[1,11]
ST.outputs[1,16] <- ST.outputs[1,8] - ST.outputs[1,12]
# Weedy
W.outputs[1,13] <- W.outputs[1,5] - W.outputs[1,9]
W.outputs[1,14] <- W.outputs[1,6] - W.outputs[1,10]
W.outputs[1,15] <- W.outputs[1,7] - W.outputs[1,11]
W.outputs[1,16] <- W.outputs[1,8] - W.outputs[1,12]

# And add the associated recruitment value.
# Competitive
C.outputs[1,17] <- CAS.IPM$R
C.outputs[1,18] <- CAT.IPM$R
C.outputs[1,19] <- CJS.IPM$R
C.outputs[1,20] <- CJT.IPM$R
# Stress-Tolerant
ST.outputs[1,17] <- STAS.IPM$R
ST.outputs[1,18] <- STAT.IPM$R
ST.outputs[1,19] <- STJS.IPM$R
ST.outputs[1,20] <- STJT.IPM$R
# Weedy
W.outputs[1,17] <- WAS.IPM$R
W.outputs[1,18] <- WAT.IPM$R
W.outputs[1,19] <- WJS.IPM$R
W.outputs[1,20] <- WJT.IPM$R

# estimate reactivity
# Competitive
try(C.outputs[1,21] <- reac(CAS.IPM$K, bound = "upper"))
try(C.outputs[1,22] <- reac(CAT.IPM$K, bound = "upper"))
try(C.outputs[1,23] <- reac(CJS.IPM$K, bound = "upper"))
try(C.outputs[1,24] <- reac(CJT.IPM$K, bound = "upper"))
# Stress-Tolerant
try(ST.outputs[1,21] <- reac(STAS.IPM$K, bound = "upper"))
try(ST.outputs[1,22] <- reac(STAT.IPM$K, bound = "upper"))
try(ST.outputs[1,23] <- reac(STJS.IPM$K, bound = "upper"))
try(ST.outputs[1,24] <- reac(STJT.IPM$K, bound = "upper"))
# Weedy
try(W.outputs[1,21] <- reac(WAS.IPM$K, bound = "upper"))
try(W.outputs[1,22] <- reac(WAT.IPM$K, bound = "upper"))
try(W.outputs[1,23] <- reac(WJS.IPM$K, bound = "upper"))
try(W.outputs[1,24] <- reac(WJT.IPM$K, bound = "upper"))

# estimate attenuation
# Competitive
try(C.outputs[1,25] <- reac(CAS.IPM$K, bound = "lower"))
try(C.outputs[1,26] <- reac(CAT.IPM$K, bound = "lower"))
try(C.outputs[1,27] <- reac(CJS.IPM$K, bound = "lower"))
try(C.outputs[1,28] <- reac(CJT.IPM$K, bound = "lower"))
# Stress-Tolerant
try(ST.outputs[1,25] <- reac(STAS.IPM$K, bound = "lower"))
try(ST.outputs[1,26] <- reac(STAT.IPM$K, bound = "lower"))
try(ST.outputs[1,27] <- reac(STJS.IPM$K, bound = "lower"))
try(ST.outputs[1,28] <- reac(STJT.IPM$K, bound = "lower"))
# Weedy
try(W.outputs[1,25] <- reac(WAS.IPM$K, bound = "lower"))
try(W.outputs[1,26] <- reac(WAT.IPM$K, bound = "lower"))
try(W.outputs[1,27] <- reac(WJS.IPM$K, bound = "lower"))
try(W.outputs[1,28] <- reac(WJT.IPM$K, bound = "lower"))

# estimate maximal amplification
# Competitive
try(C.outputs[1,29] <- maxamp(CAS.IPM$K))
try(C.outputs[1,30] <- maxamp(CAT.IPM$K))
try(C.outputs[1,31] <- maxamp(CJS.IPM$K))
try(C.outputs[1,32] <- maxamp(CJT.IPM$K))
# Stress-Tolerant
try(ST.outputs[1,29] <- maxamp(STAS.IPM$K))
try(ST.outputs[1,30] <- maxamp(STAT.IPM$K))
try(ST.outputs[1,21] <- maxamp(STJS.IPM$K))
try(ST.outputs[1,32] <- maxamp(STJT.IPM$K))
# Weedy
try(W.outputs[1,29] <- maxamp(WAS.IPM$K))
try(W.outputs[1,30] <- maxamp(WAT.IPM$K))
try(W.outputs[1,31] <- maxamp(WJS.IPM$K))
try(W.outputs[1,32] <- maxamp(WJT.IPM$K))

# and finally estimate maximal attenuation
try(C.outputs[1,33] <- maxatt(CAS.IPM$K))
try(C.outputs[1,34] <- maxatt(CAT.IPM$K))
try(C.outputs[1,35] <- maxatt(CJS.IPM$K))
try(C.outputs[1,36] <- maxatt(CJT.IPM$K))
# Stress-Tolerant
try(ST.outputs[1,33] <- maxatt(STAS.IPM$K))
try(ST.outputs[1,34] <- maxatt(STAT.IPM$K))
try(ST.outputs[1,35] <- maxatt(STJS.IPM$K))
try(ST.outputs[1,36] <- maxatt(STJT.IPM$K))
# Weedy
try(W.outputs[1,33] <- maxatt(WAS.IPM$K))
try(W.outputs[1,34] <- maxatt(WAT.IPM$K))
try(W.outputs[1,35] <- maxatt(WJS.IPM$K))
try(W.outputs[1,36] <- maxatt(WJT.IPM$K))

#######################################
# STEP 5: Run a Jack-knife to estimate variance in these various estimates.
#######################################
# create a loop that will repeat the IPM construction numerous times with 95% of the data to determine the variability in lambda and Kriess bounds.

# create storage for model outputs that will be needed for assessing the variance in vital rate sensitivities
CAS.K <- CAT.K <- CJS.K <- CJT.K <- list() # for storing the K kernels created during each re-run.
CAS.P <- CAT.P <- CJS.P <- CJT.P <- list() # storage for the P kernels
CAS.F <- CAT.F <- CJS.F <- CJT.F <- list() # storage for the F kernels
STAS.K <- STAT.K <- STJS.K <- STJT.K <- list() 
STAS.P <- STAT.P <- STJS.P <- STJT.P <- list() 
STAS.F <- STAT.F <- STJS.F <- STJT.F <- list() 
WAS.K <- WAT.K <- WJS.K <- WJT.K <- list() 
WAS.P <- WAT.P <- WJS.P <- WJT.P <- list() 
WAS.F <- WAT.F <- WJS.F <- WJT.F <- list() 
# for storing meshpoints
C.meshpts <- ST.meshpts <- W.meshpts <- list()
# for storing model parameters.
CAS.m.par <- CAT.m.par <- CJS.m.par <- CJT.m.par <- list() 
STAS.m.par <- STAT.m.par <- STJS.m.par <- STJT.m.par <- list() 
WAS.m.par <- WAT.m.par <- WJS.m.par <- WJT.m.par <- list() 

# And store the outputs from the first run through.
# Kernels
CAS.K[[1]] <- CAS.IPM$K; CAT.K[[1]] <- CAT.IPM$K; CJS.K[[1]] <- CJS.IPM$K; CJT.K[[1]] <- CJT.IPM$K 
CAS.P[[1]] <- CAS.IPM$P; CAT.P[[1]] <- CAT.IPM$P; CJS.P[[1]] <- CJS.IPM$P; CJT.P[[1]] <- CJT.IPM$P 
CAS.F[[1]] <- CAS.IPM$F; CAT.F[[1]] <- CAT.IPM$F; CJS.F[[1]] <- CJS.IPM$F; CJT.F[[1]] <- CJT.IPM$F 
STAS.K[[1]] <- STAS.IPM$K; STAT.K[[1]] <- STAT.IPM$K; STJS.K[[1]] <- STJS.IPM$K; STJT.K[[1]] <- STJT.IPM$K  
STAS.P[[1]] <- STAS.IPM$P; STAT.P[[1]] <- STAT.IPM$P; STJS.P[[1]] <- STJS.IPM$P; STJT.P[[1]] <- STJT.IPM$P  
STAS.F[[1]] <- STAS.IPM$F; STAT.F[[1]] <- STAT.IPM$F; STJS.F[[1]] <- STJS.IPM$F; STJT.F[[1]] <- STJT.IPM$F  
WAS.K[[1]] <- WAS.IPM$K; WAT.K[[1]] <- WAT.IPM$K; WJS.K[[1]] <- WJS.IPM$K; WJT.K[[1]] <- WJT.IPM$K  
WAS.P[[1]] <- WAS.IPM$P; WAT.P[[1]] <- WAT.IPM$P; WJS.P[[1]] <- WJS.IPM$P; WJT.P[[1]] <- WJT.IPM$P  
WAS.F[[1]] <- WAS.IPM$F; WAT.F[[1]] <- WAT.IPM$F; WJS.F[[1]] <- WJS.IPM$F; WJT.F[[1]] <- WJT.IPM$F 
# Meshpoints (these will remain consistent across coral types)
C.meshpts[[1]] <- CAS.IPM$meshpts 
ST.meshpts[[1]] <- STAS.IPM$meshpts
W.meshpts[[1]] <- WAS.IPM$meshpts
# Parameters
CAS.m.par[[1]] <- unlist(m.par.store[1,]); CAT.m.par[[1]] <- unlist(m.par.store[2,]); CJS.m.par[[1]] <- unlist(m.par.store[3,]); CJT.m.par[[1]] <- unlist(m.par.store[4,])  
STAS.m.par[[1]] <- unlist(m.par.store[5,]); STAT.m.par[[1]] <- unlist(m.par.store[6,]); STJS.m.par[[1]] <- unlist(m.par.store[7,]); STJT.m.par[[1]] <- unlist(m.par.store[8,])  
WAS.m.par[[1]] <- unlist(m.par.store[9,]); WAT.m.par[[1]] <- unlist(m.par.store[10,]); WJS.m.par[[1]] <- unlist(m.par.store[11,]); WJT.m.par[[1]] <- unlist(m.par.store[12,])  

# run the loop
for(ii in a:b){ # this will leave space to insert the main model details
  # print a progress marker
  print(ii)
  
  # determine the new Jackknife sample for both the main data frame and the fecundity data
  # the use of try should mean the loop will return NA's rather than breaking for errors but means that = signs can't be used.
  try(data.sample <- demographic.data[sample(nrow(demographic.data), 0.95*dim(demographic.data)[1], replace = F),]) # this takes a random 95% subset of the data.
  try(fec.sample <- fecundity.data[sample(nrow(fecundity.data), 0.95*dim(fecundity.data)[1], replace = F),])
  
  # Now force the models to follow the same regression formats as the 'original' models.
  
  ## 1. Survival -----------
  try(jsurv.mod <- glmer(surv ~ Size.t * LHS * Country * Ecoregion + (1|Colony.ID) + (Size.t|Site), family = "binomial", data = data.sample))

  ## 2. Growth -------------
  # A. Size transitions 
  try(jgrow.mod <- lmer(Size.t1 ~ Size.t + I(Size.t^2) * LHS * Country * Ecoregion + (1|Colony.ID) + (Size.t|Site), data = data.sample))
  # B. Growth variance
  try(jgrowth.sd.mod <- glmer(abs(resid(jgrow.mod)) ~ jgrow.mod@frame$Size.t * jgrow.mod@frame$LHS * jgrow.mod@frame$Country * jgrow.mod@frame$Ecoregion + (1|jgrow.mod@frame$Colony.ID) + (jgrow.mod@frame$Size.t1|jgrow.mod@frame$Site) , family = Gamma(link = "log")))
  
  ## 3. Fragmentation ------
  # Identify fragmenting colonies
  try(jcolonies.frag.TRUE <- data.sample[which(data.sample$frag == 1),]) 
  # Identify fragments produced
  try(jfrag.corals <- data.sample[which(data.sample$Fragment.t1 == "Yes"),])
  try(jfrag.corals <- fill(jfrag.corals, Size.t, .direction = "down"))
  # A. Fragmentation probability
  try(jfrag.mod.NL <- glmer(frag ~ Size.t + I(Size.t^2) * LHS * Country * Ecoregion + (1|Colony.ID) + (Size.t|Site), family = "binomial", data = data.sample))
  # B. Number of fragments
  try(jno.frag.mod <- glmer(No.frag.t1 ~ Size.t * LHS * Country + (1|Colony.ID) + (Size.t|Site), family = "poisson", data = jcolonies.frag.TRUE))
  # C. Size of fragments
  try(jfrag.size.mod.NL <- glm(Size.t1 ~ Size.t + I(Size.t^2) * LHS * Country, data = jfrag.corals))
  # D. Fragment size variance
  try(jfrag.size.sd.mod.NL <- glm(abs(resid(jfrag.size.mod.NL)) ~ jfrag.size.mod.NL$model$Size.t * jfrag.size.mod.NL$model$LHS * jfrag.size.mod.NL$model$Country, family = Gamma(link = "log")))
  
  ## 4. Recruitment --------
  # Identify recruits.
  try(jrecruit.data <- data.sample[which(data.sample$Recruit.t1 == "Yes"),])
  # A. Number of recruits
  # (This will be maintained from original data set and just allowed to vary within observed limits)
  # B. Recruit size
  try(jrecruit.size.mod2 <- lm(Size.t1 ~ LHS, data = jrecruit.data))
  
  ## 5. Fecundity ----------
  try(jfec.mod2 <- nls(Fecundity ~ a * exp(b * Size), data = fec.sample[which(fec.sample$LHS == "C"),], start = list(a = 1, b = 1)))
  try(jfec.mod3 <- nls(Fecundity ~ a * exp(b * Size), data = fec.sample[which(fec.sample$LHS == "ST"),], start = list(a = 1, b = 1)))
  try(jfec.mod4 <- nls(Fecundity ~ a * exp(b * Size), data = fec.sample[which(fec.sample$LHS == "W"),], start = list(a = 1, b = 1)))
  
  ##### Store vital rate parameters --------------------------------------------
  
  # Extract necessary mean coefficients.
  try(jsurvival.coefs <- apply(coef(jsurv.mod)$Colony.ID, 2, mean))
  try(jgrowth.coefs <- apply(coef(jgrow.mod)$Colony.ID, 2, mean))
  try(jgrowth.SD.coefs <- apply(coef(jgrowth.sd.mod)$'jgrow.mod@frame$Colony.ID', 2, mean))
  try(jgrowth.SD.coefs <- jgrowth.SD.coefs[2:25])
  try(jfrag.coefs <- apply(coef(jfrag.mod.NL)$Colony.ID, 2, mean))
  try(jfrag.no.coefs <- apply(coef(jno.frag.mod)$Colony.ID, 2, mean))
  
  ## Extract coefficients for each coral population
  
  ########################################################## Competitive ########
  
  # Competitive (Australia, Subtropical)
  try(jcoefs.C.A.S <- c(
    # Survival probability
    surv.int               =  jsurvival.coefs[1],
    surv.slope             =  jsurvival.coefs[2],
    # Growth transitions
    grow.int               =  jgrowth.coefs[1],
    grow.slope             =  jgrowth.coefs[2],
    grow.slope2            =  jgrowth.coefs[3],
    # Growth variance
    grow.sd.int            =  jgrowth.SD.coefs[1],
    grow.sd.slope          =  jgrowth.SD.coefs[2],
    # Fragmentation probability
    frag.int               =  jfrag.coefs[1],
    frag.slope             =  jfrag.coefs[2],
    frag.slope2            =  jfrag.coefs[3],
    # Number of frags produced  
    no.frag.int            =  jfrag.no.coefs[1],
    no.frag.slope          =  jfrag.no.coefs[2],
    # Fragment size
    frag.s.int             =  coef(jfrag.size.mod.NL)[1],
    frag.s.slope           =  coef(jfrag.size.mod.NL)[2],
    frag.s.slope2          =  coef(jfrag.size.mod.NL)[3],
    # Variance in fragment size
    frag.sd.int            =  coef(jfrag.size.sd.mod.NL)[1],
    frag.sd.slope          =  coef(jfrag.size.sd.mod.NL)[2],
    # Colony fecundity 
    fec.a                  =  coef(jfec.mod2)[1],
    fec.b                  =  coef(jfec.mod2)[2],
    # Number of recruits (this will allow the model to vary between observed limits)
    rec.max                =  m.par.store[1, "Rmax"],
    rec.min                =  m.par.store[1, "Rmin"],
    # recruit size coefficients
    rec.size.int           =  coefficients(jrecruit.size.mod2)[1],
    rec.size.sd            =  sigma(jrecruit.size.mod2)))
  
  # Competitive (Australia, Tropical)
  try(jcoefs.C.A.T <- c(
    # Survival probability
    surv.int               =  jsurvival.coefs[1] + jsurvival.coefs[6],
    surv.slope             =  jsurvival.coefs[2] + jsurvival.coefs[12],
    # Growth transitions
    grow.int               =  jgrowth.coefs[1] + jgrowth.coefs[7],
    grow.slope             =  jgrowth.coefs[2],
    grow.slope2            =  jgrowth.coefs[3] + jgrowth.coefs[13],
    # Growth variance
    grow.sd.int            =  jgrowth.SD.coefs[1] + jgrowth.SD.coefs[6],
    grow.sd.slope          =  jgrowth.SD.coefs[2] + jgrowth.SD.coefs[12],
    # Fragmentation probability
    frag.int               =  jfrag.coefs[1] + jfrag.coefs[7],
    frag.slope             =  jfrag.coefs[2],
    frag.slope2            =  jfrag.coefs[3] + jfrag.coefs[13],
    # Number of frags produced  
    no.frag.int            =  jfrag.no.coefs[1],
    no.frag.slope          =  jfrag.no.coefs[2],
    # Fragment size
    frag.s.int             =  coef(jfrag.size.mod.NL)[1],
    frag.s.slope           =  coef(jfrag.size.mod.NL)[2],
    frag.s.slope2          =  coef(jfrag.size.mod.NL)[3],
    # Variance in fragment size
    frag.sd.int            =  coef(jfrag.size.sd.mod.NL)[1],
    frag.sd.slope          =  coef(jfrag.size.sd.mod.NL)[2],
    # Colony fecundity 
    fec.a                  =  coef(jfec.mod2)[1],
    fec.b                  =  coef(jfec.mod2)[2],
    # Number of recruits (this will allow the model to vary between observed limits)
    rec.max                =  m.par.store[2, "Rmax"],
    rec.min                =  m.par.store[2, "Rmin"]*0.5,
    # recruit size coefficients
    rec.size.int           =  coefficients(jrecruit.size.mod2)[1],
    rec.size.sd            =  sigma(jrecruit.size.mod2)))
  
  # Competitive (Japan, Subtropical)
  try(jcoefs.C.J.S <- c(
    # Survival probability
    surv.int               =  jsurvival.coefs[1] + jsurvival.coefs[5],
    surv.slope             =  jsurvival.coefs[2] + jsurvival.coefs[9],
    # Growth transitions
    grow.int               =  jgrowth.coefs[1] + jgrowth.coefs[6],
    grow.slope             =  jgrowth.coefs[2],
    grow.slope2            =  jgrowth.coefs[3] + jgrowth.coefs[10],
    # Growth variance
    grow.sd.int            =  jgrowth.SD.coefs[1] + jgrowth.SD.coefs[5],
    grow.sd.slope          =  jgrowth.SD.coefs[2] + jgrowth.SD.coefs[9],
    # Fragmentation probability
    frag.int               =  jfrag.coefs[1] + jfrag.coefs[6],
    frag.slope             =  jfrag.coefs[2],
    frag.slope2            =  jfrag.coefs[3] + jfrag.coefs[10],
    # Number of frags produced  
    no.frag.int            =  jfrag.no.coefs[1] + jfrag.no.coefs[5],
    no.frag.slope          =  jfrag.no.coefs[2] + jfrag.no.coefs[8],
    # Fragment size
    frag.s.int             =  coef(jfrag.size.mod.NL)[1] + coef(jfrag.size.mod.NL)[6],
    frag.s.slope           =  coef(jfrag.size.mod.NL)[2],
    frag.s.slope2          =  coef(jfrag.size.mod.NL)[3] + coef(jfrag.size.mod.NL)[9],
    # Variance in fragment size
    frag.sd.int            =  coef(jfrag.size.sd.mod.NL)[1] + coef(jfrag.size.sd.mod.NL)[5],
    frag.sd.slope          =  coef(jfrag.size.sd.mod.NL)[2] + coef(jfrag.size.sd.mod.NL)[8],
    # Colony fecundity 
    fec.a                  =  coef(jfec.mod2)[1],
    fec.b                  =  coef(jfec.mod2)[2],
    # Number of recruits (this will allow the model to vary between observed limits)
    rec.max                =  m.par.store[3, "Rmax"],
    rec.min                =  m.par.store[3, "Rmin"],
    # recruit size coefficients
    rec.size.int           =  coefficients(jrecruit.size.mod2)[1],
    rec.size.sd            =  sigma(jrecruit.size.mod2)))
  
  # Competitive (Japan, Tropical)
  try(jcoefs.C.J.T <- c(
    # Survival probability
    surv.int               =  jsurvival.coefs[1] + jsurvival.coefs[15],
    surv.slope             =  jsurvival.coefs[2] + jsurvival.coefs[20],
    # Growth transitions
    grow.int               =  jgrowth.coefs[1] + jgrowth.coefs[16],
    grow.slope             =  jgrowth.coefs[2],
    grow.slope2            =  jgrowth.coefs[3] + jgrowth.coefs[21],
    # Growth variance
    grow.sd.int            =  jgrowth.SD.coefs[1] + jgrowth.SD.coefs[15],
    grow.sd.slope          =  jgrowth.SD.coefs[2] + jgrowth.SD.coefs[20],
    # Fragmentation probability
    frag.int               =  jfrag.coefs[1] + jfrag.coefs[16],
    frag.slope             =  jfrag.coefs[2],
    frag.slope2            =  jfrag.coefs[3] + jfrag.coefs[21],
    # Number of frags produced  
    no.frag.int            =  jfrag.no.coefs[1] + jfrag.no.coefs[5],
    no.frag.slope          =  jfrag.no.coefs[2] + jfrag.no.coefs[8],
    # Fragment size
    frag.s.int             =  coef(jfrag.size.mod.NL)[1] + coef(jfrag.size.mod.NL)[6],
    frag.s.slope           =  coef(jfrag.size.mod.NL)[2],
    frag.s.slope2          =  coef(jfrag.size.mod.NL)[3] + coef(jfrag.size.mod.NL)[9],
    # Variance in fragment size
    frag.sd.int            =  coef(jfrag.size.sd.mod.NL)[1] + coef(jfrag.size.sd.mod.NL)[5],
    frag.sd.slope          =  coef(jfrag.size.sd.mod.NL)[2] + coef(jfrag.size.sd.mod.NL)[8],
    # Colony fecundity 
    fec.a                  =  coef(jfec.mod2)[1],
    fec.b                  =  coef(jfec.mod2)[2],
    # Number of recruits (this will allow the model to vary between observed limits)
    rec.max                =  m.par.store[4, "Rmax"],
    rec.min                =  m.par.store[4, "Rmin"],
    # recruit size coefficients
    rec.size.int           =  coefficients(jrecruit.size.mod2)[1],
    rec.size.sd            =  sigma(jrecruit.size.mod2)))
  
  ########################################################## Stress-Tolerant ########
  
  # Stress-Tolerant (Australia, Subtropical)
  try(jcoefs.ST.A.S <- c(
    # Survival probability
    surv.int               =  jsurvival.coefs[1] + jsurvival.coefs[3],
    surv.slope             =  jsurvival.coefs[2] + jsurvival.coefs[7],
    # Growth transitions
    grow.int               =  jgrowth.coefs[1] + jgrowth.coefs[4],
    grow.slope             =  jgrowth.coefs[2],
    grow.slope2            =  jgrowth.coefs[3] + jgrowth.coefs[8],
    # Growth variance
    grow.sd.int            =  jgrowth.SD.coefs[1] + jgrowth.SD.coefs[3],
    grow.sd.slope          =  jgrowth.SD.coefs[2] + jgrowth.SD.coefs[7],
    # Fragmentation probability
    frag.int               =  jfrag.coefs[1] + jfrag.coefs[4],
    frag.slope             =  jfrag.coefs[2],
    frag.slope2            =  jfrag.coefs[3] + jfrag.coefs[8],
    # Number of frags produced  
    no.frag.int            =  jfrag.no.coefs[1] + jfrag.no.coefs[3],
    no.frag.slope          =  jfrag.no.coefs[2] + jfrag.no.coefs[6],
    # Fragment size
    frag.s.int             =  coef(jfrag.size.mod.NL)[1] + coef(jfrag.size.mod.NL)[4],
    frag.s.slope           =  coef(jfrag.size.mod.NL)[2],
    frag.s.slope2          =  coef(jfrag.size.mod.NL)[3] + coef(jfrag.size.mod.NL)[7],
    # Variance in fragment size
    frag.sd.int            =  coef(jfrag.size.sd.mod.NL)[1] + coef(jfrag.size.sd.mod.NL)[3],
    frag.sd.slope          =  coef(jfrag.size.sd.mod.NL)[2] + coef(jfrag.size.sd.mod.NL)[6],
    # Colony fecundity 
    fec.a                  =  coef(jfec.mod3)[1],
    fec.b                  =  coef(jfec.mod3)[2],
    # Number of recruits (this will allow the model to vary between observed limits)
    rec.max                =  m.par.store[5, "Rmax"],
    rec.min                =  m.par.store[5, "Rmin"],
    # recruit size coefficients
    rec.size.int           =  coefficients(jrecruit.size.mod2)[1] + coefficients(jrecruit.size.mod2)[2],
    rec.size.sd            =  sigma(jrecruit.size.mod2)))
  
  # Stress-Tolerant (Australia, Tropical)
  try(jcoefs.ST.A.T <- c(
    # Survival probability
    surv.int               =  jsurvival.coefs[1] + jsurvival.coefs[13],
    surv.slope             =  jsurvival.coefs[2] + jsurvival.coefs[18],
    # Growth transitions
    grow.int               =  jgrowth.coefs[1] + jgrowth.coefs[14],
    grow.slope             =  jgrowth.coefs[2],
    grow.slope2            =  jgrowth.coefs[3] + jgrowth.coefs[19],
    # Growth variance
    grow.sd.int            =  jgrowth.SD.coefs[1] + jgrowth.SD.coefs[13],
    grow.sd.slope          =  jgrowth.SD.coefs[2] + jgrowth.SD.coefs[18],
    # Fragmentation probability
    frag.int               =  jfrag.coefs[1] + jfrag.coefs[14],
    frag.slope             =  jfrag.coefs[2],
    frag.slope2            =  jfrag.coefs[3] + jfrag.coefs[19],
    # Number of frags produced  
    no.frag.int            =  jfrag.no.coefs[1] + jfrag.no.coefs[3],
    no.frag.slope          =  jfrag.no.coefs[2] + jfrag.no.coefs[6],
    # Fragment size
    frag.s.int             =  coef(jfrag.size.mod.NL)[1] + coef(jfrag.size.mod.NL)[4],
    frag.s.slope           =  coef(jfrag.size.mod.NL)[2],
    frag.s.slope2          =  coef(jfrag.size.mod.NL)[3] + coef(jfrag.size.mod.NL)[7],
    # Variance in fragment size
    frag.sd.int            =  coef(jfrag.size.sd.mod.NL)[1] + coef(jfrag.size.sd.mod.NL)[3],
    frag.sd.slope          =  coef(jfrag.size.sd.mod.NL)[2] + coef(jfrag.size.sd.mod.NL)[6],
    # Colony fecundity 
    fec.a                  =  coef(jfec.mod3)[1],
    fec.b                  =  coef(jfec.mod3)[2],
    # Number of recruits (this will allow the model to vary between observed limits)
    rec.max                =  m.par.store[6, "Rmax"],
    rec.min                =  m.par.store[6, "Rmin"] * 0.5,
    # recruit size coefficients
    rec.size.int           =  coefficients(jrecruit.size.mod2)[1] + coefficients(jrecruit.size.mod2)[2],
    rec.size.sd            =  sigma(jrecruit.size.mod2)))
  
  # Stress-Tolerant (Japan, Subtropical)
  try(jcoefs.ST.J.S <- c(
    # Survival probability
    surv.int               =  jsurvival.coefs[1] + jsurvival.coefs[10],
    surv.slope             =  jsurvival.coefs[2] + jsurvival.coefs[16],
    # Growth transitions
    grow.int               =  jgrowth.coefs[1] + jgrowth.coefs[11],
    grow.slope             =  jgrowth.coefs[2],
    grow.slope2            =  jgrowth.coefs[3] + jgrowth.coefs[17],
    # Growth variance
    grow.sd.int            =  jgrowth.SD.coefs[1] + jgrowth.SD.coefs[10],
    grow.sd.slope          =  jgrowth.SD.coefs[2] + jgrowth.SD.coefs[16],
    # Fragmentation probability
    frag.int               =  jfrag.coefs[1] + jfrag.coefs[11],
    frag.slope             =  jfrag.coefs[2],
    frag.slope2            =  jfrag.coefs[3] + jfrag.coefs[17],
    # Number of frags produced  
    no.frag.int            =  jfrag.no.coefs[1] + jfrag.no.coefs[9],
    no.frag.slope          =  jfrag.no.coefs[2] + jfrag.no.coefs[11],
    # Fragment size
    frag.s.int             =  coef(jfrag.size.mod.NL)[1] + coef(jfrag.size.mod.NL)[10],
    frag.s.slope           =  coef(jfrag.size.mod.NL)[2],
    frag.s.slope2          =  coef(jfrag.size.mod.NL)[3] + coef(jfrag.size.mod.NL)[12],
    # Variance in fragment size
    frag.sd.int            =  coef(jfrag.size.sd.mod.NL)[1] + coef(jfrag.size.sd.mod.NL)[9],
    frag.sd.slope          =  coef(jfrag.size.sd.mod.NL)[2] + coef(jfrag.size.sd.mod.NL)[11],
    # Colony fecundity 
    fec.a                  =  coef(jfec.mod3)[1],
    fec.b                  =  coef(jfec.mod3)[2],
    # Number of recruits (this will allow the model to vary between observed limits)
    rec.max                =  m.par.store[7, "Rmax"],
    rec.min                =  m.par.store[7, "Rmin"],
    # recruit size coefficients
    rec.size.int           =  coefficients(jrecruit.size.mod2)[1] + coefficients(jrecruit.size.mod2)[2],
    rec.size.sd            =  sigma(jrecruit.size.mod2)))
  
  # Stress-Tolerant (Japan, Tropical)
  try(jcoefs.ST.J.T <- c(
    # Survival probability
    surv.int               =  jsurvival.coefs[1] + jsurvival.coefs[21],
    surv.slope             =  jsurvival.coefs[2] + jsurvival.coefs[23],
    # Growth transitions
    grow.int               =  jgrowth.coefs[1] + jgrowth.coefs[22],
    grow.slope             =  jgrowth.coefs[2],
    grow.slope2            =  jgrowth.coefs[3] + jgrowth.coefs[24],
    # Growth variance
    grow.sd.int            =  jgrowth.SD.coefs[1] + jgrowth.SD.coefs[21],
    grow.sd.slope          =  jgrowth.SD.coefs[2] + jgrowth.SD.coefs[23],
    # Fragmentation probability
    frag.int               =  jfrag.coefs[1] + jfrag.coefs[22],
    frag.slope             =  jfrag.coefs[2],
    frag.slope2            =  jfrag.coefs[3] + jfrag.coefs[24],
    # Number of frags produced  
    no.frag.int            =  jfrag.no.coefs[1] + jfrag.no.coefs[9],
    no.frag.slope          =  jfrag.no.coefs[2] + jfrag.no.coefs[11],
    # Fragment size
    frag.s.int             =  coef(jfrag.size.mod.NL)[1] + coef(jfrag.size.mod.NL)[10],
    frag.s.slope           =  coef(jfrag.size.mod.NL)[2],
    frag.s.slope2          =  coef(jfrag.size.mod.NL)[3] + coef(jfrag.size.mod.NL)[12],
    # Variance in fragment size
    frag.sd.int            =  coef(jfrag.size.sd.mod.NL)[1] + coef(jfrag.size.sd.mod.NL)[9],
    frag.sd.slope          =  coef(jfrag.size.sd.mod.NL)[2] + coef(jfrag.size.sd.mod.NL)[11],
    # Colony fecundity 
    fec.a                  =  coef(jfec.mod3)[1],
    fec.b                  =  coef(jfec.mod3)[2],
    # Number of recruits (this will allow the model to vary between observed limits)
    rec.max                =  m.par.store[8, "Rmax"],
    rec.min                =  m.par.store[8, "Rmin"],
    # recruit size coefficients
    rec.size.int           =  coefficients(jrecruit.size.mod2)[1] + coefficients(jrecruit.size.mod2)[2],
    rec.size.sd            =  sigma(jrecruit.size.mod2)))
  
  ########################################################## Weedy ########
  
  # Weedy (Australia, Subtropical)
  try(jcoefs.W.A.S <- c(
    # Survival probability
    surv.int               =  jsurvival.coefs[1] + jsurvival.coefs[4],
    surv.slope             =  jsurvival.coefs[2] + jsurvival.coefs[8],
    # Growth transitions
    grow.int               =  jgrowth.coefs[1] + jgrowth.coefs[5],
    grow.slope             =  jgrowth.coefs[2],
    grow.slope2            =  jgrowth.coefs[3] + jgrowth.coefs[9],
    # Growth variance
    grow.sd.int            =  jgrowth.SD.coefs[1] + jgrowth.SD.coefs[4],
    grow.sd.slope          =  jgrowth.SD.coefs[2] + jgrowth.SD.coefs[8],
    # Fragmentation probability
    frag.int               =  jfrag.coefs[1] + jfrag.coefs[5],
    frag.slope             =  jfrag.coefs[2],
    frag.slope2            =  jfrag.coefs[3] + jfrag.coefs[9],
    # Number of frags produced  
    no.frag.int            =  jfrag.no.coefs[1] + jfrag.no.coefs[4],
    no.frag.slope          =  jfrag.no.coefs[2] + jfrag.no.coefs[7],
    # Fragment size
    frag.s.int             =  coef(jfrag.size.mod.NL)[1] + coef(jfrag.size.mod.NL)[5],
    frag.s.slope           =  coef(jfrag.size.mod.NL)[2],
    frag.s.slope2          =  coef(jfrag.size.mod.NL)[3] + coef(jfrag.size.mod.NL)[8],
    # Variance in fragment size
    frag.sd.int            =  coef(jfrag.size.sd.mod.NL)[1] + coef(jfrag.size.sd.mod.NL)[4],
    frag.sd.slope          =  coef(jfrag.size.sd.mod.NL)[2] + coef(jfrag.size.sd.mod.NL)[7],
    # Colony fecundity 
    fec.a                  =  coef(jfec.mod4)[1],
    fec.b                  =  coef(jfec.mod4)[2],
    # Number of recruits (this will allow the model to vary between observed limits)
    rec.max                =  m.par.store[9, "Rmax"],
    rec.min                =  m.par.store[9, "Rmin"],
    # recruit size coefficients
    rec.size.int           =  coefficients(jrecruit.size.mod2)[1] + coefficients(jrecruit.size.mod2)[3],
    rec.size.sd            =  sigma(jrecruit.size.mod2)))
  
  # Weedy (Australia, Tropical)
  try(jcoefs.W.A.T <- c(
    # Survival probability
    surv.int               =  jsurvival.coefs[1] + jsurvival.coefs[14],
    surv.slope             =  jsurvival.coefs[2] + jsurvival.coefs[19],
    # Growth transitions
    grow.int               =  jgrowth.coefs[1] + jgrowth.coefs[15],
    grow.slope             =  jgrowth.coefs[2],
    grow.slope2            =  jgrowth.coefs[3] + jgrowth.coefs[20],
    # Growth variance
    grow.sd.int            =  jgrowth.SD.coefs[1] + jgrowth.SD.coefs[14],
    grow.sd.slope          =  jgrowth.SD.coefs[2] + jgrowth.SD.coefs[19],
    # Fragmentation probability
    frag.int               =  jfrag.coefs[1] + jfrag.coefs[15],
    frag.slope             =  jfrag.coefs[2],
    frag.slope2            =  jfrag.coefs[3] + jfrag.coefs[20],
    # Number of frags produced  
    no.frag.int            =  jfrag.no.coefs[1] + jfrag.no.coefs[4],
    no.frag.slope          =  jfrag.no.coefs[2] + jfrag.no.coefs[7],
    # Fragment size
    frag.s.int             =  coef(jfrag.size.mod.NL)[1] + coef(jfrag.size.mod.NL)[5],
    frag.s.slope           =  coef(jfrag.size.mod.NL)[2],
    frag.s.slope2          =  coef(jfrag.size.mod.NL)[3] + coef(jfrag.size.mod.NL)[8],
    # Variance in fragment size
    frag.sd.int            =  coef(jfrag.size.sd.mod.NL)[1] + coef(jfrag.size.sd.mod.NL)[4],
    frag.sd.slope          =  coef(jfrag.size.sd.mod.NL)[2] + coef(jfrag.size.sd.mod.NL)[7],
    # Colony fecundity 
    fec.a                  =  coef(jfec.mod4)[1],
    fec.b                  =  coef(jfec.mod4)[2],
    # Number of recruits (this will allow the model to vary between observed limits)
    rec.max                =  m.par.store[10, "Rmax"],
    rec.min                =  m.par.store[10, "Rmin"] * 0.5,
    # recruit size coefficients
    rec.size.int           =  coefficients(jrecruit.size.mod2)[1] + coefficients(jrecruit.size.mod2)[3],
    rec.size.sd            =  sigma(jrecruit.size.mod2)))
  
  # Weedy (Japan, Subtropical)
  try(jcoefs.W.J.S <- c(
    # Survival probability
    surv.int               =  jsurvival.coefs[1] + jsurvival.coefs[11],
    surv.slope             =  jsurvival.coefs[2] + jsurvival.coefs[17],
    # Growth transitions
    grow.int               =  jgrowth.coefs[1] + jgrowth.coefs[12],
    grow.slope             =  jgrowth.coefs[2],
    grow.slope2            =  jgrowth.coefs[3] + jgrowth.coefs[18],
    # Growth variance
    grow.sd.int            =  jgrowth.SD.coefs[1] + jgrowth.SD.coefs[11],
    grow.sd.slope          =  jgrowth.SD.coefs[2] + jgrowth.SD.coefs[17],
    # Fragmentation probability
    frag.int               =  jfrag.coefs[1] + jfrag.coefs[12],
    frag.slope             =  jfrag.coefs[2],
    frag.slope2            =  jfrag.coefs[3] + jfrag.coefs[18],
    # Number of frags produced  
    no.frag.int            =  jfrag.no.coefs[1] + jfrag.no.coefs[10],
    no.frag.slope          =  jfrag.no.coefs[2] + jfrag.no.coefs[12],
    # Fragment size
    frag.s.int             =  coef(jfrag.size.mod.NL)[1] + coef(jfrag.size.mod.NL)[11],
    frag.s.slope           =  coef(jfrag.size.mod.NL)[2],
    frag.s.slope2          =  coef(jfrag.size.mod.NL)[3] + coef(jfrag.size.mod.NL)[13],
    # Variance in fragment size
    frag.sd.int            =  coef(jfrag.size.sd.mod.NL)[1] + coef(jfrag.size.sd.mod.NL)[10],
    frag.sd.slope          =  coef(jfrag.size.sd.mod.NL)[2] + coef(jfrag.size.sd.mod.NL)[12],
    # Colony fecundity 
    fec.a                  =  coef(jfec.mod4)[1],
    fec.b                  =  coef(jfec.mod4)[2],
    # Number of recruits (this will allow the model to vary between observed limits)
    rec.max                =  m.par.store[11, "Rmax"],
    rec.min                =  m.par.store[11, "Rmin"],
    # recruit size coefficients
    rec.size.int           =  coefficients(jrecruit.size.mod2)[1] + coefficients(jrecruit.size.mod2)[3],
    rec.size.sd            =  sigma(jrecruit.size.mod2)))
  
  # Weedy (Japan, Tropical)
  try(jcoefs.W.J.T <- c(
    # Survival probability
    surv.int               =  jsurvival.coefs[1] + jsurvival.coefs[22],
    surv.slope             =  jsurvival.coefs[2] + jsurvival.coefs[24],
    # Growth transitions
    grow.int               =  jgrowth.coefs[1] + jgrowth.coefs[23],
    grow.slope             =  jgrowth.coefs[2],
    grow.slope2            =  jgrowth.coefs[3] + jgrowth.coefs[25],
    # Growth variance
    grow.sd.int            =  jgrowth.SD.coefs[1] + jgrowth.SD.coefs[22],
    grow.sd.slope          =  jgrowth.SD.coefs[2] + jgrowth.SD.coefs[24],
    # Fragmentation probability
    frag.int               =  jfrag.coefs[1] + jfrag.coefs[23],
    frag.slope             =  jfrag.coefs[2],
    frag.slope2            =  jfrag.coefs[3] + jfrag.coefs[25],
    # Number of frags produced  
    no.frag.int            =  jfrag.no.coefs[1] + jfrag.no.coefs[10],
    no.frag.slope          =  jfrag.no.coefs[2] + jfrag.no.coefs[12],
    # Fragment size
    frag.s.int             =  coef(jfrag.size.mod.NL)[1] + coef(jfrag.size.mod.NL)[11],
    frag.s.slope           =  coef(jfrag.size.mod.NL)[2],
    frag.s.slope2          =  coef(jfrag.size.mod.NL)[3] + coef(jfrag.size.mod.NL)[13],
    # Variance in fragment size
    frag.sd.int            =  coef(jfrag.size.sd.mod.NL)[1] + coef(jfrag.size.sd.mod.NL)[10],
    frag.sd.slope          =  coef(jfrag.size.sd.mod.NL)[2] + coef(jfrag.size.sd.mod.NL)[12],
    # Colony fecundity 
    fec.a                  =  coef(jfec.mod4)[1],
    fec.b                  =  coef(jfec.mod4)[2],
    # Number of recruits (this will allow the model to vary between observed limits)
    rec.max                =  m.par.store[12, "Rmax"],
    rec.min                =  m.par.store[12, "Rmin"],
    # recruit size coefficients
    rec.size.int           =  coefficients(jrecruit.size.mod2)[1] + coefficients(jrecruit.size.mod2)[3],
    rec.size.sd            =  sigma(jrecruit.size.mod2)))
  
  ## Name and store variables
  # Combine together all regression parameters
  jm.par <- matrix(c(
    C_A_S = jcoefs.C.A.S,
    C_A_T = jcoefs.C.A.T,
    C_J_S = jcoefs.C.J.S,
    C_J_T = jcoefs.C.J.T,
    ST_A_S = jcoefs.ST.A.S,
    ST_A_T = jcoefs.ST.A.T,
    ST_J_S = jcoefs.ST.J.S,
    ST_J_T = jcoefs.ST.J.T,
    W_A_S = jcoefs.W.A.S,
    W_A_T = jcoefs.W.A.T,
    W_J_S = jcoefs.W.J.S,
    W_J_T = jcoefs.W.J.T), nrow = 12, ncol = 23, byrow = T)
  
  # and name columns and rows
  rownames(jm.par) <- c("C_A_S", "C_A_T", "C_J_S", "C_J_T", 
                        "ST_A_S", "ST_A_T", "ST_J_S", "ST_J_T",
                        "W_A_S", "W_A_T", "W_J_S", "W_J_T")
  
  colnames(jm.par) <-  c("surv.int", "surv.slope", # coefficients of survival
                         "grow.int", "grow.slope", "grow.slope2", # growth
                         "grow.sd.int", "grow.sd.slope", # variance in growth
                         "frag.int", "frag.slope", "frag.slope2", # fragmentation 
                         "frag.no.int", "frag.no.slope", # number of fragments
                         "frag.size.int", "frag.size.slope", "frag.size.slope2", # size of fragments
                         "frag.sd.int", "frag.sd.slope", # variance in fragment sizes
                         "fec.a", "fec.b", # colony fecundity
                         "Rmax", "Rmin", # recruit appearance
                         "rec.size", "rec.size.sd") # recruit sizes
  
  ##### Estimate model dimensions ----------------------------------------------
  
  # Competitive
  try(jC.L <- 1.1 * min(data.sample[which(data.sample$LHS == "C"), c("Size.t","Size.t1")], na.rm = T)) 
  try(jC.U <- 1.1 * max(data.sample[which(data.sample$LHS == "C"), c("Size.t","Size.t1")], na.rm = T)) 
  # Stress-Tolerant
  try(jST.L <- 1.1 * min(data.sample[which(data.sample$LHS == "ST"), c("Size.t","Size.t1")], na.rm = T)) 
  try(jST.U <- 1.1 * max(data.sample[which(data.sample$LHS == "ST"), c("Size.t","Size.t1")], na.rm = T)) 
  # Weedy
  try(jW.L <- 1.1 * min(data.sample[which(data.sample$LHS == "W"), c("Size.t","Size.t1")], na.rm = T)) 
  try(jW.U <- 1.1 * max(data.sample[which(data.sample$LHS == "W"), c("Size.t","Size.t1")], na.rm = T)) 
  # And store the parameters together
  try(jL.store <- c(jC.L, jST.L, jW.L))
  try(jU.store <- c(jC.U, jST.U, jW.U))
  
  ##### Create Jack-knifed IPMs ------------------------------------------------
  
  # Competitive
  try(jCAS.IPM <- mk_K(m = m, m.par = unlist(jm.par[1,]), L = jL.store[1], U = jU.store[1]))
  try(jCAT.IPM <- mk_K(m = m, m.par = unlist(jm.par[2,]), L = jL.store[1], U = jU.store[1]))
  try(jCJS.IPM <- mk_K(m = m, m.par = unlist(jm.par[3,]), L = jL.store[1], U = jU.store[1]))
  try(jCJT.IPM <- mk_K(m = m, m.par = unlist(jm.par[4,]), L = jL.store[1], U = jU.store[1]))
  # Stress Tolerant
  try(jSTAS.IPM <- mk_K(m = m, m.par = unlist(jm.par[5,]), L = jL.store[2], U = jU.store[2]))
  try(jSTAT.IPM <- mk_K(m = m, m.par = unlist(jm.par[6,]), L = jL.store[2], U = jU.store[2]))
  try(jSTJS.IPM <- mk_K(m = m, m.par = unlist(jm.par[7,]), L = jL.store[2], U = jU.store[2]))
  try(jSTJT.IPM <- mk_K(m = m, m.par = unlist(jm.par[8,]), L = jL.store[2], U = jU.store[2]))
  # Weedy
  try(jWAS.IPM <- mk_K(m = m, m.par = unlist(jm.par[9,]), L = jL.store[3], U = jU.store[3]))
  try(jWAT.IPM <- mk_K(m = m, m.par = unlist(jm.par[10,]), L = jL.store[3], U = jU.store[3]))
  try(jWJS.IPM <- mk_K(m = m, m.par = unlist(jm.par[11,]), L = jL.store[3], U = jU.store[3]))
  try(jWJT.IPM <- mk_K(m = m, m.par = unlist(jm.par[12,]), L = jL.store[3], U = jU.store[3]))
   
  ##### Store Model outputs ----------------------------------------------------
  
  ### Lambda 
  # Competitive
  try(C.outputs[ii,1] <- Re(eigen(jCAS.IPM$K)$values)[1])
  try(C.outputs[ii,2] <- Re(eigen(jCAT.IPM$K)$values)[1])
  try(C.outputs[ii,3] <- Re(eigen(jCJS.IPM$K)$values)[1])
  try(C.outputs[ii,4] <- Re(eigen(jCJT.IPM$K)$values)[1])
  # Stress Tolerant
  try(ST.outputs[ii,1] <- Re(eigen(jSTAS.IPM$K)$values)[1])
  try(ST.outputs[ii,2] <- Re(eigen(jSTAT.IPM$K)$values)[1])
  try(ST.outputs[ii,3] <- Re(eigen(jSTJS.IPM$K)$values)[1])
  try(ST.outputs[ii,4] <- Re(eigen(jSTJT.IPM$K)$values)[1])
  # Weedy
  try(W.outputs[ii,1] <- Re(eigen(jWAS.IPM$K)$values)[1])
  try(W.outputs[ii,2] <- Re(eigen(jWAT.IPM$K)$values)[1])
  try(W.outputs[ii,3] <- Re(eigen(jWJS.IPM$K)$values)[1])
  try(W.outputs[ii,4] <- Re(eigen(jWJT.IPM$K)$values)[1])
  
  ### Kreiss Bounds
  # Competitive
  try(C.outputs[ii,5] <- Kreiss(jCAS.IPM$K, bound = "upper")); try(C.outputs[ii,9] <- Kreiss(jCAS.IPM$K, bound = "lower"))
  try(C.outputs[ii,6] <- Kreiss(jCAT.IPM$K, bound = "upper")); try(C.outputs[ii,10] <- Kreiss(jCAT.IPM$K, bound = "lower"))
  try(C.outputs[ii,7] <- Kreiss(jCJS.IPM$K, bound = "upper")); try(C.outputs[ii,11] <- Kreiss(jCJS.IPM$K, bound = "lower"))
  try(C.outputs[ii,8] <- Kreiss(jCJT.IPM$K, bound = "upper")); try(C.outputs[ii,12] <- Kreiss(jCJT.IPM$K, bound = "lower"))
  # Stress-Tolerant
  try(ST.outputs[ii,5] <- Kreiss(jSTAS.IPM$K, bound = "upper")); try(ST.outputs[ii,9] <- Kreiss(jSTAS.IPM$K, bound = "lower"))
  try(ST.outputs[ii,6] <- Kreiss(jSTAT.IPM$K, bound = "upper")); try(ST.outputs[ii,10] <- Kreiss(jSTAT.IPM$K, bound = "lower"))
  try(ST.outputs[ii,7] <- Kreiss(jSTJS.IPM$K, bound = "upper")); try(ST.outputs[ii,11] <- Kreiss(jSTJS.IPM$K, bound = "lower"))
  try(ST.outputs[ii,8] <- Kreiss(jSTJT.IPM$K, bound = "upper")); try(ST.outputs[ii,12] <- Kreiss(jSTJT.IPM$K, bound = "lower"))
  # Weedy
  try(W.outputs[ii,5] <- Kreiss(jWAS.IPM$K, bound = "upper")); try(W.outputs[ii,9] <- Kreiss(jWAS.IPM$K, bound = "lower"))
  try(W.outputs[ii,6] <- Kreiss(jWAT.IPM$K, bound = "upper")); try(W.outputs[ii,10] <- Kreiss(jWAT.IPM$K, bound = "lower"))
  try(W.outputs[ii,7] <- Kreiss(jWJS.IPM$K, bound = "upper")); try(W.outputs[ii,11] <- Kreiss(jWJS.IPM$K, bound = "lower"))
  try(W.outputs[ii,8] <- Kreiss(jWJT.IPM$K, bound = "upper")); try(W.outputs[ii,12] <- Kreiss(jWJT.IPM$K, bound = "lower"))
  
  ### Transient envelopes 
  # Competitive
  try(C.outputs[ii,13] <- C.outputs[ii,5] - C.outputs[ii,9])
  try(C.outputs[ii,14] <- C.outputs[ii,6] - C.outputs[ii,10])
  try(C.outputs[ii,15] <- C.outputs[ii,7] - C.outputs[ii,11])
  try(C.outputs[ii,16] <- C.outputs[ii,8] - C.outputs[ii,12])
  # Stress-Tolerant
  try(ST.outputs[ii,13] <- ST.outputs[ii,5] - ST.outputs[ii,9])
  try(ST.outputs[ii,14] <- ST.outputs[ii,6] - ST.outputs[ii,10])
  try(ST.outputs[ii,15] <- ST.outputs[ii,7] - ST.outputs[ii,11])
  try(ST.outputs[ii,16] <- ST.outputs[ii,8] - ST.outputs[ii,12])
  # Weedy
  try(W.outputs[ii,13] <- W.outputs[ii,5] - W.outputs[ii,9])
  try(W.outputs[ii,14] <- W.outputs[ii,6] - W.outputs[ii,10])
  try(W.outputs[ii,15] <- W.outputs[ii,7] - W.outputs[ii,11])
  try(W.outputs[ii,16] <- W.outputs[ii,8] - W.outputs[ii,12])
  
  # Associated recruitment value.
  # Competitive
  try(C.outputs[ii,17] <- jCAS.IPM$R)
  try(C.outputs[ii,18] <- jCAT.IPM$R)
  try(C.outputs[ii,19] <- jCJS.IPM$R)
  try(C.outputs[ii,20] <- jCJT.IPM$R)
  # Stress-Tolerant
  try(ST.outputs[ii,17] <- jSTAS.IPM$R)
  try(ST.outputs[ii,18] <- jSTAT.IPM$R)
  try(ST.outputs[ii,19] <- jSTJS.IPM$R)
  try(ST.outputs[ii,20] <- jSTJT.IPM$R)
  # Weedy
  try(W.outputs[ii,17] <- jWAS.IPM$R)
  try(W.outputs[ii,18] <- jWAT.IPM$R)
  try(W.outputs[ii,19] <- jWJS.IPM$R)
  try(W.outputs[ii,20] <- jWJT.IPM$R)
  
  # estimate reactivity
  # Competitive
  try(C.outputs[ii,21] <- reac(jCAS.IPM$K, bound = "upper"))
  try(C.outputs[ii,22] <- reac(jCAT.IPM$K, bound = "upper"))
  try(C.outputs[ii,23] <- reac(jCJS.IPM$K, bound = "upper"))
  try(C.outputs[ii,24] <- reac(jCJT.IPM$K, bound = "upper"))
  # Stress-Tolerant
  try(ST.outputs[ii,21] <- reac(jSTAS.IPM$K, bound = "upper"))
  try(ST.outputs[ii,22] <- reac(jSTAT.IPM$K, bound = "upper"))
  try(ST.outputs[ii,23] <- reac(jSTJS.IPM$K, bound = "upper"))
  try(ST.outputs[ii,24] <- reac(jSTJT.IPM$K, bound = "upper"))
  # Weedy
  try(W.outputs[ii,21] <- reac(jWAS.IPM$K, bound = "upper"))
  try(W.outputs[ii,22] <- reac(jWAT.IPM$K, bound = "upper"))
  try(W.outputs[ii,23] <- reac(jWJS.IPM$K, bound = "upper"))
  try(W.outputs[ii,24] <- reac(jWJT.IPM$K, bound = "upper"))
  
  # estimate attenuation
  # Competitive
  try(C.outputs[ii,25] <- reac(jCAS.IPM$K, bound = "lower"))
  try(C.outputs[ii,26] <- reac(jCAT.IPM$K, bound = "lower"))
  try(C.outputs[ii,27] <- reac(jCJS.IPM$K, bound = "lower"))
  try(C.outputs[ii,28] <- reac(jCJT.IPM$K, bound = "lower"))
  # Stress-Tolerant
  try(ST.outputs[ii,25] <- reac(jSTAS.IPM$K, bound = "lower"))
  try(ST.outputs[ii,26] <- reac(jSTAT.IPM$K, bound = "lower"))
  try(ST.outputs[ii,27] <- reac(jSTJS.IPM$K, bound = "lower"))
  try(ST.outputs[ii,28] <- reac(jSTJT.IPM$K, bound = "lower"))
  # Weedy
  try(W.outputs[ii,25] <- reac(jWAS.IPM$K, bound = "lower"))
  try(W.outputs[ii,26] <- reac(jWAT.IPM$K, bound = "lower"))
  try(W.outputs[ii,27] <- reac(jWJS.IPM$K, bound = "lower"))
  try(W.outputs[ii,28] <- reac(jWJT.IPM$K, bound = "lower"))
  
  # estimate maximal amplification
  # Competitive
  try(C.outputs[ii,29] <- maxamp(jCAS.IPM$K))
  try(C.outputs[ii,30] <- maxamp(jCAT.IPM$K))
  try(C.outputs[ii,31] <- maxamp(jCJS.IPM$K))
  try(C.outputs[ii,32] <- maxamp(jCJT.IPM$K))
  # Stress-Tolerant
  try(ST.outputs[ii,29] <- maxamp(jSTAS.IPM$K))
  try(ST.outputs[ii,30] <- maxamp(jSTAT.IPM$K))
  try(ST.outputs[ii,21] <- maxamp(jSTJS.IPM$K))
  try(ST.outputs[ii,32] <- maxamp(jSTJT.IPM$K))
  # Weedy
  try(W.outputs[ii,29] <- maxamp(jWAS.IPM$K))
  try(W.outputs[ii,30] <- maxamp(jWAT.IPM$K))
  try(W.outputs[ii,31] <- maxamp(jWJS.IPM$K))
  try(W.outputs[ii,32] <- maxamp(jWJT.IPM$K))
  
  # and finally estimate maximal attenuation
  try(C.outputs[ii,33] <- maxatt(jCAS.IPM$K))
  try(C.outputs[ii,34] <- maxatt(jCAT.IPM$K))
  try(C.outputs[ii,35] <- maxatt(jCJS.IPM$K))
  try(C.outputs[ii,36] <- maxatt(jCJT.IPM$K))
  # Stress-Tolerant
  try(ST.outputs[ii,33] <- maxatt(jSTAS.IPM$K))
  try(ST.outputs[ii,34] <- maxatt(jSTAT.IPM$K))
  try(ST.outputs[ii,35] <- maxatt(jSTJS.IPM$K))
  try(ST.outputs[ii,36] <- maxatt(jSTJT.IPM$K))
  # Weedy
  try(W.outputs[ii,33] <- maxatt(jWAS.IPM$K))
  try(W.outputs[ii,34] <- maxatt(jWAT.IPM$K))
  try(W.outputs[ii,35] <- maxatt(jWJS.IPM$K))
  try(W.outputs[ii,36] <- maxatt(jWJT.IPM$K))
  
  # Kernels
  try(CAS.K[[ii]] <- jCAS.IPM$K); try(CAT.K[[ii]] <- jCAT.IPM$K); try(CJS.K[[ii]] <- jCJS.IPM$K); try(CJT.K[[ii]] <- jCJT.IPM$K) 
  try(CAS.P[[ii]] <- jCAS.IPM$P); try(CAT.P[[ii]] <- jCAT.IPM$P); try(CJS.P[[ii]] <- jCJS.IPM$P); try(CJT.P[[ii]] <- jCJT.IPM$P) 
  try(CAS.F[[ii]] <- jCAS.IPM$F); try(CAT.F[[ii]] <- jCAT.IPM$F); try(CJS.F[[ii]] <- jCJS.IPM$F); try(CJT.F[[ii]] <- jCJT.IPM$F) 
  try(STAS.K[[ii]] <- jSTAS.IPM$K); try(STAT.K[[ii]] <- jSTAT.IPM$K); try(STJS.K[[ii]] <- jSTJS.IPM$K); try(STJT.K[[ii]] <- jSTJT.IPM$K)  
  try(STAS.P[[ii]] <- jSTAS.IPM$P); try(STAT.P[[ii]] <- jSTAT.IPM$P); try(STJS.P[[ii]] <- jSTJS.IPM$P); try(STJT.P[[ii]] <- jSTJT.IPM$P)  
  try(STAS.F[[ii]] <- jSTAS.IPM$F); try(STAT.F[[ii]] <- jSTAT.IPM$F); try(STJS.F[[ii]] <- jSTJS.IPM$F); try(STJT.F[[ii]] <- jSTJT.IPM$F)  
  try(WAS.K[[ii]] <- jWAS.IPM$K); try(WAT.K[[ii]] <- jWAT.IPM$K); try(WJS.K[[ii]] <- jWJS.IPM$K); try(WJT.K[[ii]] <- jWJT.IPM$K)  
  try(WAS.P[[ii]] <- jWAS.IPM$P); try(WAT.P[[ii]] <- jWAT.IPM$P); try(WJS.P[[ii]] <- jWJS.IPM$P); try(WJT.P[[ii]] <- jWJT.IPM$P)  
  try(WAS.F[[ii]] <- jWAS.IPM$F); try(WAT.F[[ii]] <- jWAT.IPM$F); try(WJS.F[[ii]] <- jWJS.IPM$F); try(WJT.F[[ii]] <- jWJT.IPM$F) 
  
  # Meshpoints
  try(C.meshpts[[ii]] <- jCAS.IPM$meshpts) 
  try(ST.meshpts[[ii]] <- jSTAS.IPM$meshpts)
  try(W.meshpts[[ii]] <- jWAS.IPM$meshpts)
  
  # Vital rate parameters
  try(CAS.m.par[[ii]] <- unlist(jm.par[1,])); try(CAT.m.par[[ii]] <- unlist(jm.par[2,])); try(CJS.m.par[[ii]] <- unlist(jm.par[3,])); try(CJT.m.par[[ii]] <- unlist(jm.par[4,]))  
  try(STAS.m.par[[ii]] <- unlist(jm.par[5,])); try(STAT.m.par[[ii]] <- unlist(jm.par[6,])); try(STJS.m.par[[ii]] <- unlist(jm.par[7,])); try(STJT.m.par[[ii]] <- unlist(jm.par[8,]))  
  try(WAS.m.par[[ii]] <- unlist(jm.par[9,])); try(WAT.m.par[[ii]] <- unlist(jm.par[10,])); try(WJS.m.par[[ii]] <- unlist(jm.par[11,])); try(WJT.m.par[[ii]] <- unlist(jm.par[12,]))  
  
  # And tidy up
  try(rm(data.sample, fec.sample, jfrag.corals, jcolonies.frag.TRUE, jsurv.mod, jgrow.mod, jgrow.sd.mod, jfrag.mod.NL, jno.frag.mod, jfrag.size.mod.NL, jfrag.size.sd.mod.NL, jfec.mod2, jfec.mod3, jfec.mod4, jrecruit.size.mod2, jrecruit.data,
         jsurvival.coefs, jgrowth.coefs, jgrowth.SD.coefs, jfrag.coefs, jfrag.no.coefs, jm.par,
         jcoefs.C.A.S, jcoefs.C.A.T, jcoefs.C.J.S, jcoefs.C.J.T, jcoefs.ST.A.S, jcoefs.ST.A.T, jcoefs.ST.J.S, jcoefs.ST.J.T, jcoefs.W.A.S, jcoefs.W.A.T, jcoefs.W.J.S, jcoefs.W.J.T,
         jC.L, jC.U, jST.L, jST.U, jW.L, jW.U, jL.store, jU.store,
         jCAS.IPM, jCAT.IPM, jCJS.IPM, jCJT.IPM, jSTAS.IPM, jSTAT.IPM, jSTJS.IPM, jSTJT.IPM, jWAS.IPM, jWAT.IPM, jWJS.IPM, jWJT.IPM))
}

# Concatenate the list outputs into dataframes for extraction
# Vital rate parameters
CAS.m.par <- ldply(CAS.m.par, rbind); CAT.m.par <- ldply(CAT.m.par, rbind); CJS.m.par <- ldply(CJS.m.par, rbind); CJT.m.par <- ldply(CJT.m.par, rbind)  
STAS.m.par <- ldply(STAS.m.par, rbind); STAT.m.par <- ldply(STAT.m.par, rbind); STJS.m.par <- ldply(STJS.m.par, rbind); STJT.m.par <- ldply(STJT.m.par, rbind)  
WAS.m.par <- ldply(WAS.m.par, rbind); WAT.m.par <- ldply(WAT.m.par, rbind); WJS.m.par <- ldply(WJS.m.par, rbind); WJT.m.par <- ldply(WJT.m.par, rbind)  

# Kernels
CAS.K <- ldply(CAS.K, rbind); CAT.K <- ldply(CAT.K, rbind); CJS.K <- ldply(CJS.K, rbind); CJT.K <- ldply(CJT.K, rbind) 
CAS.P <- ldply(CAS.P, rbind); CAT.P <- ldply(CAT.P, rbind); CJS.P <- ldply(CJS.P, rbind); CJT.P <- ldply(CJT.P, rbind) 
CAS.F <- ldply(CAS.F, rbind); CAT.F <- ldply(CAT.F, rbind); CJS.F <- ldply(CJS.F, rbind); CJT.F <- ldply(CJT.F, rbind) 
STAS.K <- ldply(STAS.K, rbind); STAT.K <- ldply(STAT.K, rbind); STJS.K <- ldply(STJS.K, rbind); STJT.K <- ldply(STJT.K, rbind)  
STAS.P <- ldply(STAS.P, rbind); STAT.P <- ldply(STAT.P, rbind); STJS.P <- ldply(STJS.P, rbind); STJT.P <- ldply(STJT.P, rbind)  
STAS.F <- ldply(STAS.F, rbind); STAT.F <- ldply(STAT.F, rbind); STJS.F <- ldply(STJS.F, rbind); STJT.F <- ldply(STJT.F, rbind)  
WAS.K <- ldply(WAS.K, rbind); WAT.K <- ldply(WAT.K, rbind); WJS.K <- ldply(WJS.K, rbind); WJT.K <- ldply(WJT.K, rbind)  
WAS.P <- ldply(WAS.P, rbind); WAT.P <- ldply(WAT.P, rbind); WJS.P <- ldply(WJS.P, rbind); WJT.P <- ldply(WJT.P, rbind)  
WAS.F <- ldply(WAS.F, rbind); WAT.F <- ldply(WAT.F, rbind); WJS.F <- ldply(WJS.F, rbind); WJT.F <- ldply(WJT.F, rbind) 

# Meshpoints
C.meshpts <- ldply(C.meshpts, rbind) 
ST.meshpts <- ldply(ST.meshpts, rbind)
W.meshpts <- ldply(W.meshpts, rbind)

# And extract
# Key data
write.csv(C.outputs, file = paste0("C outputs_", a, "-", b, ".csv", sep = ""), row.names = FALSE)
write.csv(ST.outputs, file = paste0("ST outputs_", a, "-", b, ".csv", sep = ""), row.names = FALSE)
write.csv(W.outputs, file = paste0("W outputs_", a, "-", b, ".csv", sep = ""), row.names = FALSE)

# Kernels
write.csv(CAS.K, file = paste0("CAS K kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(CAT.K, file = paste0("CAT K kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(CJS.K, file = paste0("CJS K kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(CJT.K, file = paste0("CJT K kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE) 
write.csv(CAS.P, file = paste0("CAS P kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(CAT.P, file = paste0("CAT P kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(CJS.P, file = paste0("CJS P kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(CJT.P, file = paste0("CJT P kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE) 
write.csv(CAS.F, file = paste0("CAS F kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(CAT.F, file = paste0("CAT F kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(CJS.F, file = paste0("CJS F kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(CJT.F, file = paste0("CJT F kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE) 
write.csv(STAS.K, file = paste0("STAS K kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(STAT.K, file = paste0("STAT K kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(STJS.K, file = paste0("STJS K kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(STJT.K, file = paste0("STJT K kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE) 
write.csv(STAS.P, file = paste0("STAS P kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(STAT.P, file = paste0("STAT P kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(STJS.P, file = paste0("STJS P kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(STJT.P, file = paste0("STJT P kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE) 
write.csv(STAS.F, file = paste0("STAS F kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(STAT.F, file = paste0("STAT F kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(STJS.F, file = paste0("STJS F kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(STJT.F, file = paste0("STJT F kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE) 
write.csv(WAS.K, file = paste0("WAS K kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(WAT.K, file = paste0("WAT K kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(WJS.K, file = paste0("WJS K kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(WJT.K, file = paste0("WJT K kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE) 
write.csv(WAS.P, file = paste0("WAS P kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(WAT.P, file = paste0("WAT P kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(WJS.P, file = paste0("WJS P kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(WJT.P, file = paste0("WJT P kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE) 
write.csv(WAS.F, file = paste0("WAS F kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(WAT.F, file = paste0("WAT F kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(WJS.F, file = paste0("WJS F kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(WJT.F, file = paste0("WJT F kernel_", a, "-", b, ".csv", sep = ""), row.names = FALSE) 

# Vital rate parameters
write.csv(CAS.m.par, file = paste0("CAS par_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(CAT.m.par, file = paste0("CAT par_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(CJS.m.par, file = paste0("CJS par_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(CJT.m.par, file = paste0("CJT par_", a, "-", b, ".csv", sep = ""), row.names = FALSE)  
write.csv(STAS.m.par, file = paste0("STAS par_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(STAT.m.par, file = paste0("STAT par_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(STJS.m.par, file = paste0("STJS par_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(STJT.m.par, file = paste0("STJT par_", a, "-", b, ".csv", sep = ""), row.names = FALSE)  
write.csv(WAS.m.par, file = paste0("WAS par_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(WAT.m.par, file = paste0("WAT par_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(WJS.m.par, file = paste0("WJS par_", a, "-", b, ".csv", sep = ""), row.names = FALSE); write.csv(WJT.m.par, file = paste0("WJT par_", a, "-", b, ".csv", sep = ""), row.names = FALSE)  

# Meshpoints
write.csv(C.meshpts, file = paste0("C meshpts_", a, "-", b, ".csv", sep = ""), rownames = FALSE)
write.csv(ST.meshpts, file = paste0("ST meshpts_", a, "-", b, ".csv", sep = ""), rownames = FALSE)
write.csv(W.meshpts, file = paste0("W meshpts_", a, "-", b, ".csv", sep = ""), rownames = FALSE)

########################################################## End of Code ##############################################################