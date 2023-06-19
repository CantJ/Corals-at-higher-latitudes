# This script is evaluating and plotting Generation time, as a function of lambda and transient potential.
# This script relies on many on the outputs generated using the other scripts contained within this repository.
# A key part of this script is the calculation of generation time estimates for each of the focal coral assemblages.

# Author: James Cant (jic2@st-andrews.ac.uk)

# ----------------------------------------------------------------------------------

# clear workspace
rm(list=ls(all=TRUE))

# Load required packages
library(Rage)
library(caret)
library(lmodel2)
library(fields)
library(popdemo)
library(plyr)
library(stats)
library(tidyr)
library(lme4)
library(ggplot2)
library(Rmisc)

############################################
# STEP 1: Extract previously calculated lambda and transient potential estimates
############################################

# reset working directory for data extraction
setwd("data_file_location")

# define an indexing vector for condensing all the data below.
index_start <- seq(1,1000, by = 20)
index_end <- index_start + 19

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

#### and remove un-needed objects to save memory space
rm(data_list, data_list1, data_list2, x, start, end, ii, index, index_end, index_start) 

# Now stitch together the relevant data
GT_data <- data.frame("LHS" = rep(c("C", "ST", "W"), each = 4000),
                      "Country" = rep(rep(c("Australia","Japan"), each = 2000), length.out = 12000),
                      "Region" = rep(rep(c("Subtropical","Tropical"), each = 1000), length.out = 12000),
                      "GT" = rep(NA, 12000),
                      "Lambda" = c(C.outputs[,"AS Lambda"], C.outputs[,"AT Lambda"], C.outputs[,"JS Lambda"], C.outputs[,"JT Lambda"],
                                   ST.outputs[,"AS Lambda"], ST.outputs[,"AT Lambda"], ST.outputs[,"JS Lambda"], ST.outputs[,"JT Lambda"],
                                   W.outputs[,"AS Lambda"], W.outputs[,"AT Lambda"], W.outputs[,"JS Lambda"], W.outputs[,"JT Lambda"]),
                      "TE" = c(C.outputs[,"AS TE"], C.outputs[,"AT TE"], C.outputs[,"JS TE"], C.outputs[,"JT TE"],
                               ST.outputs[,"AS TE"], ST.outputs[,"AT TE"], ST.outputs[,"JS TE"], ST.outputs[,"JT TE"],
                               W.outputs[,"AS TE"], W.outputs[,"AT TE"], W.outputs[,"JS TE"], W.outputs[,"JT TE"]))

#######################################
# STEP 2: Calculate Generation time
#######################################

# reset working directory
setwd("data_file_location")

# Load required data
m.par.store <- read.csv("vital rate parameters.csv", row.names = 1)
demographic.data <- read.csv("Formatted data.csv", stringsAsFactors = TRUE)
fecundity.data <- read.csv("Fecundity data.csv", stringsAsFactors = TRUE)

# set random number seed and define indexing for the loop below
set.seed(43400)
a = 2
b = 1000

# 1. Define function for discretising and combining the survival + growth (G) recruitment (F) kernels
# Here I am only recreating the models to calculate net reproductive rate and so I am only interested in survival, growth and reproduction.
mk_GF <- function(m, m.par, L, U) {
  # define the mesh points 
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  # determine the number of new larvae settling 
  recruit.no <- recruits(m.par)
  # Define the survival and growth kernel
  G <- h * (outer(meshpts, meshpts, G_z1z, m.par))
  # Define the recruitment kernel
  F <- h * (outer(meshpts, meshpts, F_z1z, m.par, recruit.no))
  return(list(meshpts = meshpts, G = G, F = F, R = recruit.no))
} # The model is also returning the number of recruits added so that I can 
# store them and have them associated with corresponding lambda/Kreiss values.

# 2. Define the survival/growth kernel 
# For this kernel I am only interested in the survival and growth of colonies.
G_z1z <- function (Size.t1, Size.t, m.par) {
  return( s_z(Size.t, m.par) * g_z1z(Size.t1, Size.t, m.par) ) 
}             # survive            # and grow/shrink


# 3. Define the production of new recruits into the populations.
F_z1z <- function (Size.t1, Size.t, m.par, recruit.no) {
  R <- settleR(Size.t, m.par, recruit.no) # estimate the proportional fecundity of colonies
  return( R * c0_z1(Size.t1, m.par) )
}     # Number of recruits # Size of successful recruits


# define vital rate functions -----------------------------------------------

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


# 3. A function that randomly selects a number of settling larvae (from within observed limits)
recruits <- function(m.par){
  R1 <- m.par["Rmax"]
  R2 <- m.par["Rmin"]
  rawR <- runif(1, min = R2, max = R1) # generate a random estimate from with the estimated range
  roundR <- round(rawR) # ensures we are dealing with a whole number of recruits.
  return(roundR)
}

# 4.  This function converts size-specific fecundity into a proportional recruit contribution (this function completes the IPM loop between 
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

# 5. Recruit Size
# the pdf of size z1 recruits next summer - 
c0_z1 <- function(Size.t1, m.par){
  mean <- m.par["rec.size"]
  sd <- m.par["rec.size.sd"]
  p_den_rcsz <- dnorm(Size.t1, mean = mean, sd = sd) 
  return(p_den_rcsz)
}

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

# Now define the number of bins 
# The number of bins used needs to be slightly higher than optimum to account for the extended size range for integration.
# And will be standardized across all models  
m <- 200 

### Implement IPM models and estimate Net reproductive rate ------------------------------------
# Generate output storage.
C.R0 <- matrix(NA, nrow = 1000, ncol = 4)
ST.R0 <- matrix(NA, nrow = 1000, ncol = 4)
W.R0 <- matrix(NA, nrow = 1000, ncol = 4)
colnames(C.R0) <- colnames(ST.R0) <- colnames(W.R0) <- c("AS R0", "AT R0", "JS R0", "JT R0")

### 1. Construct each IPM and display lambda to check the models have worked
# Competitive
CAS.IPM <- mk_GF(m = m, m.par = unlist(m.par.store[1,]), L = L.store[1], U = U.store[1])
CAS.R0 <- net_repro_rate(CAS.IPM$G, CAS.IPM$F); CAS.R0; C.R0[1,1] <- CAS.R0

CAT.IPM <- mk_GF(m = m, m.par = unlist(m.par.store[2,]), L = L.store[1], U = U.store[1])
CAT.R0 <- net_repro_rate(CAT.IPM$G, CAT.IPM$F); CAT.R0; C.R0[1,2] <- CAT.R0

CJS.IPM <- mk_GF(m = m, m.par = unlist(m.par.store[3,]), L = L.store[1], U = U.store[1])
CJS.R0 <- net_repro_rate(CJS.IPM$G, CJS.IPM$F); CJS.R0; C.R0[1,3] <- CJS.R0

CJT.IPM <- mk_GF(m = m, m.par = unlist(m.par.store[4,]), L = L.store[1], U = U.store[1])
CJT.R0 <- net_repro_rate(CJT.IPM$G, CJT.IPM$F); CJT.R0; C.R0[1,4] <- CJT.R0

# Stress Tolerant
STAS.IPM <- mk_GF(m = m, m.par = unlist(m.par.store[5,]), L = L.store[2], U = U.store[2])
STAT.IPM <- mk_GF(m = m, m.par = unlist(m.par.store[6,]), L = L.store[2], U = U.store[2])
STJS.IPM <- mk_GF(m = m, m.par = unlist(m.par.store[7,]), L = L.store[2], U = U.store[2])
STJT.IPM <- mk_GF(m = m, m.par = unlist(m.par.store[8,]), L = L.store[2], U = U.store[2])

STAS.R0<- net_repro_rate(STAS.IPM$G, STAS.IPM$F); STAS.R0; ST.R0[1,1] <- STAS.R0
STAT.R0<- net_repro_rate(STAT.IPM$G, STAT.IPM$F); STAT.R0; ST.R0[1,2] <- STAT.R0
STJS.R0<- net_repro_rate(STJS.IPM$G, STJS.IPM$F); STJS.R0; ST.R0[1,3] <- STJS.R0
STJT.R0<- net_repro_rate(STJT.IPM$G, STJT.IPM$F); STJT.R0; ST.R0[1,4] <- STJT.R0

# Weedy
WAS.IPM <- mk_GF(m = m, m.par = unlist(m.par.store[9,]), L = L.store[3], U = U.store[3])
WAT.IPM <- mk_GF(m = m, m.par = unlist(m.par.store[10,]), L = L.store[3], U = U.store[3])
WJS.IPM <- mk_GF(m = m, m.par = unlist(m.par.store[11,]), L = L.store[3], U = U.store[3])
WJT.IPM <- mk_GF(m = m, m.par = unlist(m.par.store[12,]), L = L.store[3], U = U.store[3])

WAS.R0 <- net_repro_rate(WAS.IPM$G, WAS.IPM$F); WAS.R0; W.R0[1,1] <- WAS.R0
WAT.R0 <- net_repro_rate(WAT.IPM$G, WAT.IPM$F); WAT.R0; W.R0[1,2] <- WAT.R0
WJS.R0 <- net_repro_rate(WJS.IPM$G, WJS.IPM$F); WJS.R0; W.R0[1,3] <- WJS.R0
WJT.R0 <- net_repro_rate(WJT.IPM$G, WJT.IPM$F); WJT.R0; W.R0[1,4] <- WJT.R0

#######################################
# STEP 5: Run a Jack-knife to estimate variance in Generation time
#######################################

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
  
  ## 3. Recruitment --------
  # Identify recruits.
  try(jrecruit.data <- data.sample[which(data.sample$Recruit.t1 == "Yes"),])
  # A. Number of recruits
  # (This will be maintained from original data set and just allowed to vary within observed limits)
  # B. Recruit size
  try(jrecruit.size.mod2 <- lm(Size.t1 ~ LHS, data = jrecruit.data))
  
  ## 4. Fecundity ----------
  try(jfec.mod2 <- nls(Fecundity ~ a * exp(b * Size), data = fec.sample[which(fec.sample$LHS == "C"),], start = list(a = 1, b = 1)))
  try(jfec.mod3 <- nls(Fecundity ~ a * exp(b * Size), data = fec.sample[which(fec.sample$LHS == "ST"),], start = list(a = 1, b = 1)))
  try(jfec.mod4 <- nls(Fecundity ~ a * exp(b * Size), data = fec.sample[which(fec.sample$LHS == "W"),], start = list(a = 1, b = 1)))
  
  ##### Store vital rate parameters --------------------------------------------
  
  # Extract necessary mean coefficients.
  try(jsurvival.coefs <- apply(coef(jsurv.mod)$Colony.ID, 2, mean))
  try(jgrowth.coefs <- apply(coef(jgrow.mod)$Colony.ID, 2, mean))
  try(jgrowth.SD.coefs <- apply(coef(jgrowth.sd.mod)$'jgrow.mod@frame$Colony.ID', 2, mean))
  try(jgrowth.SD.coefs <- jgrowth.SD.coefs[2:25])
   
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
    frag.int               =  m.par.store[1, 8],
    frag.slope             =  m.par.store[1, 9],
    frag.slope2            =  m.par.store[1, 10],
    # Number of frags produced
    no.frag.int            =  m.par.store[1, 11],
    no.frag.slope          =  m.par.store[1, 12],
    # Fragment size
    frag.s.int             =  m.par.store[1, 13],
    frag.s.slope           =  m.par.store[1, 14],
    frag.s.slope2          =  m.par.store[1, 15],
    # Variance in fragment size
    frag.sd.int            =  m.par.store[1, 16],
    frag.sd.slope          =  m.par.store[1, 17],
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
    frag.int               =  m.par.store[2, 8],
    frag.slope             =  m.par.store[2, 9],
    frag.slope2            =  m.par.store[2, 10],
    # Number of frags produced
    no.frag.int            =  m.par.store[2, 11],
    no.frag.slope          =  m.par.store[2, 12],
    # Fragment size
    frag.s.int             =  m.par.store[2, 13],
    frag.s.slope           =  m.par.store[2, 14],
    frag.s.slope2          =  m.par.store[2, 15],
    # Variance in fragment size
    frag.sd.int            =  m.par.store[2, 16],
    frag.sd.slope          =  m.par.store[2, 17],
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
    frag.int               =  m.par.store[3, 8],
    frag.slope             =  m.par.store[3, 9],
    frag.slope2            =  m.par.store[3, 10],
    # Number of frags produced
    no.frag.int            =  m.par.store[3, 11],
    no.frag.slope          =  m.par.store[3, 12],
    # Fragment size
    frag.s.int             =  m.par.store[3, 13],
    frag.s.slope           =  m.par.store[3, 14],
    frag.s.slope2          =  m.par.store[3, 15],
    # Variance in fragment size
    frag.sd.int            =  m.par.store[3, 16],
    frag.sd.slope          =  m.par.store[3, 17],
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
    frag.int               =  m.par.store[4, 8],
    frag.slope             =  m.par.store[4, 9],
    frag.slope2            =  m.par.store[4, 10],
    # Number of frags produced
    no.frag.int            =  m.par.store[4, 11],
    no.frag.slope          =  m.par.store[4, 12],
    # Fragment size
    frag.s.int             =  m.par.store[4, 13],
    frag.s.slope           =  m.par.store[4, 14],
    frag.s.slope2          =  m.par.store[4, 15],
    # Variance in fragment size
    frag.sd.int            =  m.par.store[4, 16],
    frag.sd.slope          =  m.par.store[4, 17],
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
    frag.int               =  m.par.store[5, 8],
    frag.slope             =  m.par.store[5, 9],
    frag.slope2            =  m.par.store[5, 10],
    # Number of frags produced
    no.frag.int            =  m.par.store[5, 11],
    no.frag.slope          =  m.par.store[5, 12],
    # Fragment size
    frag.s.int             =  m.par.store[5, 13],
    frag.s.slope           =  m.par.store[5, 14],
    frag.s.slope2          =  m.par.store[5, 15],
    # Variance in fragment size
    frag.sd.int            =  m.par.store[5, 16],
    frag.sd.slope          =  m.par.store[5, 17],
    # Colony Fecundity
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
    frag.int               =  m.par.store[6, 8],
    frag.slope             =  m.par.store[6, 9],
    frag.slope2            =  m.par.store[6, 10],
    # Number of frags produced
    no.frag.int            =  m.par.store[6, 11],
    no.frag.slope          =  m.par.store[6, 12],
    # Fragment size
    frag.s.int             =  m.par.store[6, 13],
    frag.s.slope           =  m.par.store[6, 14],
    frag.s.slope2          =  m.par.store[6, 15],
    # Variance in fragment size
    frag.sd.int            =  m.par.store[6, 16],
    frag.sd.slope          =  m.par.store[6, 17],
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
    frag.int               =  m.par.store[7, 8],
    frag.slope             =  m.par.store[7, 9],
    frag.slope2            =  m.par.store[7, 10],
    # Number of frags produced
    no.frag.int            =  m.par.store[7, 11],
    no.frag.slope          =  m.par.store[7, 12],
    # Fragment size
    frag.s.int             =  m.par.store[7, 13],
    frag.s.slope           =  m.par.store[7, 14],
    frag.s.slope2          =  m.par.store[7, 15],
    # Variance in fragment size
    frag.sd.int            =  m.par.store[7, 16],
    frag.sd.slope          =  m.par.store[7, 17],
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
    frag.int               =  m.par.store[8, 8],
    frag.slope             =  m.par.store[8, 9],
    frag.slope2            =  m.par.store[8, 10],
    # Number of frags produced
    no.frag.int            =  m.par.store[8, 11],
    no.frag.slope          =  m.par.store[8, 12],
    # Fragment size
    frag.s.int             =  m.par.store[8, 13],
    frag.s.slope           =  m.par.store[8, 14],
    frag.s.slope2          =  m.par.store[8, 15],
    # Variance in fragment size
    frag.sd.int            =  m.par.store[8, 16],
    frag.sd.slope          =  m.par.store[8, 17],
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
    frag.int               =  m.par.store[9, 8],
    frag.slope             =  m.par.store[9, 9],
    frag.slope2            =  m.par.store[9, 10],
    # Number of frags produced
    no.frag.int            =  m.par.store[9, 11],
    no.frag.slope          =  m.par.store[9, 12],
    # Fragment size
    frag.s.int             =  m.par.store[9, 13],
    frag.s.slope           =  m.par.store[9, 14],
    frag.s.slope2          =  m.par.store[9, 15],
    # Variance in fragment size
    frag.sd.int            =  m.par.store[9, 16],
    frag.sd.slope          =  m.par.store[9, 17],
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
    frag.int               =  m.par.store[10, 8],
    frag.slope             =  m.par.store[10, 9],
    frag.slope2            =  m.par.store[10, 10],
    # Number of frags produced
    no.frag.int            =  m.par.store[10, 11],
    no.frag.slope          =  m.par.store[10, 12],
    # Fragment size
    frag.s.int             =  m.par.store[10, 13],
    frag.s.slope           =  m.par.store[10, 14],
    frag.s.slope2          =  m.par.store[10, 15],
    # Variance in fragment size
    frag.sd.int            =  m.par.store[10, 16],
    frag.sd.slope          =  m.par.store[10, 17],
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
    frag.int               =  m.par.store[11, 8],
    frag.slope             =  m.par.store[11, 9],
    frag.slope2            =  m.par.store[11, 10],
    # Number of frags produced
    no.frag.int            =  m.par.store[11, 11],
    no.frag.slope          =  m.par.store[11, 12],
    # Fragment size
    frag.s.int             =  m.par.store[11, 13],
    frag.s.slope           =  m.par.store[11, 14],
    frag.s.slope2          =  m.par.store[11, 15],
    # Variance in fragment size
    frag.sd.int            =  m.par.store[11, 16],
    frag.sd.slope          =  m.par.store[11, 17],
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
    frag.int               =  m.par.store[12, 8],
    frag.slope             =  m.par.store[12, 9],
    frag.slope2            =  m.par.store[12, 10],
    # Number of frags produced
    no.frag.int            =  m.par.store[12, 11],
    no.frag.slope          =  m.par.store[12, 12],
    # Fragment size
    frag.s.int             =  m.par.store[12, 13],
    frag.s.slope           =  m.par.store[12, 14],
    frag.s.slope2          =  m.par.store[12, 15],
    # Variance in fragment size
    frag.sd.int            =  m.par.store[12, 16],
    frag.sd.slope          =  m.par.store[12, 17],
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
  try(jCAS.IPM <- mk_GF(m = m, m.par = unlist(jm.par[1,]), L = jL.store[1], U = jU.store[1]))
  try(jCAT.IPM <- mk_GF(m = m, m.par = unlist(jm.par[2,]), L = jL.store[1], U = jU.store[1]))
  try(jCJS.IPM <- mk_GF(m = m, m.par = unlist(jm.par[3,]), L = jL.store[1], U = jU.store[1]))
  try(jCJT.IPM <- mk_GF(m = m, m.par = unlist(jm.par[4,]), L = jL.store[1], U = jU.store[1]))
  # Stress Tolerant
  try(jSTAS.IPM <- mk_GF(m = m, m.par = unlist(jm.par[5,]), L = jL.store[2], U = jU.store[2]))
  try(jSTAT.IPM <- mk_GF(m = m, m.par = unlist(jm.par[6,]), L = jL.store[2], U = jU.store[2]))
  try(jSTJS.IPM <- mk_GF(m = m, m.par = unlist(jm.par[7,]), L = jL.store[2], U = jU.store[2]))
  try(jSTJT.IPM <- mk_GF(m = m, m.par = unlist(jm.par[8,]), L = jL.store[2], U = jU.store[2]))
  # Weedy
  try(jWAS.IPM <- mk_GF(m = m, m.par = unlist(jm.par[9,]), L = jL.store[3], U = jU.store[3]))
  try(jWAT.IPM <- mk_GF(m = m, m.par = unlist(jm.par[10,]), L = jL.store[3], U = jU.store[3]))
  try(jWJS.IPM <- mk_GF(m = m, m.par = unlist(jm.par[11,]), L = jL.store[3], U = jU.store[3]))
  try(jWJT.IPM <- mk_GF(m = m, m.par = unlist(jm.par[12,]), L = jL.store[3], U = jU.store[3]))
  
  ##### Store Generation time ----------------------------------------------------
  # Competitive
  try(C.R0[ii,1] <- net_repro_rate(jCAS.IPM$G, jCAS.IPM$F))
  try(C.R0[ii,2] <- net_repro_rate(jCAT.IPM$G, jCAT.IPM$F))
  try(C.R0[ii,3] <- net_repro_rate(jCJS.IPM$G, jCJS.IPM$F))
  try(C.R0[ii,4] <- net_repro_rate(jCJT.IPM$G, jCJT.IPM$F))
  # Stress Tolerant
  try(ST.R0[ii,1] <- net_repro_rate(jSTAS.IPM$G, jSTAS.IPM$F))
  try(ST.R0[ii,2] <- net_repro_rate(jSTAT.IPM$G, jSTAT.IPM$F))
  try(ST.R0[ii,3] <- net_repro_rate(jSTJS.IPM$G, jSTJS.IPM$F))
  try(ST.R0[ii,4] <- net_repro_rate(jSTJT.IPM$G, jSTJT.IPM$F))
  # Weedy
  try(W.R0[ii,1] <- net_repro_rate(jWAS.IPM$G, jWAS.IPM$F))
  try(W.R0[ii,2] <- net_repro_rate(jWAT.IPM$G, jWAT.IPM$F))
  try(W.R0[ii,3] <- net_repro_rate(jWJS.IPM$G, jWJS.IPM$F))
  try(W.R0[ii,4] <- net_repro_rate(jWJT.IPM$G, jWJT.IPM$F))
  
  # And tidy up
  try(rm(data.sample, fec.sample, jsurv.mod, jgrow.mod, jgrowth.sd.mod, jfec.mod2, jfec.mod3, jfec.mod4, jrecruit.size.mod2, jrecruit.data,
         jsurvival.coefs, jgrowth.coefs, jgrowth.SD.coefs, jm.par,
         jcoefs.C.A.S, jcoefs.C.A.T, jcoefs.C.J.S, jcoefs.C.J.T, jcoefs.ST.A.S, jcoefs.ST.A.T, jcoefs.ST.J.S, jcoefs.ST.J.T, jcoefs.W.A.S, jcoefs.W.A.T, jcoefs.W.J.S, jcoefs.W.J.T,
         jC.L, jC.U, jST.L, jST.U, jW.L, jW.U, jL.store, jU.store,
         jCAS.IPM, jCAT.IPM, jCJS.IPM, jCJT.IPM, jSTAS.IPM, jSTAT.IPM, jSTJS.IPM, jSTJT.IPM, jWAS.IPM, jWAT.IPM, jWJS.IPM, jWJT.IPM))
}

# Save outputs as a checkpoint
write.csv(C.R0, file = "C_R0.csv", row.names = FALSE)
write.csv(ST.R0, file = "ST_R0.csv", row.names = FALSE)
write.csv(W.R0, file = "W_R0.csv", row.names = FALSE)

#reset working directory
setwd("C:/Users/rhian/OneDrive/Desktop/James/Chapter 4/Data")

## Read in the stored R0 data
C_R0 <- read.csv("C_R0.csv")
ST_R0 <- read.csv("ST_R0.csv")
W_R0 <- read.csv("W_R0.csv")

# and attached to the main dataframe being worked on here
GT_data$net_repro_rate <-  c(C_R0[,"AS.R0"], C_R0[,"AT.R0"], C_R0[,"JS.R0"], C_R0[,"JT.R0"],
                             ST_R0[,"AS.R0"], ST_R0[,"AT.R0"], ST_R0[,"JS.R0"], ST_R0[,"JT.R0"],
                             W_R0[,"AS.R0"], W_R0[,"AT.R0"], W_R0[,"JS.R0"], W_R0[,"JT.R0"])
# Because I used set seed I know these values will align with those already generated during earlier attempts

# now calculate the generation time of each population
GT_data$GT <- log(GT_data$net_repro_rate) / log(GT_data$Lambda)
hist(GT_data$GT)
## and tidy up generation time data by removing negative entries.
GT_data$GT[which(GT_data$GT < 0)] <- NA
# Now remove any outliers
# Function taken from Robs code.
removeOutliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.025, .975), na.rm = na.rm) #this determines the 25% and 75% limits of the data set based on the supplied vector (ie. where 95% of the data lies) 
  H <- 1.5 * IQR(x, na.rm = na.rm) # determines the interquartile range (variance of the dataset)
  y <- x
  y[x < (qnt[1] - H)] <- NA #any values below or above these ranges are then excluded.
  y[x > (qnt[2] + H)] <- NA
  y
}
# loop through and remove outliers
GT_data$GT <- removeOutliers(GT_data$GT, na.rm = T)
# re check distribution
hist(GT_data$GT)
hist(log(GT_data$GT))
# and log transform
GT_data$logGT <- log(GT_data$GT)

# save data as a checkpoint
write.csv(GT_data, file = "Generation time outputs.csv", row.names = FALSE)

# Little tidy up (remove weird lambdas)
GT_data[which(GT_data$Lambda > 2 | GT_data$Lambda < 0), c("GT","Lambda","TE")] <- NA
# check distributions
hist(GT_data$Lambda)
hist(GT_data$TE) # needs fixing
BoxCoxTrans(GT_data$TE, na.rm = TRUE)
GT_data$tTE <- 1/(GT_data$TE^0.1) # this is now inverted.
hist(GT_data$logGT)
# and remove NAs
GT_data <- na.omit(GT_data)
# Estimate population means
GT_data_mean <- aggregate(GT_data, list(GT_data$Country, GT_data$Region, GT_data$LHS), mean, na.rm = T)[c(1:3,7:12)]
colnames(GT_data_mean) <- c("Country", "Region", "LHS", "GT", "Lambda","TE","R0","LogGT","tTE")
# and population CI
GT_data_sd <- aggregate(GT_data, list(GT_data$Country, GT_data$Region, GT_data$LHS), CI)[,c(1:3,7:12)]
GT_sd <- as.data.frame(GT_data_sd[,8,])
L_sd <- as.data.frame(GT_data_sd[,5,])
TE_sd <- as.data.frame(GT_data_sd[,9,])
# add these to the mean files
GT_data_mean <- do.call(cbind, list(GT_data_mean, GT_sd, L_sd, TE_sd))
colnames(GT_data_mean) <- c("Country", "Region", "LHS", "GT", "Lambda", "TE", "R0", "LogGT", "tTE",
                            "GT.upper", "GT.mean", "GT.lower",
                            "L.upper", "L.mean", "L.lower",
                            "TE.upper", "TE.mean", "TE.lower")

# Run a linear regression on the full data (Type 2 regression - cannot handle interaction terms)
# first test the variance magnitude of the variables of interest to determine the type of type 2 linear regression to use.
var(GT_data$logGT);var(GT_data$tTE);var(GT_data$Lambda) # the variance is not equal but larger in x than y (TE and Lambda), the data is normal and dimensionless.
# so a ranged major axis regression will be used.
# to explore relationship between lambda
l_mod <- lmodel2(Lambda ~ logGT, data = GT_data, "interval", "interval")
# and transient potential
te_mod <- lmodel2(tTE ~ logGT, data = GT_data, "interval", "interval")
# and extract data for plotting regression lines
GT_data$predict.L.mean <- l_mod$regression.results[4,2] + (GT_data$logGT*l_mod$regression.results[4,3])
GT_data$predict.L.lower <- l_mod$confidence.intervals[4,2] + (GT_data$logGT*l_mod$confidence.intervals[4,4])
GT_data$predict.L.upper <- l_mod$confidence.intervals[4,3] + (GT_data$logGT*l_mod$confidence.intervals[4,5])
GT_data$predict.TE.mean <- te_mod$regression.results[4,2] + (GT_data$logGT*te_mod$regression.results[4,3])
GT_data$predict.TE.lower <- te_mod$confidence.intervals[4,2] + (GT_data$logGT*te_mod$confidence.intervals[4,4])
GT_data$predict.TE.upper <- te_mod$confidence.intervals[4,3] + (GT_data$logGT*te_mod$confidence.intervals[4,5])

# Plot lambda ~ generation time
ggplot(GT_data_mean, aes(x = LogGT, y = Lambda)) +
  geom_ribbon(aes(ymin = predict.L.lower, ymax = predict.L.upper, x = logGT), fill = "gray", alpha = 0.5, data = GT_data) +
  geom_errorbar(aes(ymin = L.lower, ymax = L.upper, width = 0.05), size = 0.9) +
  geom_errorbarh(aes(xmin = GT.lower, xmax = GT.upper, height = 0.012), size = 0.9) +
  geom_point(aes(col = LHS, shape = Region), size = 6) +
  geom_line(aes(y = predict.L.mean, x = logGT), col = "Black", linetype = "dashed", size = 1.1, data = GT_data) +
  # Formatting
  scale_x_continuous(expand = expansion(mult = c(0.0,0.0)), limits = c(2.3,6.5)) +
  scale_y_continuous(expand = expansion(mult = c(0.05,0.05)), limits = c(0.5,1.2)) +
  scale_color_manual(name = "", values = c("blue2", "orange2", "red2"), 
                     labels = c("Competitive","Stress-Tolerant","Weedy"),
                     guide = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 10))) +
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

# plot transient potential ~ generation time
ggplot(GT_data_mean, aes(x = LogGT, y = -tTE)) +
  geom_ribbon(aes(ymin = -predict.TE.lower, ymax = -predict.TE.upper, x = logGT), fill = "gray", alpha = 0.5, data = GT_data) +
  geom_errorbar(aes(ymin = -TE.lower, ymax = -TE.upper, width = 0.05), size = 0.9) +
  geom_errorbarh(aes(xmin = GT.lower, xmax = GT.upper, height = 0.012), size = 0.9) +
  geom_point(aes(col = LHS, shape = Region), size = 6) +
  geom_line(aes(y = -predict.TE.mean, x = logGT), col = "Black", linetype = "dashed", size = 1.1, data = GT_data) +
  # Formatting
  scale_x_continuous(expand = expansion(mult = c(0.0,0.0)), limits = c(2.3,6.5)) +
  scale_y_continuous(expand = expansion(mult = c(0.05,0.05)), limits = c(-1.2,-0.6)) +
  scale_color_manual(name = "", values = c("blue2", "orange2", "red2"), 
                     labels = c("Competitive","Stress-Tolerant","Weedy"),
                     guide = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 10))) +
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

################################################################ End code #################################################################
