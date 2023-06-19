# This script is calculating the vital rate sensitivities of the transient dynamics of tropical and subtropical croa assemblages.
# This script relies on outputs from the 'Constructing IPMs' script.
# Due to the size of the sensitivity matrices that are to be generated using this script, the process in carried out iteratively will later scipts used to condense this information back together as needed.

# Author: James Cant (jic2@st-andrews.ac.uk)
# ----------------------------------------------------------------------------------

#load/install required packages
if (!require("popdemo")) {
  install.packages("popdemo")
}

if (!require("dplyr")) {
  install.packages("dplyr")
}

# set working directory for sourcing output files
# setwd("data_file_location")

# load output files if needed
C.outputs <- read.csv("Competitive outputs.csv")
ST.outputs <- read.csv("Stress Tolerant outputs.csv")
W.outputs <- read.csv("Weedy outputs.csv")
# arrange as a list for the loop below
outputList <- list(C.outputs, ST.outputs, W.outputs)
# and tidy up
rm(C.outputs, ST.outputs, W.outputs)

############################################
# STEP 1: Define sensitivity and elasticity function
############################################

# 1. Sensitivity & Elasticity of the transient envelope - brute force method (Morris & Doak 2002)
SensTE <- function(i, matK, MAmp, j){ 
  
  # create blank dataframe
  Sens.out <- data.frame(Sens = NA,
                         Elas = NA,
                         SensMean = NA,
                         ElasMean = NA)
  
  # Brute force perturbation analysis
  matNew <- matK # re-save the original matrix to avoid changing it accidentally
  perturb <- 0.01 # makes it simple to change the size of the perturbation
  elementdiff <- matNew[i,j]*perturb
  matNew[i,j] <- matNew[i,j] + (elementdiff) # make a small perturbation to the element (the multiplication ensures that zero values remain unchanged and therefore won't return a sensitivity value)
  try(MAmp.new <- max.amp(matNew)) # estimate new maximal amplification
  try(SensK <- (MAmp.new - MAmp)/elementdiff) # using the change in the envelope and the change in element calculate the element level sensitivity - this follows equation 9.4 (Morris & Doak 2002 pg332)
  if(SensK %in% c(Inf, NA)){SensK = 0}
  try(ElasK <- ((MAmp.new - MAmp)/MAmp)/(elementdiff/matK[i,j])) # and calculate the corresponding elasticity values
  
  # Store outputs
  Sens.out <- Sens.out %>% mutate(Sens = SensK,
                                  Elas = ElasK,
                                  SensMean = NA, # just a little place holder to keep dimensions matching
                                  ElasMean = NA)
  
  # progress read out 
  print(i)
  
  # and return the sensitivity and elasticity matrices for the transient envelope of matK.
  return(Sens.out)
}

# 2. Sensitivity & Elasticity of the transient envelope - include estimation of mean matrix sensitivities
SensTEMean <- function(i, matK, MAmp, j, matK2){ 
  
  # First estimate the mean matrix and its associated transient bounds 
  matMean <- (matK + matK2)/2 # calculate mean matrix
  #estimate the kreiss bounds for the mean matrix
  try(MAmp2 <- max.amp(matMean)) # estimate new maximal amplification for the mean matrix
  
  # create blank dataframe
  Sens.out <- data.frame(Sens = NA,
                         Elas = NA,
                         SensMean = NA,
                         ElasMean = NA)
  
  # Brute force perturbation analysis of the original matrix
  matNew <- matK # re-save the original matrix to avoid changing it accidentally
  perturb <- 0.01 # makes it simple to change the size of the perturbation
  elementdiff <- matNew[i,j]*perturb
  matNew[i,j] <- matNew[i,j] + (elementdiff) # make a small perturbation to the element (the multiplication ensures that zero values remain unchanged and therefore won't return a sensitivity value)
  try(MAmp.new <- max.amp(matNew)) # estimate new maximal amplification
  try(SensK <- (MAmp.new - MAmp)/elementdiff) # using the change in the envelope and the change in element calculate the element level sensitivity - this follows equation 9.4 (Morris & Doak 2002 pg332)
  if(SensK %in% c(Inf, NA)){SensK = 0}
  try(ElasK <- ((MAmp.new - MAmp)/MAmp)/(elementdiff/matK[i,j])) # and calculate the corresponding elasticity values
  
  #Now estimate the sensitivity and elasticity of the mean matrix
  #reassign mean matrix
  matMeanNew <- matMean
  elementdiff2 <- matMeanNew[i,j]*perturb
  matMeanNew[i,j] <- matMeanNew[i,j] + (elementdiff2) # make a small perturbation to the element (the multiplication ensures that zero values remain unchanged and therefore won't return a sensitivity value)n
  try(MAmp2.new <- max.amp(matMeanNew)) # estimate new maximal amplification
  try(SensK2 <- (MAmp2.new - MAmp2)/elementdiff2) # using the change in the envelope and the change in element calculate the element level sensitivity - this follows equation 9.4 (Morris & Doak 2002 pg332)
  if(SensK2 %in% c(Inf, NA)){SensK2 = 0}
  try(ElasK2 <- ((MAmp2.new - MAmp2)/MAmp2)/(elementdiff2/matMean[i,j])) # and calculate the corresponding elasticity values
  
  # Store outputs
  Sens.out <- Sens.out %>% mutate(Sens = SensK, # original matrix sensitivities
                                  Elas = ElasK,
                                  SensMean = SensK2, # mean matrix sensitivities
                                  ElasMean = ElasK2)
  
  # progress read out 
  print(i)
  
  # and return the sensitivity and elasticity matrices for the transient envelope of matK.
  return(Sens.out)
}

# 3. Modified function to speed up analysis of maximal amplification indices
max.amp <-
  function(A){
    # check matrix dimensions
    if(any(length(dim(A)) != 2, dim(A)[1] != dim(A)[2])) stop("A must be a square matrix")
    
    # Calculate key parameters and produce storage outputs
    order <- dim(A)[1]
    maxN<-numeric(order)
    times<-numeric(order)
    M <- A
    # estimate lambda
    eigvals <- eigen(M)$values
    lmax <- which.max(Re(eigvals))
    lambda <- Re(eigvals[lmax])
    
    # remove the asymptotic characteristics of the model
    A <- M/lambda
    
    # Estimate maximal amplification
    # what is the time required for convergence?
    # this is a slow function so from previous run throughs that this converges at 45. So I will set that manually.
    maxtime <- 45
    # project the population to convergence
    projection <- project(A, time = maxtime)
    for(ii in 1:order){
      maxN[ii] <- max(projection[,ii])
      times[ii] <- which.max(projection[,ii])-1
    }
    # extract desired value
    rhomax <- max(maxN)
    t <- times[which.max(maxN)]
    
    # return maximal amplification value
    return(maxamp = rhomax)
  }

############################################
# STEP 2: Extract initial coral models and estimate element level sensitivities and elasticities
############################################

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

# Define indexing vectors
j <- 126 # what row is being dealt with
output_index <- rep(1:3, each = 4)
L_index <- rep(1:4, length.out = 12)
MAmp_index <- rep(29:32, length.out = 12)

# create output listS
sensList <- list()
elasList <- list()
sensMeanList <- list()
elasMeanList <- list()

# Run sensitivity loop
for(ii in 1:12) { # there is twelve models for this to be repeated on
  # progress read out
  cat("\n",ii,"\n")
  
  # extract lambda
  Lambda <- outputList[[output_index[ii]]][1,L_index[ii]]
  # Extract transient envelope  
  MAmp <- outputList[[output_index[ii]]][1,MAmp_index[ii]]
  # Select Kernel and standardize to distinguish transient characteristics 
  matK <- ListK[[ii]]/Lambda
  # Extract corresponding tropical kernel & standardize
  matK2 <- ListKtrop[[ii]]
  # and standardize also
  Lambda2 <- as.numeric(eigen(matK2)$values[1])
  matK2 <- matK2/Lambda2
  
  # define the dimensions of the model
  drow <- dim(matK)[1] #number of rows
  
  # run sensitivity analysis
  if(ii %in% c(1,3,5,7,9,11)) { 
    try(Sens.out <- lapply(1:drow, SensTEMean, matK = matK, MAmp = MAmp, j = j, matK2 = matK2))
  } else {
    try(Sens.out <- lapply(1:drow, SensTE, matK = matK, MAmp = MAmp, j = j))
  }
  # reformat and extract the components from the returned list
  SensK <- unlist(lapply(Sens.out, '[[', 1)) # this is a column vector of element level sensitivities
  ElasK <- unlist(lapply(Sens.out, '[[', 2)) # and this is elasticities
  SensK2 <- unlist(lapply(Sens.out, '[[', 3)) # and repeat for mean sensitivities
  ElasK2 <- unlist(lapply(Sens.out, '[[', 4)) 
  # tidy up output
  SensK[is.nan(SensK)] = 0
  ElasK[is.nan(ElasK)] = 0
  SensK2[is.nan(SensK2)] = 0
  ElasK2[is.nan(ElasK2)] = 0
  
  # and save
  sensList[[ii]] <- SensK
  elasList[[ii]] <- ElasK
  sensMeanList[[ii]] <- SensK2
  elasMeanList[[ii]] <- ElasK2
  
  # and tidy up memory
  rm(Lambda, Lambda2, MAmp, matK, matK2, Sens.out, ElasK, SensK, SensK2, ElasK2)
}

############################################
# STEP 3: Reformat and save outputs
############################################
# Unlist outputs 
output_file <- data.frame(Sens = unlist(sensList),
                          Elas = unlist(elasList),
                          SensMean = unlist(sensMeanList),
                          ElasMean = unlist(elasMeanList))
# and save
write.csv(output_file, file = paste0("Sensitivity matrices_col-",j,".csv", sep = ""), row.names = FALSE)

###################################################################################### End of Code ###################################


