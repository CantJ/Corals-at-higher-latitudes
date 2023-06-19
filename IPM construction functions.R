# Functions for creating IPM models
# -----------------------------------------------------------

# 1a. Define function for discretising and combining the survival + growth (G) fragmentation (H) (G + H = P) and recruitment (F) kernels
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

# 1b. Define mK_P function
mk_P <- function(m, m.par, L, U) {
  # define the mesh points 
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  # Define the survival and growth kernel
  G <- h * (outer(meshpts, meshpts, G_z1z, m.par))
  # Define the fragmentation kernel
  H <- h * (outer(meshpts, meshpts, H_z1z, m.par))
  P <- G + H # define adult only kernel (P and H)
  return(list(P = P, meshpts = meshpts, G = G, H = H))
}


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

