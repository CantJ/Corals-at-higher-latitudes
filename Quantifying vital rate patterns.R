# This script is for estimating the vital rates characteristics of corals from 
# the three life-history classifications, Weedy, Competitive, and Stress tolerant
# and how they vary between tropical and subtropical regions.
# These characteristics will then be used to construct IPM models evaluating how these characteristics influence bio geographic trends in coral assemblage performance
# across Australia and Japan.

# Author: James Cant (jic2@st-andrews.ac.uk)
# Date last modified: April 2020
# -----------------------------------------------------------------------------------------

# clear workspace
rm(list=ls(all=TRUE))

# set working directory
setwd("Data_file_location")

# load relevant packages.
library(dplyr)
library(lme4)
library(ggeffects)
library(car)
library(ggplot2)
library(Rmisc)
library(tidyr)
library(ggstance)

# read in the full demographic dataset
demographic.data <- read.csv("Raw Coral Demographic data.csv", stringsAsFactors = TRUE)
# and read in the fecundity dataset
# This dataset is accessible from the coral trait database.
fec.data <- read.csv("Fecundity data - Hall & Hughes 1996.csv", stringsAsFactors = TRUE)

# subset demographic datafile to focus on only the relevant populations.
demographic.data <- demographic.data[which(demographic.data$Group == "Scleractinia" & 
                                           demographic.data$LHS != "UA"),] # remove unassigned life histories

# ensure data is from desired time period
demographic.data <- demographic.data[which(demographic.data$time.t %in% c("2016", "2017", "2018")),]

###############################################
# STEP 1: Format variables, and estimate survival and fragmentation
###############################################

# remove entries corresponding with lost corals. These can't be used in estimating demographic parameters.
demographic.data <- subset(demographic.data, Size.t != "Lost" | is.na(Size.t)) #<- this removes all rows containing "Lost" in the size t section, but importantly retains any NA entries - these are really important!
demographic.data <- subset(demographic.data, Size.t1 != "Lost" | is.na(Size.t1)) #<- this does the same for Size at t+1

# convert variables from the dataframe into their required formats
demographic.data$Size.t <- as.numeric(paste(demographic.data$Size.t))
demographic.data$Size.t1 <- as.numeric(paste(demographic.data$Size.t1))
demographic.data$Tag.year <- as.factor(demographic.data$Tag.year)
demographic.data$No.frag.t1 <- as.numeric(paste(demographic.data$No.frag.t1))
demographic.data$LHS <- droplevels(demographic.data$LHS, exclude = "UA") # just to remove the Un-assigned category which has been dropped from the data.

# covert size data into cm^2
demographic.data$Size.t <- demographic.data$Size.t/100; demographic.data$Size.t1 <- demographic.data$Size.t1/100 

# check the distribution of size data and transform to achieve normality if necessary
hist(demographic.data$Size.t)
hist(demographic.data$Size.t1)
hist(log(demographic.data$Size.t))
hist(log(demographic.data$Size.t1)) #this have a better distribution
# transform size data
demographic.data$Size.t <- log(demographic.data$Size.t)
demographic.data$Size.t1 <- log(demographic.data$Size.t1) 

# estimate colony survival
demographic.data$surv <- NA
demographic.data$surv[which(!is.na(demographic.data$Size.t) & !is.na(demographic.data$Size.t1))] = 1 #<- when both Size.t0 and Size.t1 are both NOT NA then surv = 1
demographic.data$surv[which(!is.na(demographic.data$Size.t) & is.na(demographic.data$Size.t1))] = 0 #<- when Size.t is NOT NA but Size.t1 is then surv = 0
table(demographic.data$surv)
# estimate colony fragmentation
demographic.data$frag <- NA
demographic.data$frag[which(demographic.data$No.frag.t1 > 0)] = 1 #<- when a colony produced fragments fragmentation = 1
demographic.data$frag[which(demographic.data$No.frag.t1 == 0 & demographic.data$surv == 1)] = 0 #<- if the colony survived but reported no fragmemts then fragmention = 0
table(demographic.data$frag)
# The data is now ready for analyses

###############################################
# STEP 2: Format fecundity data
###############################################
# Firstly I need to assign corals to their respective life-history strategies
table(fec.data$specie_name) # following the assigning of groups used previously Acropora corals are Competitive, Goniastrea is Stress Tolerant, Stylophora is weedy.
# There is no generalist corals in this dataset although these corals are being excluded from the analysis (see below).
# Add the additional variable
fec.data$LHS <- NA
fec.data[which(fec.data$specie_id %in% c(68,80,109,115)),]$LHS <- "C"
fec.data[which(fec.data$specie_id %in% c(805)),]$LHS <- "ST"
fec.data[which(fec.data$specie_id %in% c(1441)),]$LHS <- "W"
# check that has worked as expected
table(fec.data$specie_id); table(fec.data$LHS)

# this datafile needs to rearranging due to the odd layout of the original format.
# First extract colony fecundity data
indexfecundity = which(fec.data$trait_name == "Colony fecundity") 
fecundity <- as.vector(fec.data$value[indexfecundity]) 
fecundity <- as.numeric(fecundity) # these are our gonad densities but need adjusting to cm3
fecundity <- fecundity/100
# Now extract out the colony size data
indexsize <- which(fec.data$trait_name == "Colony area") 
size <- as.vector(fec.data$value[indexsize])
size <- as.numeric(size) #this is our size data - but it needs converting to the same format as the main demographic data file
# size doesn't need adjusting as it is already in cm2
size <- log(size) #this puts it into the same format as the data we are working with
# And finally extract LHS info
LHS <- as.vector(fec.data$LHS[indexsize])
LHS <- as.factor(LHS) 
# Bring the fecundity data back together
fecundity.data <- data.frame(Fecundity = fecundity, Size = size, LHS = LHS)

###############################################
# STEP 3: Estimate demographic parameters.
###############################################

# Firstly, what species are being dealt with in each Life-history category.
species.list <- demographic.data %>% 
                      group_by(LHS, Genus, Species) %>% dplyr::summarise()

# Given the different number of intervals for which we have data across the various sites I am going to be calculating vital rates with time pooled.
# Therefore all regression analyses will need to account for the random factor if site and colony ID.

# Also how many tagged colonies are being followed each year across the different groups (Report these as cumalative totals across all years in the manuscript)
overall.details <- summarySE(demographic.data, measurevar = "Size.t", groupvars = c("Country","Ecoregion","LHS","time.t"), na.rm = TRUE)
# given the limited quantity of data across the board for generalists species I am going to omit them from the study :(
demographic.data <- demographic.data[which(demographic.data$LHS %in% c("C", "ST", "W")),]

###### 1. Survival ----------------------------------------------------------------

# Visualise raw data
plot(surv ~ Size.t, col = LHS, data = demographic.data)

# Run regression models (how does including the random effect of site and colony ID influence the quality of the model?)
surv.mod <- glmer(surv ~ Size.t * LHS * Country * Ecoregion + (1|Colony.ID) + (Size.t|Site), family = "binomial", data = demographic.data)
surv.mod2 <- glmer(surv ~ Size.t * LHS * Country * Ecoregion + (1|Colony.ID) + (1|Site), family = "binomial", data = demographic.data)
surv.mod3 <- glmer(surv ~ Size.t * LHS * Country * Ecoregion + (1|Colony.ID), family = "binomial", data = demographic.data)
# which model fits best?
anova(surv.mod, surv.mod2, surv.mod3) # it is model 1 (According to AICs)
# also a quick check to see the importance of each term within the model - i.e. can it be simplified?
car::Anova(surv.mod) # nope - all variables are necessary at some point.

# Visualise vital rate patterns across the different groups
survival <- ggpredict(surv.mod, terms = c("Size.t[all]", "Ecoregion", "LHS", "Country"))

# Now plot the trends in survival and how it varies between groups and regions.
survival.plot <- ggplot(survival, aes(x = x)) +
  geom_line(aes(y = predicted, col = facet, linetype = group), size = 1.1) +
  facet_wrap(~panel) +
  # Formatting
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0,1)) +
  scale_color_manual(name = "", values = c("blue2","orange3","red2"), 
                     labels = c("Competitive    ","Stress-Tolerant","Weedy          "), #this is to give constant space for inserting pictures
                     guide = guide_legend(override.aes = list(shape = 15, size = 10))) +
  scale_linetype_manual(name = "", values = c("solid","dashed"), 
                        breaks = c("Tropical", "Subtropical"), 
                        labels = c("Tropical", "Subtropical"), 
                        guide = FALSE) +
  xlab("") + 
  ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0.9),axis.line.y.right = element_line(colour = "black", size = 1),
        legend.background = element_blank(), legend.key = element_blank(),legend.position = "bottom",
        legend.box = "horizontal")



###### 2. Growth ------------------------------------------------------------------

# A. Size transitions ----------------- 
# Check transitions in size
plot(Size.t1 ~ Size.t, col = LHS, data = demographic.data); abline(0,1, lwd = 2)

# Run regression models (how does including the random effect of site and colony ID influence the quality of the model?) whilst also comparing linear and non-linear
grow.mod <- lmer(Size.t1 ~ Size.t + I(Size.t^2) * LHS * Country * Ecoregion + (1|Colony.ID) + (Size.t|Site), data = demographic.data)
grow.modL <- lmer(Size.t1 ~ Size.t * LHS * Country * Ecoregion + (1|Colony.ID) + (Size.t|Site), data = demographic.data)
grow.mod2 <- lmer(Size.t1 ~ Size.t + I(Size.t^2) * LHS * Country * Ecoregion + (1|Colony.ID) + (1|Site), data = demographic.data)
grow.mod3 <- lmer(Size.t1 ~ Size.t + I(Size.t^2) * LHS * Country * Ecoregion + (1|Colony.ID), data = demographic.data)
# which model fits best?
anova(grow.mod, grow.modL, grow.mod2, grow.mod3) # it is model 1 (According to AICs) and this is non-linear
# also a quick check to see the importance of each term within the model - i.e. can it be simplified?
car::Anova(grow.mod) # nope - all variables are necessary at some point.

# Visualise vital rate patterns across the different groups
growth <- ggpredict(grow.mod, terms = c("Size.t[all]", "Ecoregion", "LHS", "Country"))

# Now plot the trends in survival and how it varies between groups and regions.
growth.plot <- ggplot(growth, aes(x = x)) +
  geom_line(aes(y = predicted, col = facet, linetype = group), size = 1.1) +
  geom_abline(linetype = "dotted", size = 1.1) +
  facet_wrap(~panel) +
  # Formatting
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(-2,10)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(-2,10)) +
  scale_color_manual(name = "", values = c("blue2","orange3","red2"), 
                     labels = c("Competitive    ","Stress-Tolerant","Weedy          "), #this is to give constant space for inserting pictures
                     guide = guide_legend(override.aes = list(shape = 15, size = 10))) +
  scale_linetype_manual(name = "", values = c("solid","dashed"), 
                        breaks = c("Tropical", "Subtropical"), 
                        labels = c("Tropical", "Subtropical"), 
                        guide = FALSE) +
  xlab("") + 
  ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0.9),axis.line.y.right = element_line(colour = "black", size = 1),
        legend.background = element_blank(), legend.key = element_blank(),legend.position = "bottom",
        legend.box = "horizontal")

# B. Growth variance ------------------
# How does the variance in size.t vary with initial colony size?

# Plot variance from the growth model
plot(grow.mod@frame$Size.t, abs(resid(grow.mod)), col = grow.mod@frame$LHS,
     xlab='Size at t', ylab='residual variance') 

# Evaluate the shape of the relationship (linear versus non-linear)
growth.sd.mod = glmer(abs(resid(grow.mod)) ~ grow.mod@frame$Size.t * grow.mod@frame$LHS * grow.mod@frame$Country * grow.mod@frame$Ecoregion + (1|grow.mod@frame$Colony.ID) + (grow.mod@frame$Size.t1|grow.mod@frame$Site) , family = Gamma(link = "log")) # the gamma fit used here allows for a non-linear relationship but prevents the variance from dropping below zero
growth.sd.modL = lmer(abs(resid(grow.mod)) ~ grow.mod@frame$Size.t * grow.mod@frame$LHS * grow.mod@frame$Country * grow.mod@frame$Ecoregion + (1|grow.mod@frame$Colony.ID) + (grow.mod@frame$Size.t1|grow.mod@frame$Site))
# which models fit best
AIC(growth.sd.mod, growth.sd.modL) 
BIC(growth.sd.mod, growth.sd.modL) # non-linear relationship fits the best.
# a variance in growth plot isn't needed for the manuscript.



###### 3. Fragmentation -----------------------------------------------------------

# For some of the components that make up this vital rate I am only interested in the colonies that fragmented
# and those produced by fragmentation events.
# Fragmenting colonies
colonies.frag.TRUE <- demographic.data[which(demographic.data$frag == 1),] 
# Fragments produced
frag.corals <- demographic.data[which(demographic.data$Fragment.t1 == "Yes"),]
# and fill in the blank size.t space so that I can link the parent colony size with each fragment.
frag.corals <- fill(frag.corals, Size.t, .direction = "down")

# A. Fragmentation probability --------
# For this I am still interested in all colonies (fragmented or not)

# Visualise raw data
plot(frag ~ Size.t, col = Ecoregion, data = demographic.data) # appears to be a slight intermediate peak in fragmentation

# Run regression models (how does including the random effect of site and colony ID influence the quality of the model?)
frag.mod <- glmer(frag ~ Size.t * LHS * Country * Ecoregion + (1|Colony.ID) + (Size.t|Site), family = "binomial", data = demographic.data)
frag.mod.NL <- glmer(frag ~ Size.t + I(Size.t^2) * LHS * Country * Ecoregion + (1|Colony.ID) + (Size.t|Site), family = "binomial", data = demographic.data) # this polynomial format might help to capture the intermediate peak in fragentation
frag.mod2 <- glmer(frag ~ Size.t * LHS * Country * Ecoregion + (1|Colony.ID) + (1|Site), family = "binomial", data = demographic.data)
frag.mod3 <- glmer(frag ~ Size.t * LHS * Country * Ecoregion + (1|Colony.ID), family = "binomial", data = demographic.data)
# which model fits best?
anova(frag.mod, frag.mod.NL, frag.mod2, frag.mod3) # it is model 1 (According to AICs).
# also a quick check to see the importance of each term within the model - i.e. can it be simplified?
car::Anova(frag.mod) # nope - all variables are necessary at some point.

# Visualise vital rate patterns across the different groups
fragmentation <- ggpredict(frag.mod, terms = c("Size.t[all]", "Ecoregion", "LHS", "Country"))
fragmentation2 <- ggpredict(frag.mod.NL, terms = c("Size.t[all]", "Ecoregion", "LHS", "Country"))
anova(frag.mod, frag.mod.NL) #Although the AIC scores suggests that model 1 is the best, the non-linear model is not significantly different and visually is a better fit for the data.
# rename the country variable in the predicted dataframe for combining plots
colnames(fragmentation2) <- colnames(fragmentation) <- c("x", "predicted", "Ecoregion", "LHS", "Country")   

# Now plot the trends in survival and how it varies between groups and regions.
fragmentation.plot <- ggplot(fragmentation2, aes(x = x)) +
  geom_line(aes(y = predicted, col = LHS, linetype = Ecoregion), size = 1.1) +
  geom_point(aes(x = Size.t, y = frag, col = LHS, shape = Ecoregion),
             size = 1.1, position = position_jitter(0,0.1), data = demographic.data) +
  facet_grid(LHS~Country) +
  # Formatting
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)), limits = c(0,1)) +
  scale_color_manual(name = "", values = c("blue2","orange3","red2"), 
                     labels = c("Competitive    ","Stress-Tolerant","Weedy          "), #this is to give constant space for inserting pictures
                     guide = guide_legend(override.aes = list(shape = 15, size = 10))) +
  scale_shape_manual(breaks = c("Tropical", "Subtropical"), # open = Subtropical, closed = Tropical
                     values = c(16, 0),
                     guide = FALSE)+
  scale_linetype_manual(name = "", values = c("solid","dashed"), 
                        breaks = c("Tropical", "Subtropical"), 
                        labels = c("Tropical", "Subtropical"), 
                        guide = FALSE) +
  xlab("") + 
  ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0.9),axis.line.y.right = element_line(colour = "black", size = 1),
        legend.background = element_blank(), legend.key = element_blank(),legend.position = "bottom",
        legend.box = "horizontal")
# Non-linear fragmentation fits the data better. 

# B. Number of fragments --------------

# Visualise raw data
plot(No.frag.t1 ~ Size.t, col = LHS, data = colonies.frag.TRUE) # production definitely increases with colony size
# However I need to count how many fragmentation events there were in each group to check I have enough to model the groups separately across countries and ecoregions
# overall we observed 191 cases of fragmentation.
frag.totals <- summarySE(colonies.frag.TRUE, measurevar = "Size.t1", groupvars = c("Country","Ecoregion","LHS"), na.rm = TRUE)[,1:4] # yeah there are some sites where there isn't enough data for modeling this.
# to get around this I will group corals by country
summarySE(colonies.frag.TRUE, measurevar = "Size.t1", groupvars = c("Country","LHS"), na.rm = TRUE)[,1:4] # that is better.
# so we are making the assumption that the number of fragments produced as a function of size is not different between tropical and subtropical corals

# Linear model
no.frag.mod <- glmer(No.frag.t1 ~ Size.t * LHS * Country + (1|Colony.ID) + (Size.t|Site), family = "poisson", data = colonies.frag.TRUE) # a poisson distribution is required due to the fact that this is colony counts
# Non-linear model for comparison
no.frag.mod.NL <- glmer(No.frag.t1 ~ Size.t + I(Size.t^2) * LHS * Country + (1|Colony.ID) + (Size.t|Site), family = "poisson", data = colonies.frag.TRUE)
# which model fits the best?
anova(no.frag.mod, no.frag.mod.NL) # AIC scores suggest the linear model is best

# Now view the trend in the number of fragments
No.fragments <- ggpredict(no.frag.mod, terms = c("Size.t[all]", "Country", "LHS"))

# Plot the number of frags produced against the size of colonies.
No.fragments.plot <- ggplot(No.fragments, aes(x = x)) +
  geom_line(aes(y = predicted, col = facet), size = 1.1) +
  facet_wrap(~group) +
  # Formatting
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_color_manual(name = "", values = c("blue2","orange2","red2"), 
                     labels = c("Competitive    ","Stress-Tolerant","Weedy          "), #this is to give constant space for inserting pictures
                     guide = guide_legend(override.aes = list(shape = 15, size = 10))) +
  xlab("") + 
  ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0.9),axis.line.y.right = element_line(colour = "black", size = 1),
        legend.background = element_blank(), legend.key = element_blank(),legend.position = "bottom",
        legend.box = "horizontal")

# C. Size of fragments ----------------

# Visualise raw data
plot(Size.t1 ~ Size.t, col = LHS, data = frag.corals) # production definitely increases with colony size
# again just check how much data there is - because if there is not a lot of No.frag~size data then there is likely to be not a lot of frag.size~size data.
summarySE(frag.corals, measurevar = "Size.t1", groupvars = c("Country", "Ecoregion", "LHS"), na.rm = TRUE)
summarySE(frag.corals, measurevar = "Size.t1", groupvars = c("Country", "LHS"), na.rm = TRUE) # this is better.

# Run regression models (there is not enough repetition across individuals or sites to include random factors here)
# linear model
frag.size.mod <- glm(Size.t1 ~ Size.t * LHS * Country, data = frag.corals) # here Size.t1 is fragment size and Size t is parent colony size
# non-linear model
frag.size.mod.NL <- glm(Size.t1 ~ Size.t + I(Size.t^2) * LHS * Country, data = frag.corals)
# Which model fits best
AIC(frag.size.mod,frag.size.mod.NL)
BIC(frag.size.mod,frag.size.mod.NL) # Non-linear is the best according to the AICs

# Extract the relevant variables.
size.frag <- ggpredict(frag.size.mod.NL, terms = c("Size.t[all]", "Country", "LHS"))

# Plot mean frag size against parent colony size. Include reference line.
size.of.frags <- ggplot(size.frag, aes(x = x)) +
  geom_line(aes(y = predicted, col = facet), size = 1.1) +
  geom_abline(linetype = "dotted", size = 1.1) +
  facet_wrap(~group) +
  # Formatting
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_color_manual(name = "", values = c("blue2","orange2","red2"), 
                     labels = c("Competitive    ","Stress-Tolerant","Weedy          "), #this is to give constant space for inserting pictures
                     guide = guide_legend(override.aes = list(shape = 15, size = 10))) +
  xlab("") + 
  ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0.9),axis.line.y.right = element_line(colour = "black", size = 1),
        legend.background = element_blank(), legend.key = element_blank(),legend.position = "bottom",
        legend.box = "horizontal")
# There is a case where small colonies produce fragments larger than themselves but this will be restricted by the limited likelihood of small corals from fragmenting

# How does the variance in fragment sizes vary with parent colony size.
plot(frag.size.mod.NL$model$Size.t, abs(resid(frag.size.mod.NL)), col = frag.size.mod.NL$model$LHS,
     xlab='Size at t',ylab='residual variance') 

# what is the shape of this relationship
frag.size.sd.mod.NL = glm(abs(resid(frag.size.mod.NL)) ~ frag.size.mod.NL$model$Size.t * frag.size.mod.NL$model$LHS * frag.size.mod.NL$model$Country, family = Gamma(link = "log")) # the gamma fit used here allows for a non-linear relationship but prevents the variance from dropping below zero
frag.size.sd.mod = lm(abs(resid(frag.size.mod.NL)) ~ frag.size.mod.NL$model$Size.t * frag.size.mod.NL$model$LHS * frag.size.mod.NL$model$Country)
# which models fit best
AIC(frag.size.sd.mod.NL, frag.size.sd.mod) 
BIC(frag.size.sd.mod.NL, frag.size.sd.mod) # non-linear relationship fits the best.
# a variance in fragment size plot isn't needed for the manuscript.



###### 4. Recruitment -------------------------------------------------------------
# Here I am only interested in recorded recruits.
recruit.data <- demographic.data[which(demographic.data$Recruit.t1 == "Yes"),]

# A. Number of recruits ---------------

# How many recruits from each life-history category were observed across the years in each location.
recruits.counts <- summarySE(recruit.data, measurevar = "Size.t1", groupvars = c("Country","Ecoregion","LHS","time.t"), na.rm = TRUE)[,1:5]

# B. Recruit size ---------------------

# Check the overall recruit size distribution
hist(recruit.data$Size.t1)

# Run a regression model to extra the required size parameters.
# Since the models being constructed are open models, recruit size is to be modeled independently of parent colony details
recruit.size.mod <- lm(Size.t1 ~ LHS * Country * Ecoregion, data = recruit.data) 
# just to the lack of recruits observed across some coral groups in particular years I will model recruit size independent of year.
car::Anova(recruit.size.mod)

# Calculate the size distribution of juvenile colony sizes for plotting.
# Competitive
C_S_A <- density(recruit.data[which(recruit.data$LHS == "C" & recruit.data$Country == "Australia" & recruit.data$Ecoregion == "Subtropical"),]$Size.t1)
C_T_A <- density(recruit.data[which(recruit.data$LHS == "C" & recruit.data$Country == "Australia" & recruit.data$Ecoregion == "Tropical"),]$Size.t1)
C_S_J <- density(recruit.data[which(recruit.data$LHS == "C" & recruit.data$Country == "Japan" & recruit.data$Ecoregion == "Subtropical"),]$Size.t1)
C_T_J <- density(recruit.data[which(recruit.data$LHS == "C" & recruit.data$Country == "Japan" & recruit.data$Ecoregion == "Tropical"),]$Size.t1)
# Stress-Tolerant
ST_S_A <- density(recruit.data[which(recruit.data$LHS == "ST" & recruit.data$Country == "Australia" & recruit.data$Ecoregion == "Subtropical"),]$Size.t1)
ST_T_A <- density(recruit.data[which(recruit.data$LHS == "ST" & recruit.data$Country == "Australia" & recruit.data$Ecoregion == "Tropical"),]$Size.t1)
ST_S_J <- density(recruit.data[which(recruit.data$LHS == "ST" & recruit.data$Country == "Japan" & recruit.data$Ecoregion == "Subtropical"),]$Size.t1)
ST_T_J <- density(recruit.data[which(recruit.data$LHS == "ST" & recruit.data$Country == "Japan" & recruit.data$Ecoregion == "Tropical"),]$Size.t1)
# Weedy
W_S_A <- density(recruit.data[which(recruit.data$LHS == "W" & recruit.data$Country == "Australia" & recruit.data$Ecoregion == "Subtropical"),]$Size.t1)
W_T_A <- density(recruit.data[which(recruit.data$LHS == "W" & recruit.data$Country == "Australia" & recruit.data$Ecoregion == "Tropical"),]$Size.t1)
W_S_J <- density(recruit.data[which(recruit.data$LHS == "W" & recruit.data$Country == "Japan" & recruit.data$Ecoregion == "Subtropical"),]$Size.t1)
W_T_J <- density(recruit.data[which(recruit.data$LHS == "W" & recruit.data$Country == "Japan" & recruit.data$Ecoregion == "Tropical"),]$Size.t1)
# and condense the data together.
recruit.sizes <- data.frame(Country = rep(c("Australia","Japan","Australia","Japan","Australia","Japan"), each = 512),
                            LHS = rep(c("C","ST","W"), each = 1024),
                            S.X = c(C_S_A$x,C_S_J$x,ST_S_A$x,ST_S_J$x,W_S_A$x,W_S_J$x),
                            S.Y = c(C_S_A$y,C_S_J$y,ST_S_A$y,ST_S_J$y,W_S_A$y,W_S_J$y),
                            T.X = c(C_T_A$x,C_T_J$x,ST_T_A$x,ST_T_J$x,W_T_A$x,W_T_J$x),
                            T.Y = c(C_T_A$y,C_T_J$y,ST_T_A$y,ST_T_J$y,W_T_A$y,W_T_J$y))
                            
# plot the densities together to visualise the similarity across regions and years.
recruit.size.frequency <- ggplot(recruit.sizes) +
  geom_area(aes(x = S.X, y = S.Y, fill = LHS, colour = LHS), linetype = "dashed", 
            alpha = 0.6, size = 1.1) +
  geom_area(aes(x = T.X, y = T.Y, fill = LHS, colour = LHS), linetype = "solid", 
            alpha = 0.8, size = 1.1) +
  facet_grid(LHS ~ Country) +
  # Formatting
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(-5,7)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0,0.5)) +
  scale_fill_manual(name = "", values = c("blue2","orange3","red2"), 
                    labels = c("Competitive    ","Stress-Tolerant","Weedy          "), #this is to give constant space for inserting pictures
                    guide = guide_legend(override.aes = list(shape = 15, size = 10))) +
  scale_color_manual(name = "", values = c("blue2","orange3","red2"), 
                    labels = c("Competitive    ","Stress-Tolerant","Weedy          "), #this is to give constant space for inserting pictures
                    guide = FALSE) +
  xlab("") + 
  ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0.9),axis.line.y.right = element_line(colour = "black", size = 1),
        legend.background = element_blank(), legend.key = element_blank(),legend.position = "bottom",
        legend.box = "horizontal", panel.spacing.y = unit(1.5, "lines"))

# The recruit size distributions highlights remarkable similarities across the countries and ecoregions for each coral group.
# Therefore I will model recruit size distributions independent of country and ecoregion - this also helps to deal with the lack of recruits observed for the generalists in subtropical Australia.
recruit.size.mod2 <- lm(Size.t1 ~ LHS, data = recruit.data)



###### 5. Fecundity ---------------------------------------------------------------
# This vital rate is purely so that I can combine the dynamics of existing colonies with recruitment patterns.
# I just need to relationship between colony size and its fecundity. During parameterisation the IPM kernel will constrain fecundity within observed limits.

# Visualise the data
plot(Fecundity ~ Size, col = LHS, data = fecundity.data) # very exponential relationship across groups

# estimate the relationship
fec.mod <- glm(Fecundity ~ Size + I(Size^2) * LHS, data = fecundity.data)
fec.mod2 <- nls(Fecundity ~ a * exp(b * Size), data = fecundity.data[which(fecundity.data$LHS == "C"),], start = list(a = 1, b = 1))
fec.mod3 <- nls(Fecundity ~ a * exp(b * Size), data = fecundity.data[which(fecundity.data$LHS == "ST"),], start = list(a = 1, b = 1))
fec.mod4 <- nls(Fecundity ~ a * exp(b * Size), data = fecundity.data[which(fecundity.data$LHS == "W"),], start = list(a = 1, b = 1))

# Check the model predictions
fecundity <- ggpredict(fec.mod, terms = c("Size[all]", "LHS"))
fecundity2 <- data.frame(x = fecundity.data$Size, 
                         y = predict(fec.mod2, list(Size = fecundity.data$Size))) # Competitive
fecundity3 <- data.frame(x = fecundity.data$Size, 
                         y = predict(fec.mod3, list(Size = fecundity.data$Size))) # Stress tolerant
fecundity4 <- data.frame(x = fecundity.data$Size, 
                         y = predict(fec.mod4, list(Size = fecundity.data$Size))) # Weedy

# visually check 
ggplot(fecundity.data, aes(x = Size)) +
  geom_point(aes(y = Fecundity, col = LHS)) +
  geom_line(aes(x = x, y = predicted, col = group), size = 1.1, data = fecundity) +
  geom_line(aes(x = x, y = y), col = "red", linetype = "dashed", data = fecundity2) +
  geom_line(aes(x = x, y = y), col = "green", linetype = "dashed", data = fecundity3) +
  geom_line(aes(x = x, y = y), col = "blue", linetype = "dashed", data = fecundity4)
# the nls regression models offer the best prediction.

# Plot proportional colony fecundity against the size of colonies.
fecundity.plot <- ggplot(fecundity2, aes(x = Size)) +
  geom_line(aes(x = x, y = y, col = "blue2"), size = 1.1) +  
  geom_line(aes(x = x, y = y, col = "orange2"), size = 1.1, data = fecundity3) +
  geom_line(aes(x = x, y = y, col = "red2"), size = 1.1, data = fecundity4) +
  # Formatting
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_color_manual(name = "", values = c("blue2","orange2","red2"), 
                     labels = c("Competitive    ","Stress-Tolerant","Weedy          "), #this is to give constant space for inserting pictures
                     guide = guide_legend(override.aes = list(shape = 15, size = 10))) +
  xlab("") + 
  ylab("") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        panel.border = element_blank(), axis.ticks = element_line(size = 0.9),axis.line.y.right = element_line(colour = "black", size = 1),
        legend.background = element_blank(), legend.key = element_blank(),legend.position = "bottom",
        legend.box = "horizontal")


###############################################
# STEP 4: Store demographic parameters.
###############################################
# Store the relevant coefficients for use in constructing the sub kernels necessary for comparing the dynamics of the different coral groups.
# The IPMs I will be constructing will consist of G (Survival and Growth), and H (Fragmentation) and F (recruitment) subkernels which will be condensed in to a single K kernel.

# Because of the way the IPM functions are defined in the 'IPM construction' script I need to extract coefficients for each group separately.
# firstly because of the inclusion of random effects of individuals in some of the vital rate regression models I need to extract the relevant mean coefficients.
survival.coefs <- apply(coef(surv.mod)$Colony.ID, 2, mean)
growth.coefs <- apply(coef(grow.mod)$Colony.ID, 2, mean)
growth.SD.coefs <- apply(coef(growth.sd.mod)$'grow.mod@frame$Colony.ID', 2, mean); growth.SD.coefs <- growth.SD.coefs[2:25]
frag.coefs <- apply(coef(frag.mod.NL)$Colony.ID, 2, mean)
frag.no.coefs <- apply(coef(no.frag.mod)$Colony.ID, 2, mean)

# Extract coefficients for each coral population

########################################################## Competitive ########

# Competitive (Australia, Subtropical)
coefs.C.A.S <- c(
  # Survival probability
  surv.int               =  survival.coefs[1],
  surv.slope             =  survival.coefs[2],
  # Growth transitions
  grow.int               =  growth.coefs[1],
  grow.slope             =  growth.coefs[2],
  grow.slope2            =  growth.coefs[3],
  # Growth variance
  grow.sd.int            =  growth.SD.coefs[1],
  grow.sd.slope          =  growth.SD.coefs[2],
  # Fragmentation probability
  frag.int               =  frag.coefs[1],
  frag.slope             =  frag.coefs[2],
  frag.slope2            =  frag.coefs[3],
  # Number of frags produced  
  no.frag.int            =  frag.no.coefs[1],
  no.frag.slope          =  frag.no.coefs[2],
  # Fragment size
  frag.s.int             =  coef(frag.size.mod.NL)[1],
  frag.s.slope           =  coef(frag.size.mod.NL)[2],
  frag.s.slope2          =  coef(frag.size.mod.NL)[3],
  # Variance in fragment size
  frag.sd.int            =  coef(frag.size.sd.mod.NL)[1],
  frag.sd.slope          =  coef(frag.size.sd.mod.NL)[2],
  # Colony fecundity 
  fec.a                  =  coef(fec.mod2)[1],
  fec.b                  =  coef(fec.mod2)[2],
  # Number of recruits (this will allow the model to vary between observed limits)
  rec.max                =  max(recruits.counts[which(recruits.counts$Country == "Australia" &
                                                        recruits.counts$Ecoregion == "Subtropical" &
                                                        recruits.counts$LHS == "C"), "N"]),
  rec.min                =  min(recruits.counts[which(recruits.counts$Country == "Australia" &
                                                        recruits.counts$Ecoregion == "Subtropical" &
                                                        recruits.counts$LHS == "C"), "N"]),
  # recruit size coefficients
  rec.size.int           =  coefficients(recruit.size.mod2)[1],
  rec.size.sd            =  sigma(recruit.size.mod2))

# Competitive (Australia, Tropical)
coefs.C.A.T <- c(
  # Survival probability
  surv.int               =  survival.coefs[1] + survival.coefs[6],
  surv.slope             =  survival.coefs[2] + survival.coefs[12],
  # Growth transitions
  grow.int               =  growth.coefs[1] + growth.coefs[7],
  grow.slope             =  growth.coefs[2],
  grow.slope2            =  growth.coefs[3] + growth.coefs[13],
  # Growth variance
  grow.sd.int            =  growth.SD.coefs[1] + growth.SD.coefs[6],
  grow.sd.slope          =  growth.SD.coefs[2] + growth.SD.coefs[12],
  # Fragmentation probability
  frag.int               =  frag.coefs[1] + frag.coefs[7],
  frag.slope             =  frag.coefs[2],
  frag.slope2            =  frag.coefs[3] + frag.coefs[13],
  # Number of frags produced  
  no.frag.int            =  frag.no.coefs[1],
  no.frag.slope          =  frag.no.coefs[2],
  # Fragment size
  frag.s.int             =  coef(frag.size.mod.NL)[1],
  frag.s.slope           =  coef(frag.size.mod.NL)[2],
  frag.s.slope2          =  coef(frag.size.mod.NL)[3],
  # Variance in fragment size
  frag.sd.int            =  coef(frag.size.sd.mod.NL)[1],
  frag.sd.slope          =  coef(frag.size.sd.mod.NL)[2],
  # Colony fecundity 
  fec.a                  =  coef(fec.mod2)[1],
  fec.b                  =  coef(fec.mod2)[2],
  # Number of recruits (this will allow the model to vary between observed limits)
  rec.max                =  max(recruits.counts[which(recruits.counts$Country == "Australia" &
                                                        recruits.counts$Ecoregion == "Tropical" &
                                                        recruits.counts$LHS == "C"), "N"]),
  rec.min                =  min(recruits.counts[which(recruits.counts$Country == "Australia" &
                                                        recruits.counts$Ecoregion == "Tropical" &
                                                        recruits.counts$LHS == "C"), "N"]),
  # recruit size coefficients
  rec.size.int           =  coefficients(recruit.size.mod2)[1],
  rec.size.sd            =  sigma(recruit.size.mod2))

# Competitive (Japan, Subtropical)
coefs.C.J.S <- c(
  # Survival probability
  surv.int               =  survival.coefs[1] + survival.coefs[5],
  surv.slope             =  survival.coefs[2] + survival.coefs[9],
  # Growth transitions
  grow.int               =  growth.coefs[1] + growth.coefs[6],
  grow.slope             =  growth.coefs[2],
  grow.slope2            =  growth.coefs[3] + growth.coefs[10],
  # Growth variance
  grow.sd.int            =  growth.SD.coefs[1] + growth.SD.coefs[5],
  grow.sd.slope          =  growth.SD.coefs[2] + growth.SD.coefs[9],
  # Fragmentation probability
  frag.int               =  frag.coefs[1] + frag.coefs[6],
  frag.slope             =  frag.coefs[2],
  frag.slope2            =  frag.coefs[3] + frag.coefs[10],
  # Number of frags produced  
  no.frag.int            =  frag.no.coefs[1] + frag.no.coefs[5],
  no.frag.slope          =  frag.no.coefs[2] + frag.no.coefs[8],
  # Fragment size
  frag.s.int             =  coef(frag.size.mod.NL)[1] + coef(frag.size.mod.NL)[6],
  frag.s.slope           =  coef(frag.size.mod.NL)[2],
  frag.s.slope2          =  coef(frag.size.mod.NL)[3] + coef(frag.size.mod.NL)[9],
  # Variance in fragment size
  frag.sd.int            =  coef(frag.size.sd.mod.NL)[1] + coef(frag.size.sd.mod.NL)[5],
  frag.sd.slope          =  coef(frag.size.sd.mod.NL)[2] + coef(frag.size.sd.mod.NL)[8],
  # Colony fecundity 
  fec.a                  =  coef(fec.mod2)[1],
  fec.b                  =  coef(fec.mod2)[2],
  # Number of recruits (this will allow the model to vary between observed limits)
  rec.max                =  max(recruits.counts[which(recruits.counts$Country == "Japan" &
                                                        recruits.counts$Ecoregion == "Subtropical" &
                                                        recruits.counts$LHS == "C"), "N"]),
  rec.min                =  min(recruits.counts[which(recruits.counts$Country == "Japan" &
                                                        recruits.counts$Ecoregion == "Subtropical" &
                                                        recruits.counts$LHS == "C"), "N"]),
  # recruit size coefficients
  rec.size.int           =  coefficients(recruit.size.mod2)[1],
  rec.size.sd            =  sigma(recruit.size.mod2))

# Competitive (Japan, Tropical)
coefs.C.J.T <- c(
  # Survival probability
  surv.int               =  survival.coefs[1] + survival.coefs[15],
  surv.slope             =  survival.coefs[2] + survival.coefs[20],
  # Growth transitions
  grow.int               =  growth.coefs[1] + growth.coefs[16],
  grow.slope             =  growth.coefs[2],
  grow.slope2            =  growth.coefs[3] + growth.coefs[21],
  # Growth variance
  grow.sd.int            =  growth.SD.coefs[1] + growth.SD.coefs[15],
  grow.sd.slope          =  growth.SD.coefs[2] + growth.SD.coefs[20],
  # Fragmentation probability
  frag.int               =  frag.coefs[1] + frag.coefs[16],
  frag.slope             =  frag.coefs[2],
  frag.slope2            =  frag.coefs[3] + frag.coefs[21],
  # Number of frags produced  
  no.frag.int            =  frag.no.coefs[1] + frag.no.coefs[5],
  no.frag.slope          =  frag.no.coefs[2] + frag.no.coefs[8],
  # Fragment size
  frag.s.int             =  coef(frag.size.mod.NL)[1] + coef(frag.size.mod.NL)[6],
  frag.s.slope           =  coef(frag.size.mod.NL)[2],
  frag.s.slope2          =  coef(frag.size.mod.NL)[3] + coef(frag.size.mod.NL)[9],
  # Variance in fragment size
  frag.sd.int            =  coef(frag.size.sd.mod.NL)[1] + coef(frag.size.sd.mod.NL)[5],
  frag.sd.slope          =  coef(frag.size.sd.mod.NL)[2] + coef(frag.size.sd.mod.NL)[8],
  # Colony fecundity 
  fec.a                  =  coef(fec.mod2)[1],
  fec.b                  =  coef(fec.mod2)[2],
  # Number of recruits (this will allow the model to vary between observed limits)
  rec.max                =  max(recruits.counts[which(recruits.counts$Country == "Japan" &
                                                        recruits.counts$Ecoregion == "Tropical" &
                                                        recruits.counts$LHS == "C"), "N"]),
  rec.min                =  min(recruits.counts[which(recruits.counts$Country == "Japan" &
                                                        recruits.counts$Ecoregion == "Tropical" &
                                                        recruits.counts$LHS == "C"), "N"]),
  # recruit size coefficients
  rec.size.int           =  coefficients(recruit.size.mod2)[1],
  rec.size.sd            =  sigma(recruit.size.mod2))

########################################################## Stress-Tolerant ########

# Stress-Tolerant (Australia, Subtropical)
coefs.ST.A.S <- c(
  # Survival probability
  surv.int               =  survival.coefs[1] + survival.coefs[3],
  surv.slope             =  survival.coefs[2] + survival.coefs[7],
  # Growth transitions
  grow.int               =  growth.coefs[1] + growth.coefs[4],
  grow.slope             =  growth.coefs[2],
  grow.slope2            =  growth.coefs[3] + growth.coefs[8],
  # Growth variance
  grow.sd.int            =  growth.SD.coefs[1] + growth.SD.coefs[3],
  grow.sd.slope          =  growth.SD.coefs[2] + growth.SD.coefs[7],
  # Fragmentation probability
  frag.int               =  frag.coefs[1] + frag.coefs[4],
  frag.slope             =  frag.coefs[2],
  frag.slope2            =  frag.coefs[3] + frag.coefs[8],
  # Number of frags produced  
  no.frag.int            =  frag.no.coefs[1] + frag.no.coefs[3],
  no.frag.slope          =  frag.no.coefs[2] + frag.no.coefs[6],
  # Fragment size
  frag.s.int             =  coef(frag.size.mod.NL)[1] + coef(frag.size.mod.NL)[4],
  frag.s.slope           =  coef(frag.size.mod.NL)[2],
  frag.s.slope2          =  coef(frag.size.mod.NL)[3] + coef(frag.size.mod.NL)[7],
  # Variance in fragment size
  frag.sd.int            =  coef(frag.size.sd.mod.NL)[1] + coef(frag.size.sd.mod.NL)[3],
  frag.sd.slope          =  coef(frag.size.sd.mod.NL)[2] + coef(frag.size.sd.mod.NL)[6],
  # Colony fecundity 
  fec.a                  =  coef(fec.mod3)[1],
  fec.b                  =  coef(fec.mod3)[2],
  # Number of recruits (this will allow the model to vary between observed limits)
  rec.max                =  max(recruits.counts[which(recruits.counts$Country == "Australia" &
                                                        recruits.counts$Ecoregion == "Subtropical" &
                                                        recruits.counts$LHS == "ST"), "N"]),
  rec.min                =  min(recruits.counts[which(recruits.counts$Country == "Australia" &
                                                        recruits.counts$Ecoregion == "Subtropical" &
                                                        recruits.counts$LHS == "ST"), "N"]),
  # recruit size coefficients
  rec.size.int           =  coefficients(recruit.size.mod2)[1] + coefficients(recruit.size.mod2)[2],
  rec.size.sd            =  sigma(recruit.size.mod2))

# Stress-Tolerant (Australia, Tropical)
coefs.ST.A.T <- c(
  # Survival probability
  surv.int               =  survival.coefs[1] + survival.coefs[13],
  surv.slope             =  survival.coefs[2] + survival.coefs[18],
  # Growth transitions
  grow.int               =  growth.coefs[1] + growth.coefs[14],
  grow.slope             =  growth.coefs[2],
  grow.slope2            =  growth.coefs[3] + growth.coefs[19],
  # Growth variance
  grow.sd.int            =  growth.SD.coefs[1] + growth.SD.coefs[13],
  grow.sd.slope          =  growth.SD.coefs[2] + growth.SD.coefs[18],
  # Fragmentation probability
  frag.int               =  frag.coefs[1] + frag.coefs[14],
  frag.slope             =  frag.coefs[2],
  frag.slope2            =  frag.coefs[3] + frag.coefs[19],
  # Number of frags produced  
  no.frag.int            =  frag.no.coefs[1] + frag.no.coefs[3],
  no.frag.slope          =  frag.no.coefs[2] + frag.no.coefs[6],
  # Fragment size
  frag.s.int             =  coef(frag.size.mod.NL)[1] + coef(frag.size.mod.NL)[4],
  frag.s.slope           =  coef(frag.size.mod.NL)[2],
  frag.s.slope2          =  coef(frag.size.mod.NL)[3] + coef(frag.size.mod.NL)[7],
  # Variance in fragment size
  frag.sd.int            =  coef(frag.size.sd.mod.NL)[1] + coef(frag.size.sd.mod.NL)[3],
  frag.sd.slope          =  coef(frag.size.sd.mod.NL)[2] + coef(frag.size.sd.mod.NL)[6],
  # Colony fecundity 
  fec.a                  =  coef(fec.mod3)[1],
  fec.b                  =  coef(fec.mod3)[2],
  # Number of recruits (this will allow the model to vary between observed limits)
  rec.max                =  max(recruits.counts[which(recruits.counts$Country == "Australia" &
                                                        recruits.counts$Ecoregion == "Tropical" &
                                                        recruits.counts$LHS == "ST"), "N"]),
  rec.min                =  min(recruits.counts[which(recruits.counts$Country == "Australia" &
                                                        recruits.counts$Ecoregion == "Tropical" &
                                                        recruits.counts$LHS == "ST"), "N"]),
  # recruit size coefficients
  rec.size.int           =  coefficients(recruit.size.mod2)[1] + coefficients(recruit.size.mod2)[2],
  rec.size.sd            =  sigma(recruit.size.mod2))

# Stress-Tolerant (Japan, Subtropical)
coefs.ST.J.S <- c(
  # Survival probability
  surv.int               =  survival.coefs[1] + survival.coefs[10],
  surv.slope             =  survival.coefs[2] + survival.coefs[16],
  # Growth transitions
  grow.int               =  growth.coefs[1] + growth.coefs[11],
  grow.slope             =  growth.coefs[2],
  grow.slope2            =  growth.coefs[3] + growth.coefs[17],
  # Growth variance
  grow.sd.int            =  growth.SD.coefs[1] + growth.SD.coefs[10],
  grow.sd.slope          =  growth.SD.coefs[2] + growth.SD.coefs[16],
  # Fragmentation probability
  frag.int               =  frag.coefs[1] + frag.coefs[11],
  frag.slope             =  frag.coefs[2],
  frag.slope2            =  frag.coefs[3] + frag.coefs[17],
  # Number of frags produced  
  no.frag.int            =  frag.no.coefs[1] + frag.no.coefs[9],
  no.frag.slope          =  frag.no.coefs[2] + frag.no.coefs[11],
  # Fragment size
  frag.s.int             =  coef(frag.size.mod.NL)[1] + coef(frag.size.mod.NL)[10],
  frag.s.slope           =  coef(frag.size.mod.NL)[2],
  frag.s.slope2          =  coef(frag.size.mod.NL)[3] + coef(frag.size.mod.NL)[12],
  # Variance in fragment size
  frag.sd.int            =  coef(frag.size.sd.mod.NL)[1] + coef(frag.size.sd.mod.NL)[9],
  frag.sd.slope          =  coef(frag.size.sd.mod.NL)[2] + coef(frag.size.sd.mod.NL)[11],
  # Colony fecundity 
  fec.a                  =  coef(fec.mod3)[1],
  fec.b                  =  coef(fec.mod3)[2],
  # Number of recruits (this will allow the model to vary between observed limits)
  rec.max                =  max(recruits.counts[which(recruits.counts$Country == "Japan" &
                                                        recruits.counts$Ecoregion == "Subtropical" &
                                                        recruits.counts$LHS == "ST"), "N"]),
  rec.min                =  min(recruits.counts[which(recruits.counts$Country == "Japan" &
                                                        recruits.counts$Ecoregion == "Subtropical" &
                                                        recruits.counts$LHS == "ST"), "N"]),
  # recruit size coefficients
  rec.size.int           =  coefficients(recruit.size.mod2)[1] + coefficients(recruit.size.mod2)[2],
  rec.size.sd            =  sigma(recruit.size.mod2))

# Stress-Tolerant (Japan, Tropical)
coefs.ST.J.T <- c(
  # Survival probability
  surv.int               =  survival.coefs[1] + survival.coefs[21],
  surv.slope             =  survival.coefs[2] + survival.coefs[23],
  # Growth transitions
  grow.int               =  growth.coefs[1] + growth.coefs[22],
  grow.slope             =  growth.coefs[2],
  grow.slope2            =  growth.coefs[3] + growth.coefs[24],
  # Growth variance
  grow.sd.int            =  growth.SD.coefs[1] + growth.SD.coefs[21],
  grow.sd.slope          =  growth.SD.coefs[2] + growth.SD.coefs[23],
  # Fragmentation probability
  frag.int               =  frag.coefs[1] + frag.coefs[22],
  frag.slope             =  frag.coefs[2],
  frag.slope2            =  frag.coefs[3] + frag.coefs[24],
  # Number of frags produced  
  no.frag.int            =  frag.no.coefs[1] + frag.no.coefs[9],
  no.frag.slope          =  frag.no.coefs[2] + frag.no.coefs[11],
  # Fragment size
  frag.s.int             =  coef(frag.size.mod.NL)[1] + coef(frag.size.mod.NL)[10],
  frag.s.slope           =  coef(frag.size.mod.NL)[2],
  frag.s.slope2          =  coef(frag.size.mod.NL)[3] + coef(frag.size.mod.NL)[12],
  # Variance in fragment size
  frag.sd.int            =  coef(frag.size.sd.mod.NL)[1] + coef(frag.size.sd.mod.NL)[9],
  frag.sd.slope          =  coef(frag.size.sd.mod.NL)[2] + coef(frag.size.sd.mod.NL)[11],
  # Colony fecundity 
  fec.a                  =  coef(fec.mod3)[1],
  fec.b                  =  coef(fec.mod3)[2],
  # Number of recruits (this will allow the model to vary between observed limits)
  rec.max                =  max(recruits.counts[which(recruits.counts$Country == "Japan" &
                                                        recruits.counts$Ecoregion == "Tropical" &
                                                        recruits.counts$LHS == "ST"), "N"]),
  rec.min                =  min(recruits.counts[which(recruits.counts$Country == "Japan" &
                                                        recruits.counts$Ecoregion == "Tropical" &
                                                        recruits.counts$LHS == "ST"), "N"]),
  # recruit size coefficients
  rec.size.int           =  coefficients(recruit.size.mod2)[1] + coefficients(recruit.size.mod2)[2],
  rec.size.sd            =  sigma(recruit.size.mod2))

########################################################## Weedy ########

# Weedy (Australia, Subtropical)
coefs.W.A.S <- c(
  # Survival probability
  surv.int               =  survival.coefs[1] + survival.coefs[4],
  surv.slope             =  survival.coefs[2] + survival.coefs[8],
  # Growth transitions
  grow.int               =  growth.coefs[1] + growth.coefs[5],
  grow.slope             =  growth.coefs[2],
  grow.slope2            =  growth.coefs[3] + growth.coefs[9],
  # Growth variance
  grow.sd.int            =  growth.SD.coefs[1] + growth.SD.coefs[4],
  grow.sd.slope          =  growth.SD.coefs[2] + growth.SD.coefs[8],
  # Fragmentation probability
  frag.int               =  frag.coefs[1] + frag.coefs[5],
  frag.slope             =  frag.coefs[2],
  frag.slope2            =  frag.coefs[3] + frag.coefs[9],
  # Number of frags produced  
  no.frag.int            =  frag.no.coefs[1] + frag.no.coefs[4],
  no.frag.slope          =  frag.no.coefs[2] + frag.no.coefs[7],
  # Fragment size
  frag.s.int             =  coef(frag.size.mod.NL)[1] + coef(frag.size.mod.NL)[5],
  frag.s.slope           =  coef(frag.size.mod.NL)[2],
  frag.s.slope2          =  coef(frag.size.mod.NL)[3] + coef(frag.size.mod.NL)[8],
  # Variance in fragment size
  frag.sd.int            =  coef(frag.size.sd.mod.NL)[1] + coef(frag.size.sd.mod.NL)[4],
  frag.sd.slope          =  coef(frag.size.sd.mod.NL)[2] + coef(frag.size.sd.mod.NL)[7],
  # Colony fecundity 
  fec.a                  =  coef(fec.mod4)[1],
  fec.b                  =  coef(fec.mod4)[2],
  # Number of recruits (this will allow the model to vary between observed limits)
  rec.max                =  max(recruits.counts[which(recruits.counts$Country == "Australia" &
                                                        recruits.counts$Ecoregion == "Subtropical" &
                                                        recruits.counts$LHS == "W"), "N"]),
  rec.min                =  min(recruits.counts[which(recruits.counts$Country == "Australia" &
                                                        recruits.counts$Ecoregion == "Subtropical" &
                                                        recruits.counts$LHS == "W"), "N"]),
  # recruit size coefficients
  rec.size.int           =  coefficients(recruit.size.mod2)[1] + coefficients(recruit.size.mod2)[3],
  rec.size.sd            =  sigma(recruit.size.mod2))

# Weedy (Australia, Tropical)
coefs.W.A.T <- c(
  # Survival probability
  surv.int               =  survival.coefs[1] + survival.coefs[14],
  surv.slope             =  survival.coefs[2] + survival.coefs[19],
  # Growth transitions
  grow.int               =  growth.coefs[1] + growth.coefs[15],
  grow.slope             =  growth.coefs[2],
  grow.slope2            =  growth.coefs[3] + growth.coefs[20],
  # Growth variance
  grow.sd.int            =  growth.SD.coefs[1] + growth.SD.coefs[14],
  grow.sd.slope          =  growth.SD.coefs[2] + growth.SD.coefs[19],
  # Fragmentation probability
  frag.int               =  frag.coefs[1] + frag.coefs[15],
  frag.slope             =  frag.coefs[2],
  frag.slope2            =  frag.coefs[3] + frag.coefs[20],
  # Number of frags produced  
  no.frag.int            =  frag.no.coefs[1] + frag.no.coefs[4],
  no.frag.slope          =  frag.no.coefs[2] + frag.no.coefs[7],
  # Fragment size
  frag.s.int             =  coef(frag.size.mod.NL)[1] + coef(frag.size.mod.NL)[5],
  frag.s.slope           =  coef(frag.size.mod.NL)[2],
  frag.s.slope2          =  coef(frag.size.mod.NL)[3] + coef(frag.size.mod.NL)[8],
  # Variance in fragment size
  frag.sd.int            =  coef(frag.size.sd.mod.NL)[1] + coef(frag.size.sd.mod.NL)[4],
  frag.sd.slope          =  coef(frag.size.sd.mod.NL)[2] + coef(frag.size.sd.mod.NL)[7],
  # Colony fecundity 
  fec.a                  =  coef(fec.mod4)[1],
  fec.b                  =  coef(fec.mod4)[2],
  # Number of recruits (this will allow the model to vary between observed limits)
  rec.max                =  max(recruits.counts[which(recruits.counts$Country == "Australia" &
                                                        recruits.counts$Ecoregion == "Tropical" &
                                                        recruits.counts$LHS == "W"), "N"]),
  rec.min                =  min(recruits.counts[which(recruits.counts$Country == "Australia" &
                                                        recruits.counts$Ecoregion == "Tropical" &
                                                        recruits.counts$LHS == "W"), "N"]),
  # recruit size coefficients
  rec.size.int           =  coefficients(recruit.size.mod2)[1] + coefficients(recruit.size.mod2)[3],
  rec.size.sd            =  sigma(recruit.size.mod2))

# Weedy (Japan, Subtropical)
coefs.W.J.S <- c(
  # Survival probability
  surv.int               =  survival.coefs[1] + survival.coefs[11],
  surv.slope             =  survival.coefs[2] + survival.coefs[17],
  # Growth transitions
  grow.int               =  growth.coefs[1] + growth.coefs[12],
  grow.slope             =  growth.coefs[2],
  grow.slope2            =  growth.coefs[3] + growth.coefs[18],
  # Growth variance
  grow.sd.int            =  growth.SD.coefs[1] + growth.SD.coefs[11],
  grow.sd.slope          =  growth.SD.coefs[2] + growth.SD.coefs[17],
  # Fragmentation probability
  frag.int               =  frag.coefs[1] + frag.coefs[12],
  frag.slope             =  frag.coefs[2],
  frag.slope2            =  frag.coefs[3] + frag.coefs[18],
  # Number of frags produced  
  no.frag.int            =  frag.no.coefs[1] + frag.no.coefs[10],
  no.frag.slope          =  frag.no.coefs[2] + frag.no.coefs[12],
  # Fragment size
  frag.s.int             =  coef(frag.size.mod.NL)[1] + coef(frag.size.mod.NL)[11],
  frag.s.slope           =  coef(frag.size.mod.NL)[2],
  frag.s.slope2          =  coef(frag.size.mod.NL)[3] + coef(frag.size.mod.NL)[13],
  # Variance in fragment size
  frag.sd.int            =  coef(frag.size.sd.mod.NL)[1] + coef(frag.size.sd.mod.NL)[10],
  frag.sd.slope          =  coef(frag.size.sd.mod.NL)[2] + coef(frag.size.sd.mod.NL)[12],
  # Colony fecundity 
  fec.a                  =  coef(fec.mod4)[1],
  fec.b                  =  coef(fec.mod4)[2],
  # Number of recruits (this will allow the model to vary between observed limits)
  rec.max                =  max(recruits.counts[which(recruits.counts$Country == "Japan" &
                                                        recruits.counts$Ecoregion == "Subtropical" &
                                                        recruits.counts$LHS == "W"), "N"]),
  rec.min                =  min(recruits.counts[which(recruits.counts$Country == "Japan" &
                                                        recruits.counts$Ecoregion == "Subtropical" &
                                                        recruits.counts$LHS == "W"), "N"]),
  # recruit size coefficients
  rec.size.int           =  coefficients(recruit.size.mod2)[1] + coefficients(recruit.size.mod2)[3],
  rec.size.sd            =  sigma(recruit.size.mod2))

# Weedy (Japan, Tropical)
coefs.W.J.T <- c(
  # Survival probability
  surv.int               =  survival.coefs[1] + survival.coefs[22],
  surv.slope             =  survival.coefs[2] + survival.coefs[24],
  # Growth transitions
  grow.int               =  growth.coefs[1] + growth.coefs[23],
  grow.slope             =  growth.coefs[2],
  grow.slope2            =  growth.coefs[3] + growth.coefs[25],
  # Growth variance
  grow.sd.int            =  growth.SD.coefs[1] + growth.SD.coefs[22],
  grow.sd.slope          =  growth.SD.coefs[2] + growth.SD.coefs[24],
  # Fragmentation probability
  frag.int               =  frag.coefs[1] + frag.coefs[23],
  frag.slope             =  frag.coefs[2],
  frag.slope2            =  frag.coefs[3] + frag.coefs[25],
  # Number of frags produced  
  no.frag.int            =  frag.no.coefs[1] + frag.no.coefs[10],
  no.frag.slope          =  frag.no.coefs[2] + frag.no.coefs[12],
  # Fragment size
  frag.s.int             =  coef(frag.size.mod.NL)[1] + coef(frag.size.mod.NL)[11],
  frag.s.slope           =  coef(frag.size.mod.NL)[2],
  frag.s.slope2          =  coef(frag.size.mod.NL)[3] + coef(frag.size.mod.NL)[13],
  # Variance in fragment size
  frag.sd.int            =  coef(frag.size.sd.mod.NL)[1] + coef(frag.size.sd.mod.NL)[10],
  frag.sd.slope          =  coef(frag.size.sd.mod.NL)[2] + coef(frag.size.sd.mod.NL)[12],
  # Colony fecundity 
  fec.a                  =  coef(fec.mod4)[1],
  fec.b                  =  coef(fec.mod4)[2],
  # Number of recruits (this will allow the model to vary between observed limits)
  rec.max                =  max(recruits.counts[which(recruits.counts$Country == "Japan" &
                                                        recruits.counts$Ecoregion == "Tropical" &
                                                        recruits.counts$LHS == "W"), "N"]),
  rec.min                =  min(recruits.counts[which(recruits.counts$Country == "Japan" &
                                                        recruits.counts$Ecoregion == "Tropical" &
                                                        recruits.counts$LHS == "W"), "N"]),
  # recruit size coefficients
  rec.size.int           =  coefficients(recruit.size.mod2)[1] + coefficients(recruit.size.mod2)[3],
  rec.size.sd            =  sigma(recruit.size.mod2))

########## Name and store variables

# Combine together all regression parameters
m.par <- matrix(c(
  C_A_S = coefs.C.A.S,
  C_A_T = coefs.C.A.T,
  C_J_S = coefs.C.J.S,
  C_J_T = coefs.C.J.T,
  ST_A_S = coefs.ST.A.S,
  ST_A_T = coefs.ST.A.T,
  ST_J_S = coefs.ST.J.S,
  ST_J_T = coefs.ST.J.T,
  W_A_S = coefs.W.A.S,
  W_A_T = coefs.W.A.T,
  W_J_S = coefs.W.J.S,
  W_J_T = coefs.W.J.T), nrow = 12, ncol = 23, byrow = T)

# and name columns and rows
rownames(m.par) <- c("C_A_S", "C_A_T", "C_J_S", "C_J_T", 
                     "ST_A_S", "ST_A_T", "ST_J_S", "ST_J_T",
                     "W_A_S", "W_A_T", "W_J_S", "W_J_T")

colnames(m.par) <-  c("surv.int", "surv.slope", # coefficients of survival
                      "grow.int", "grow.slope", "grow.slope2", # growth
                      "grow.sd.int", "grow.sd.slope", # variance in growth
                      "frag.int", "frag.slope", "frag.slope2", # fragmentation 
                      "frag.no.int", "frag.no.slope", # number of fragments
                      "frag.size.int", "frag.size.slope", "frag.size.slope2", # size of fragments
                      "frag.sd.int", "frag.sd.slope", # variance in fragment sizes
                      "fec.a", "fec.b", # colony fecundity
                      "Rmax", "Rmin", # recruit appearance
                      "rec.size", "rec.size.sd") # recruit sizes

# and extract relevant files for later!
write.csv(m.par, file = "vital rate parameters.csv")
write.csv(demographic.data , file = "Formatted data.csv" , row.names = FALSE)
write.csv(fecundity.data, file = "Fecundity data.csv", row.names = FALSE)

############################################################# End of Code#############################################################