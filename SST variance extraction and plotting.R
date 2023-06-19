# Script for extracting and plotting abiotic variability map across all survey sites using temperature as a representative abiotic variable 

# Author: James Cant (jic2@st-andrews.ac.uk)
# Date last modified: April 2021

# -----------------------------------------------------------------------------------------

# set working directory
setwd("Data_file_location")

# Load required packages
require(raster)
require(ncdf4)
source("Scale bar functions.R") # Functions derived from https://github.com/3wen/legendMap
require(reshape2)
require(tis)
require(ggplot2)
require(ggspatial)
require(maptools)
require(grid)
require(maps)
require(colorednoise)
require(stats)
require(viridis)

################################
# Step 1: Extract SST data
################################

# create a storage list
sst_data <- list()

# define key GPS coordinates
GPS_coords <- data.frame(Lat = c(-30.3, -23.4385315792, 32.76887, 26.5463),
                         Lon = c(153.143, 151.908396366, 132.63558, 128.086))
rownames(GPS_coords) <- c("AS","AT","JS","JT")

# define GPS boundaries for plotting
latmin <- as.numeric(-40.0)
latmax <- as.numeric(40.0)
lonmin <- as.numeric(120.0)
lonmax <- as.numeric(160.0)

# Generate series of month and year combinations for the time period over which I am interested in.
YearStart <- 1950
YearEnd <- 2019

# create time series datafile
date <- data.frame(Year = rep(YearStart:YearEnd, each = 12),
                   Month = rep(c("01","02","03","04","05","06","07","08","09","10","11","12"), length.out = 840)) 

# loop through available urls to extract temporal sst data
for (ii in 1:dim(date)[1]) {
  # progress read out
  print(ii)
  
  # define url
  the_url <- paste0("https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdHadISST.nc?sst%5B(",
                    date[ii, "Year"], "-", date[ii, "Month"], "-16T12:00:00Z)%5D%5B(", 
                    latmax, "):(", latmin, ")%5D%5B(", lonmin, "):(", lonmax, ")%5D&.draw=surface&.vars=longitude%7Clatitude%7Csst&.colorBar=%7C%7C%7C%7C%7C&.bgColor=0xffccccff", sep = "")
  # and define file name
  the_file <- paste0("SST_", date[ii, "Year"], "-", date[ii, "Month"], ".nc", sep = "")

  # download the required file
  download.file(the_url, the_file, mode='wb')
  
  # open the downloaded file
  rawSST <- nc_open(the_file)

  # now, grab stuff out of the nc file.
  temp <- list()
  temp$lats <- ncvar_get(rawSST, "latitude")
  temp$lon <- ncvar_get(rawSST, "longitude")
  temp$sst <- ncvar_get(rawSST, "sst")
  # Quick memory save
  nc_close(rawSST)
  
  # remove the downloaded nc file to save memory space
  file.remove(the_file)  
  
  # reshape the extracted sst data
  dimnames(temp$sst) <- list(long = temp$lon, lat = temp$lats)
  temp <- melt(temp$sst, value.name = "sst")
  
  # store extracted data
  sst_data[[ii]] <- temp
}

# Now to analyze the extracted sst data it needs reformatting out of a list
SST <- array(as.numeric(unlist(sst_data)), dim=c(3321, 3, 840))

# extract the mean temperature, cv, and variance for each pixel
SST_stats <- data.frame(Lat = SST[,2,1], # attach GPS
                        Lon = SST[,1,1],
                        # Mean temperature
                        SSTmean = RowMeans(SST[,3,], na.rm = TRUE),
                        # Coefficient of variation
                        CV  = apply(SST[,3,], 1, raster::cv, na.rm = TRUE),
                        # Autocorrelation
                        SSTAuto = rep(NA, dim(SST)[1]), # these will need completing in a loop below
                        # Variance frequency
                        SSTFreq = rep(NA, dim(SST)[1]))

# run a loop estimating the autocorrelation and frequency of monthly temperature series.
for(ii in 1:dim(SST_stats)[1]) { # repeat the loop for every location
  #progress read out
  print(ii)
  
  # reformat each locations temperature series as a time series
  try(temp_series <- ts(SST[ii,3,], start = 1, frequency = 1))
  
  # estimate the autocorrelation of the time series
  try(SST_stats$SSTAuto[ii] <- autocorrelation(temp_series, biasCorrection = TRUE))
  
  # determine the spectrum of the time series and extract the frequency and spectrum vectors
  try(temp_spectra <- spectrum(temp_series, plot = FALSE))
  try(spectra.x <- as.numeric(temp_spectra[["freq"]]))
  try(spectra.y <- as.numeric(temp_spectra[["spec"]]))
  # run a regression model between the two vectors
  try(spectra.mod <- lm(log(spectra.y)~log(spectra.x)))
  # the frequency of the environmental variance is the slope coefficient of this regression model
  try(SST_stats$SSTFreq[ii] <- as.numeric(coef(spectra.mod)[2]))
  
  # and tidy up
  rm(temp_series, temp_spectra, spectra.x, spectra.y, spectra.mod)
}

# quick tidy up
SST_stats[is.nan(SST_stats[,3]), "SSTmean"] <- NA
SST_stats[is.nan(SST_stats[,4]), "CV"] <- NA
SST_stats[is.nan(SST_stats[,5]), "SSTAuto"] <- NA
SST_stats[is.nan(SST_stats[,6]), "SSTFreq"] <- NA

# And write out the datafile
write.csv(SST_stats, file = "Temporal SST data.csv")

################################
# Step 2: Plot data for map figure
################################

# Load datafile
SST_stats <- read.csv("Temporal SST data.csv", row.names = 1)

# Source world map polygons
world <- map_data("world2", ".")

ggplot(data = SST_stats, aes(x = Lon, y = Lat, fill = SSTmean)) + 
  geom_raster(interpolate = TRUE) +
  labs(y = NULL, # add axis labels
       x = NULL) +
  scale_x_continuous(breaks = c(120,140,160),
                     labels = c("120.0°E","140.0°E","160.0°E"),
                     expand = c(0,0)) + 
  scale_y_continuous(breaks = c(-35,-20,0,20,35),
                     labels = c("35.0°S","20.0°S","0.0°","20.0°N","35.0°N"),
                     expand = c(0,0)) + 
  scale_fill_gradientn(colours = viridis(5, option = "B"),
                       na.value = "white",
                       guide = guide_colorbar(ticks = F, title = NULL, # these control the colour scale bar
                                              direction = "vertical", reverse = F, label = F,
                                              barwidth = 1, barheight = 15)) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), color = "gray50", fill = "gray50") +
  annotation_north_arrow(height = unit(1, "cm"), width = unit(0.6, "cm"), 
                         pad_x = unit(2.55, "cm"), pad_y = unit(4.0, "cm"),
                         style = north_arrow_orienteering(fill = c("black","black"))) + 
  scale_bar(lon = 125, lat = -30, distance_lon = 800, distance_lat = 200, distance_legend = 300, 
            legend_size = 4, dist_unit = "km", orientation = FALSE) + 
  coord_fixed(ratio = 1.3, xlim = c(lonmin, lonmax), ylim = c(latmin, latmax)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(size = 20, colour = "black"), axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
        legend.text = element_text(), legend.title = element_blank(), panel.border = element_rect(fill = NA))

# This map works best on the dimensions 978*796

################################
# Step 3: Extract key variables for the populations
################################

# Add space for variables to coordinate dataframe
GPS_coords$Mean <- GPS_coords$CV <- GPS_coords$Auto <- GPS_coords$Freq <- NA 

# run loop to extract required info.
for (ii in 1:dim(GPS_coords)[1]){
  
  # find closest matching GPS points and
  index <- which((abs(SST_stats$Lat - GPS_coords$Lat[ii]) == min(abs(SST_stats$Lat - GPS_coords$Lat[ii]))) & (abs(SST_stats$Lon - GPS_coords$Lon[ii]) == min(abs(SST_stats$Lon - GPS_coords$Lon[ii]))))
  
  # and extract
  GPS_coords[ii, "Mean"] <- SST_stats[index, "SSTmean"]
  GPS_coords[ii, "CV"] <- SST_stats[index, "CV"]
  GPS_coords[ii, "Auto"] <- SST_stats[index, "SSTAuto"]
  GPS_coords[ii, "Freq"] <- SST_stats[index, "SSTFreq"]
  
  # a little tidy up 
  rm(index)
}

# And write out the datafile
write.csv(GPS_coords, file = "Site SST info.csv")

####################################################### End of code #################################################

