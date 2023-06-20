# Corals-at-higher-latitudes
---
This repository contains R scripts detailing the analyses described in **Cant et al. (2023) Coral assemblages at higher latitudes favour short-term potential over long-term performance** published in Ecology.
For further information please contact me (James Cant) at *jic2@st-andrews.ac.uk*.

## Script details
---
Although users are welcome to download any number of the scripts contained within this repository, please be aware that each script has been designed for use in sequence (as displayed below), with each relying on outputs from previous scripts.

***Quantifying vital rate patterns:*** This script carries out the initial process of cleaning and reformating raw data describing the annual survival, and size transitions of tagged coral colonies. The raw data used for the manuscript associated with this code repository is hosted on *Dryad* and can be download at *doi:10.5061/dryad.w0vt4b8xd*. Following data cleaning, this script then estimates the survival, growth, fragmentation, and recruitment dynamics, of the focal coral assembalages before extracting all nessecary vital rate descriptors. 

***Constructing IPMs:*** Using the vital rate parameters extracted in the previous script, this R script parameterises a series of Integral Projection Models (IPMs) outlining the annual dynamics of the focal coral assemblages. This script also contains code for carrying out JackKnife resampling in order to generate a series of alternative IPMs to enable subsequent sensitivity analyses.

***SST variance extraction and plotting:*** Taking a series of GPS coordinates, this script accesses temperature records housed in the [NOAA environmental data repository] (https://coastwatch.pfeg.noaa.gov/), to calculate measures of monthly sea surface temperature variability (mean, coefficient of variance, autocorrelation, and variance frequency). This script also outlines the code needed for reproducing maps visualising spatial patterns in temperature variability.

***Evaluating spatial patterns:*** This script processes the IPMs parameterised using a previous script to evaluate for taxonomic and spatial variability in various estimates of transient (*short-term*) and asymptotic (*long-term*) demographic properties. 

Together the ***Computing sensitivities*** & ***Sensitivity analyses*** scripts are used for iteratively processing each JackKnifed IPM to assess the vital rate sensitivities of the transient and long-term demographic metrics that formed the basis of the manuscript associated with these code repository.

***IPM construction functions:*** This script contains a series of functions used in parameterising IPMs, and is called by the ***Sensitivity analyses*** script.

***Generation time analysis:*** After estimating a measure of generation time for each parameterised IPM this script works to evaluate how the generation time of coral assmeblages correlates with their their demographic characterisitics of long-term performance and resilience.
