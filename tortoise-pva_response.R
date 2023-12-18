############################################################################
############################################################################
##########    Modeling population viability of gopher tortoises   ##########
##########                across the species range                ##########
##########                B Folt, March 2022                      ##########
############################################################################
############################################################################

# Clear workspace
rm(list=ls())

# Set options
options(scipen=999)
set.seed(999)

###KJL specify which version of the code you want to run
this.version = 3 #choose 1-3 (see below)
version<-c("original","dd","3percent")[this.version] #specify the version of metapopulation growth
dir.create(version,showWarnings=F)#make output folder for this version
versionwd<-paste0("./",version)
print(sprintf("this version is: %s",version))

# Objective: use matrix population projection models to estimate
# population growth and extinction risk for gopher tortoises across the
# species range under different scenarios of climate and habitat change

# Methods:
# (i) estimate demographic parameters for tortoises (e.g., survival,fecundity, age of maturity, etc.)
#     and future changes in threats to tortoises, all which vary across the species range.
# (ii) model how future scenarios of climate and habitat change will influence tortoise demography,
#     population resiliency, and population redundancy
# (iii) visualize the results in tables and graphs

### Packages
# This project requires a number of R packages
#   Load the required packages; or, if they are not present in the user's default
#   library, they will be downloaded
packages = c("installr", "raster", "sp", "elevatr", "lognorm", "plyr", "plotrix",
             "clipr", "ggmap", "ggplot2", "maptools", "dplyr", "ggsn",
             "mapdata", "reshape", "ggrepel", "patchwork", "png") # Specify the packages of interest

package.check <- lapply( # Load or install+load all packages
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Rtools must also be installed: https://cran.rstudio.com/bin/windows/Rtools/
#   Check to see if Rtools is installed already
#install.Rtools()

# Install the package 'legendMap' using devtools
#devtools::install_github("3wen/legendMap")
library(legendMap)



##################################################################################
########## Section 1)                                                   ##########
########## Future scenarios of climate and habitat change               ##########
##################################################################################

# Use data describing tortoise population estimates throughout the species range
# and incorporate data describing scenarios of future climate and habitat change
datum = read.csv("tortoise-surveys.csv", header=TRUE)
head(datum)

# Change analysis units ("Units") to numeric
datum[datum$Units == "Unit 1",]$Units = 1
datum[datum$Units == "Unit 2",]$Units = 2
datum[datum$Units == "Unit 3",]$Units = 3
datum[datum$Units == "Unit 4",]$Units = 4
datum[datum$Units == "Unit 5",]$Units = 5
datum[,"Units"] = as.data.frame(apply(as.data.frame(datum[,"Units"]),2,as.numeric))
str(datum$Units)
summary(datum$Units)

# Plot data to make sure that representative genetic units do not overlap;
# the genetic populations identified by Gaillard et al. 2017 J. Fish Wildl. Management
plot(datum$Long, datum$Lat, pch=datum$Analysis.Unit, col=datum$Units)

# Round things out by changing the representative unit designations
# for a few instances
fix1 = row.names(subset(datum, datum$Long < -82.205 & datum$Long > -82.215)) #row 513
fix2 = row.names(subset(datum, datum$Units == "4" & datum$Long < -83.9)) #rows 438, 491
fix3 = row.names(subset(datum, datum$Units == "2" & datum$Long < -88)) #rows 54

datum[fix1,]$Units = 4 #change row to unit 4
datum[fix2[1],]$Units = 3 #change row to unit 3
datum[fix2[2],]$Units = 3 #change row to unit 3
datum[fix3,]$Units = 1 #change row to unit 1

# Plot data again to verify representative units do not overlap
plot(datum$Long, datum$Lat, pch=datum$Analysis.Unit, col=datum$Units)
#We seem to have done a decent job at lining these up w/ approximate
#genetic populations identified by Gaillard et al. 2017


### Some populations do not have any estimate of population size. Remove these
datum = droplevels(subset(datum, datum$Pop_Estimate_SSA > 0))

### Some populations do not have spatial GIS data. Remove these
datum = droplevels(subset(datum, datum$Ac_600m != "NA"))

### Save the cleaned datafile
saveRDS(datum, "clean-tortoise-data.rds")

##### Scenario 1
##### Global warming drives increases in mean annual temperatures

# Global warming is predicted to drive increase in mean annual temperature
# throughout Southeastern North America. To account for this, extract data for
# historical Mean Annual Temperatures (MAT) for each site, then model
# different increases in MAT over the next 100 years
# (e.g., 1.0, 2.0, and 4.0 degree increases)

# Extract WorldClim data for mean annual temp and precipitation
#   using the packages raster and sp
r = getData("worldclim", var="bio", res=10)
r = r[[c(1,12)]] #temperature, annual precipitation
names(r) = c("MAT_C","Precip_mm")
# There appears to be a way to use this package & the WorldClim data
# to predict future temperatures at sites under RCP 4.5 and 8.5 scenarios
# Might be more concise than the current manual functions below

# Identify coordinates for tortoise localities
pops = SpatialPoints(cbind(datum$Long, datum$Lat))
str(pops)

# Extract climate values for tortoise populations
values = extract(r,pops)

# Rescale temperature values
values[,1] = values[,1]/10

# Bind temperature and precipitation data to tortoise dataframe
MAT_C = values[,c("MAT_C")]
datum = cbind(datum, MAT_C)



##### Scenario 2
##### Sea-level rise may destroy habitat and inundate coastal populations in ocean

# Warming global temperatures will cause polar ice to melt and causing sea-level
# rise (SLR) worldwide, threatening the existence of coastal habitat. Evaluate the
# potential effects of SLR on tortoise populations by using NOAA sea-level rise
# scenarios of "intermediate-high", "high", and "extreme" SLR in the southeastern
# and observing how it may influence the existence of tortoise habitat

# Estimate elevation (meters) above sea leave (asl) for each tortoise population
ll_prj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #uses library(elevatr)
e = get_elev_point(pops, prj=ll_prj, src="aws", units="meters")
Elev_m = as.data.frame(e$elevation)

## NOAA projects sea-level rise to influence coastal areas broadly using different
## scenarios, such as the intermediate-high (IH), high (H), and extreme (EXT) at 40, 60, and 80 year projections.
## Using these estimates of rise, we can model how elevation of tortoise populations will decrease
## E.g., Fort Myers (FM) station
IH_40 = 0.78
IH_60 = 1.26
IH_80 = 1.83

H_40 = 1.23
H_60 = 1.71
H_80 = 2.55

EXT_40 = 1.23
EXT_60 = 2.08
EXT_80 = 3.16

# Bind elevation to sites and rename it
elevation = data.frame(c(as.data.frame(e$elevation)))
names(elevation) = c("Elev_m")
datum = cbind(datum, elevation)
head(datum)

## Sea-level rise will also disrupt connectivity among local populations within metapopulations.

# For each landscape population, we have the current area of habitat (datum$Acres_2_5k) and the
# predicted future area of habitat lost to SLR in 2060, 2080, and 2100 with three scenarios:
# intermediate-high, high, and extreme
head(datum$Acres_2_5k) # Current area

# E.g., NOAA SLR predictions for 2060
head(datum$NOAA_SLR_IntHigh_2060_ac) # Intermediate-high
head(datum$NOAA_SLR_High_2060_ac) # High
head(datum$NOAA_SLR_Ext_2060_ac)  # Extreme

# Add in zeroes for cells with NAs where no SLR effects are predicted
SLR = datum[,c("NOAA_SLR_IntHigh_2060_ac","NOAA_SLR_High_2060_ac","NOAA_SLR_Ext_2060_ac",
               "NOAA_SLR_IntHigh_2080_ac","NOAA_SLR_High_2080_ac","NOAA_SLR_Ext_2080_ac",
               "NOAA_SLR_IntHigh_2100_ac","NOAA_SLR_High_2100_ac","NOAA_SLR_Ext_2100_ac")] # Subset to columns of interest
SLR[is.na(SLR)] = 0 # Replace NAs with zeros
drops = c("NOAA_SLR_IntHigh_2060_ac","NOAA_SLR_High_2060_ac","NOAA_SLR_Ext_2060_ac",
          "NOAA_SLR_IntHigh_2080_ac","NOAA_SLR_High_2080_ac","NOAA_SLR_Ext_2080_ac",
          "NOAA_SLR_IntHigh_2100_ac","NOAA_SLR_High_2100_ac","NOAA_SLR_Ext_2100_ac") # Identify old columns to drop
datum = datum[, !(names(datum) %in% drops)] # Drop old columns with NAs
datum = cbind(datum,SLR) # Bind in columns with zeroes

# Estimate the proportion of landscape remaining after SLR in the future under the different scenarios
SLR_2060_IH = (datum$Acres_2_5k-datum$NOAA_SLR_IntHigh_2060_ac)/datum$Acres_2_5k
SLR_2060_H = (datum$Acres_2_5k-datum$NOAA_SLR_High_2060_ac)/datum$Acres_2_5k
SLR_2060_EXT = (datum$Acres_2_5k-datum$NOAA_SLR_Ext_2060_ac)/datum$Acres_2_5k
SLR_2080_IH = (datum$Acres_2_5k-datum$NOAA_SLR_IntHigh_2080_ac)/datum$Acres_2_5k
SLR_2080_H = (datum$Acres_2_5k-datum$NOAA_SLR_High_2080_ac)/datum$Acres_2_5k
SLR_2080_EXT = (datum$Acres_2_5k-datum$NOAA_SLR_Ext_2080_ac)/datum$Acres_2_5k
SLR_2100_IH = (datum$Acres_2_5k-datum$NOAA_SLR_IntHigh_2100_ac)/datum$Acres_2_5k
SLR_2100_H = (datum$Acres_2_5k-datum$NOAA_SLR_High_2100_ac)/datum$Acres_2_5k
SLR_2100_EXT = (datum$Acres_2_5k-datum$NOAA_SLR_Ext_2100_ac)/datum$Acres_2_5k

# Predicted rate of change in habitat loss to SLR
dSLR_2060_IH = SLR_2060_IH-1
dSLR_2060_H = SLR_2060_H-1
dSLR_2060_EXT = SLR_2060_EXT-1
dSLR_2080_IH = SLR_2080_IH-1
dSLR_2080_H = SLR_2080_H-1
dSLR_2080_EXT = SLR_2080_EXT-1
dSLR_2100_IH = SLR_2100_IH-1
dSLR_2100_H = SLR_2100_H-1
dSLR_2100_EXT = SLR_2100_EXT-1

# Replace predicted sea-level rise acreage in datafile with rates of habitat loss
drops = c("NOAA_SLR_IntHigh_2060_ac","NOAA_SLR_High_2060_ac","NOAA_SLR_Ext_2060_ac",
                 "NOAA_SLR_IntHigh_2080_ac","NOAA_SLR_High_2080_ac","NOAA_SLR_Ext_2080_ac",
                 "NOAA_SLR_IntHigh_2100_ac","NOAA_SLR_High_2100_ac","NOAA_SLR_Ext_2100_ac")
datum = datum[ , !(names(datum) %in% drops)]
datum = cbind(datum, dSLR_2060_IH, dSLR_2060_H, dSLR_2060_EXT,
              dSLR_2080_IH, dSLR_2080_H, dSLR_2080_EXT,
              dSLR_2100_IH, dSLR_2100_H, dSLR_2100_EXT)

##### Scenario 3
##### Land-use change may influence dispersal dynamics among populations

# Tortoise populations we have data on are from protected conservation lands, but
# changes in regional land-use/urbanization might alter properties adjacent to conservation lands
# in ways that might influence patterns of dispersal and metapopulation dynamics.
# We used the SLEUTH model (https://databasin.org/datasets/e5860ced8b4844e88431cdbefe425e1a) to
# estimate predicted urbanization at different probability thresholds for
# landscape meta-populations of gopher tortoises
head(datum)

## For each landscape population, we estimated the current area of urban habitat using SLEUTH current conditions
### and the predicted future area of urban habitat using SLEUTH in 2060, 2080, and 2100 at three probability
## thresholds: V200, V500, and V950
head(datum$Urban_SLEUTH_Current_ac) # Current urbanization around landscape populations

# E.g., 2060
head(datum$Urban_SLEUTH_2060_V200_ac) # Relaxed threshold (much urbanization)
head(datum$Urban_SLEUTH_2060_V500_ac) # Medium threshold (mean predicted urbanization)
head(datum$Urban_SLEUTH_2060_V950_ac) # Conservative (low predicted urbanization)

mean(na.omit(datum$Urban_SLEUTH_2060_V200_ac)) # Mean high urbanization threshold
mean(na.omit(datum$Urban_SLEUTH_2060_V500_ac)) # Mean intermediate urbanization threshold
mean(na.omit(datum$Urban_SLEUTH_2060_V950_ac)) # Mean low urbanization threshold

# Add in zeroes for cells with NAs where no urbanization effects are predicted
URB = datum[,c("Urban_SLEUTH_Current_ac",
               "Urban_SLEUTH_2060_V200_ac","Urban_SLEUTH_2060_V500_ac","Urban_SLEUTH_2060_V950_ac",
               "Urban_SLEUTH_2080_V200_ac","Urban_SLEUTH_2080_V500_ac","Urban_SLEUTH_2080_V950_ac",
               "Urban_SLEUTH_2100_V200_ac","Urban_SLEUTH_2100_V500_ac","Urban_SLEUTH_2100_V950_ac")] # Subset to columns of interest
URB[is.na(URB)] = 0 # Replace NAs with zeros
drops = c("Urban_SLEUTH_Current_ac",
          "Urban_SLEUTH_2060_V200_ac","Urban_SLEUTH_2060_V500_ac","Urban_SLEUTH_2060_V950_ac",
          "Urban_SLEUTH_2080_V200_ac","Urban_SLEUTH_2080_V500_ac","Urban_SLEUTH_2080_V950_ac",
          "Urban_SLEUTH_2100_V200_ac","Urban_SLEUTH_2100_V500_ac","Urban_SLEUTH_2100_V950_ac") # Identify old columns to drop
datum = datum[, !(names(datum) %in% drops)] # Drop old columns with NAs
datum = cbind(datum,URB) # Bind in columns with zeroes

# Manipulate the SLEUTH predictions so that they appropriately summarize
# predicted changes in urbanization for landscape populations for each timestep and threshold

# Current (initial) proportion of dispersal (non-urban) habitat for all landscape populations
PDHi = (datum$Acres_2_5k-datum$Urban_SLEUTH_Current_ac)/datum$Acres_2_5k

# Future proportion of dispersal (non-urban) habitat in the future
PDH_2060_200 = (datum$Acres_2_5k-datum$Urban_SLEUTH_2060_V200_ac)/datum$Acres_2_5k
PDH_2060_500 = (datum$Acres_2_5k-datum$Urban_SLEUTH_2060_V500_ac)/datum$Acres_2_5k
PDH_2060_950 = (datum$Acres_2_5k-datum$Urban_SLEUTH_2060_V950_ac)/datum$Acres_2_5k

PDH_2080_200 = (datum$Acres_2_5k-datum$Urban_SLEUTH_2080_V200_ac)/datum$Acres_2_5k
PDH_2080_500 = (datum$Acres_2_5k-datum$Urban_SLEUTH_2080_V500_ac)/datum$Acres_2_5k
PDH_2080_950 = (datum$Acres_2_5k-datum$Urban_SLEUTH_2080_V950_ac)/datum$Acres_2_5k

PDH_2100_200 = (datum$Acres_2_5k-datum$Urban_SLEUTH_2100_V200_ac)/datum$Acres_2_5k
PDH_2100_500 = (datum$Acres_2_5k-datum$Urban_SLEUTH_2100_V500_ac)/datum$Acres_2_5k
PDH_2100_950 = (datum$Acres_2_5k-datum$Urban_SLEUTH_2100_V950_ac)/datum$Acres_2_5k

# Predicted change in proportion of dispersal habitat in future
# delta PDH; dPDH
dPDH_2060_200 = PDH_2060_200-PDHi
dPDH_2060_500 = PDH_2060_500-PDHi
dPDH_2060_950 = PDH_2060_950-PDHi

dPDH_2080_200 = PDH_2080_200-PDHi
dPDH_2080_500 = PDH_2080_500-PDHi
dPDH_2080_950 = PDH_2080_950-PDHi

dPDH_2100_200 = PDH_2100_200-PDHi
dPDH_2100_500 = PDH_2100_500-PDHi
dPDH_2100_950 = PDH_2100_950-PDHi

# Replace predicted urbanization habitat loss acreage in datafile with rates of habitat loss
drops = c("Urban_SLEUTH_2060_V200_ac","Urban_SLEUTH_2060_V500_ac","Urban_SLEUTH_2060_V950_ac",
          "Urban_SLEUTH_2080_V200_ac","Urban_SLEUTH_2080_V500_ac","Urban_SLEUTH_2080_V950_ac",
          "Urban_SLEUTH_2100_V200_ac","Urban_SLEUTH_2100_V500_ac","Urban_SLEUTH_2100_V950_ac")
datum = datum[ , !(names(datum) %in% drops)]
datum = cbind(datum, dPDH_2060_200, dPDH_2060_500, dPDH_2060_950,
              dPDH_2080_200, dPDH_2080_500, dPDH_2080_950,
              dPDH_2100_200, dPDH_2100_500, dPDH_2100_950)

# Number of populations per metapopulation
PopsPerMetapop = matrix(NA, length(datum[,1]), 1)
colnames(PopsPerMetapop) = c("PopsPerMetapop")

for (q in 1:length(datum[,1])){
  id = datum[q,]$Landscape_pop_ID
  pops = subset(datum,datum$Landscape_pop_ID == id)
  PopsPerMetapop[q,] = length(pops[,1])
}

# Save these objects to datafiles
datum = cbind(datum, PDHi, PopsPerMetapop)

## Land managers of local populations may have a harder time using prescribed fire to
## manage tortoise habitat as urbanization encroaches on populations. We estimated the
## current and future projected distance of local populations to urban areas.
## Estimate the predicted change of urbanization distance for each local population.
head(datum)

# Current distance to urban area
head(datum$Urban_SLEUTH_Current_Distance_km)

# E.g., 2060
head(datum$Urban_SLEUTH_Distance_2060_V200_km) # distance to urban predicted by SLEUTH P=0.20
head(datum$Urban_SLEUTH_Distance_2060_V500_km) # distance to urban predicted by SLEUTH P=0.50
head(datum$Urban_SLEUTH_Distance_2060_V950_km) # distance to urban predicted by SLEUTH P=0.95

# Predicted change (delta; d) in Distance to Urban Area (dDUA)
dDUA_2060_200 = datum$Urban_SLEUTH_Distance_2060_V200_km-datum$Urban_SLEUTH_Current_Distance_km
dDUA_2060_500 = datum$Urban_SLEUTH_Distance_2060_V500_km-datum$Urban_SLEUTH_Current_Distance_km
dDUA_2060_950 = datum$Urban_SLEUTH_Distance_2060_V950_km-datum$Urban_SLEUTH_Current_Distance_km

dDUA_2080_200 = datum$Urban_SLEUTH_Distance_2080_V200_km-datum$Urban_SLEUTH_Current_Distance_km
dDUA_2080_500 = datum$Urban_SLEUTH_Distance_2080_V500_km-datum$Urban_SLEUTH_Current_Distance_km
dDUA_2080_950 = datum$Urban_SLEUTH_Distance_2080_V950_km-datum$Urban_SLEUTH_Current_Distance_km

dDUA_2100_200 = datum$Urban_SLEUTH_Distance_2100_V200_km-datum$Urban_SLEUTH_Current_Distance_km
dDUA_2100_500 = datum$Urban_SLEUTH_Distance_2100_V500_km-datum$Urban_SLEUTH_Current_Distance_km
dDUA_2100_950 = datum$Urban_SLEUTH_Distance_2100_V950_km-datum$Urban_SLEUTH_Current_Distance_km

datum = cbind(datum,
              dDUA_2060_200,dDUA_2060_500,dDUA_2060_950,
              dDUA_2080_200,dDUA_2080_500,dDUA_2080_950,
              dDUA_2100_200,dDUA_2100_500,dDUA_2100_950)


dev.off()
par(mfrow=c(2,2))#, oma=c(5,4,0,0)+0.5, mar=c(1,1,1,1)+0.5)
hist(datum$Urban_SLEUTH_Current_Distance_km, xlab="Distance to urban area (km)", ylab="No. of populations",main="NLCD current")
hist(datum$Urban_SLEUTH_Distance_2060_V500_km, xlab="Distance to urban area (km)", ylab="No. of populations",main="SLEUTH 2060 V500")
hist(datum$Urban_SLEUTH_Distance_2080_V500_km, xlab="Distance to urban area (km)", ylab="No. of populations",main="SLEUTH 2080 V500")
hist(datum$Urban_SLEUTH_Distance_2080_V500_km, xlab="Distance to urban area (km)", ylab="No. of populations",main="SLEUTH 2100 V500")



##### Scenario 4
##### Climate warming may decrease habitat management with prescribed fire
##### and overall habitat quality for tortoises in the longleaf pine ecosystem

# From meeting with land manager experts, a consensus emerged that land managers
# are worried about how climate warming will decrease their ability to use prescribed fire,
# due to increasing temperatures and decreased ability to safely use prescribed fire.
# Indeed, a recent paper (Kupfer et al. 2020 International Journal of Wildfire) found that
# the 'burn window' in the southeastern US is projected to decreases substantially during the
# summer, 'growing-season' burn conditions. According to this paper, to burn in the
# Southeastern US historically, 76.6% of days were available for burning during
# the winter (January-February), 80% of days were burnable in the spring transitionary period (March, April, May), and
# 64.5% of days were burnable in the summer (June-July).

# However, after considering two climate models (RCP 4.5 and RCP 8.5), those authors
# estimated that the number of burnable days will change. Specify those changes with objects:

# Winter - 76.6% historical:
#     RCP 4.5: increase to 78.2% by 2100
#     RCP 8.5: increase to 79.6% by 2100
BurnW_45 = 0.016
BurnW_45_SD = 0.0025

BurnW_85 = 0.03
BurnW_85_SD = 0.005

DaysW = 59

# Spring - 80.0% historical:
#     RCP 4.5: –4.0% (range: –5.9 to –0.6%) by 2100
#     RCP 8.5: –10.5% (range: –17.0 to –5.6%) by 2100
BurnSp_45 = -0.04
BurnSp_45_SD = 0.0075

BurnSp_85 = -0.105
BurnSp_85_SD = 0.02

DaysSp = 92

# Summer - 64.5% historical:
#     RCP 4.5: -23.9% by 2100
#     RCP 8.5: -43.6 by 2100
BurnSu_45 = -0.239
BurnSu_45_SD = 0.03

BurnSu_85 = -0.435
BurnSu_85_SD = 0.06

DaysSu = 61



########################################################################
########## Section 2)                                         ##########
########## Estimating demographic rates for populations       ##########
########################################################################

demorates = read.csv("demographic-rates.csv", header=TRUE)
head(demorates)

# Use WorldClim data from above to extract mean annual temperature and
# precipitation for each site with demographic data
r

# Identify coordinates for tortoise localities
sites = SpatialPoints(cbind(demorates$Long, demorates$Lat)) # uses packages raster and sp
str(sites)

# Extract climate values for tortoise populations
values = extract(r,sites)

# Rescale temperature values
values[,1] = values[,1]/10

# Bind temperature and precipitation data to tortoise dataframe
demorates = cbind(demorates,values[,c("MAT_C","Precip_mm")])
head(demorates)

#detach(package:raster); detach(package:sp)


##### 1A)
##### Age of sexual maturity

# Subset to age of sexual maturity parameters
maturityF = droplevels(subset(demorates, demorates$Parameter == "Maturity_F"))
maturityM = droplevels(subset(demorates, demorates$Parameter == "Maturity_M"))
maturity = droplevels(subset(rbind(maturityF, maturityM), Estimate != "NA"))
length(maturity$Estimate)

# Simple regression of age of maturity by mean-annual temperature; varying by sex
shapes = c(1,2)
shapes = shapes[as.numeric(as.factor(maturity$Parameter))]
plot(maturity$MAT_C, maturity$Estimate, pch=shapes,
     ylab="Age of sexual maturity", xlab="Mean annual temperature (C)")
legend("topright", inset=0.05, legend=c("Females","Males"),
       pch=c(1,2), cex=1)

# Build regression models estimating effects of MAT and sex on
# age of sexual maturity
mod1 = lm(Estimate ~ MAT_C, data=maturity)
summary(mod1); confint(mod1)
AIC(mod1)

mod2 = lm(Estimate ~ MAT_C + Parameter, data=maturity)
summary(mod2); confint(mod2)
AIC(mod2)

mod3 = lm(Estimate ~ MAT_C*Parameter, data=maturity)
summary(mod3); confint(mod3)
AIC(mod3)

## The most parsimonious model (mod2; lowest AIC and all parameters significant)
## supported MAT_C and sex effects on age of maturity

## Create linear model of just the female data
## and extract/save the intercept and slope
females = lm(Estimate ~ MAT_C, data=maturityF)
summary(females); confint(females)

Mint = summary(females)$coefficients[1,1]
MintSE = summary(females)$coefficients[1,2]
MintSD = MintSE*length(na.omit(maturityF$Estimate))^0.5
Mslope = summary(females)$coefficients[2,1]
MslopeSE = summary(females)$coefficients[2,2]
MslopeSD = MslopeSE*length(na.omit(maturityF$Estimate))^0.5
# For each deg C increase in MAT, female age of maturity (MA) decreased
# by 1.41 years (0.19-2.62, 95% CI), which was significant (P=0.029).
# Intercept = 43.5

# Check the assumptions of linear regression model:
#   continuous variables, linear, normally-distribution residuals,
#   homoschedasticity, autocorrelation
plot(females)
#   Continuous variables: both variables are continuous
#   Linearity: Residuals vs. fitted plot indicates a horizontal line, which is indicative of
#     a linear relationship
#   Normally-distributed residuals: the Normal Q-Q plot shows residuals on the straight-
#     dashed line, indicating normality of residuals
plot(maturityF$MAT_C[-c(5)],residuals(females)) # non-increasing trend, suggesting homoscedasticity
acf(residuals(females)[order(maturityF$MAT_C[-c(5)])]) # autocorrelation absent among data

## Age of maturity can be modeled with random draws based on MAT of population
## M - maturity age of adult females; e.g., at a site with MAT of 21 degrees C
mat = 21
muM = Mint + Mslope*mat
sdM = MslopeSD
rlnorm(100,getParmsLognormForMoments(muM,sdM^2)[1],getParmsLognormForMoments(muM,sdM^2)[2]) # uses package lognorm
hist(rlnorm(100,getParmsLognormForMoments(muM,sdM^2)[1],getParmsLognormForMoments(muM,sdM^2)[2]))


##### 1B)
##### Fecundity

# Analyze data to produce mean estimates of fecundity
# for tortoise populations' across the species range

fecundity = droplevels(subset(demorates, demorates$Parameter == "Fecundity"))
length(fecundity$Estimate)

# Explore fecundity data with box-and-whiskers plot of genetic populations
fecundity$Genetic.Population <- factor(fecundity$Genetic.Population,
                                       levels=c("Western", "Central", "West Georgia", "East Georgia", "Florida"))
plot(fecundity$Genetic.Population, fecundity$Estimate,
     ylab="Mean clutch size (no. eggs)", xlab="Genetic population")

# Regression plot
colors=c(1,2,3,4,5)
colors=colors[as.numeric(fecundity$Genetic.Population)]
plot(fecundity$MAT_C, fecundity$Estimate,
     ylab="Mean clutch size (no. eggs)", xlab="Mean annual temperature (C)",
     col=colors)
legend("topleft", legend=levels(fecundity$Genetic.Population),
       col=c(1,2,3,4,5,6), pch=1, cex=0.8)
mod = lm(fecundity$Estimate ~ fecundity$MAT_C)
summary(mod); confint(mod)
abline(mod)
# For each degree increase in MAT, population-level fecundity
# increases by 0.52 eggs/clutch; intercept = -4.38; P < 0.001

# Saving linear model parameters
Fint = summary(mod)$coefficients[1,1]
FintSE = summary(mod)$coefficients[1,2]
FintSD = FintSE*length(na.omit(fecundity$Estimate))^0.5
Fslope = summary(mod)$coefficients[2,1]
FslopeSE = summary(mod)$coefficients[2,2]
FslopeSD = FslopeSE*length(na.omit(fecundity$Estimate))^0.5

# Check the assumptions of linear regression model:
#   continuous variables, linear, normally-distribution residuals,
#   homoschedasticity, autocorrelation
plot(mod)
#   Continuous variables: both variables are continuous
#   Linearity: Residuals vs. fitted plot indicates a horizontal line, which is indicative of
#     a linear relationship
#   Normally-distributed residuals: the Normal Q-Q plot shows residuals on the straight-
#     dashed line, indicating normality of residuals
#   Homoscedasticity: Scale-Location graph is horizontal, indicating homosce
plot(fecundity$MAT_C,residuals(mod)) # non-increasing trend, suggesting homoscedasticity
acf(residuals(mod)[order(fecundity$MAT_C)]) # autocorrelation absent among data


# For example, at a site with MAT=21 deg C...
mat = 21
muF = Fint + Fslope*mat
sdF = FslopeSD
rlnorm(100,getParmsLognormForMoments(muF,sdF^2)[1],getParmsLognormForMoments(muF,sdF^2)[2])
hist(rlnorm(100,getParmsLognormForMoments(muF,sdF^2)[1],getParmsLognormForMoments(muF,sdF^2)[2]))


##### 1C)
##### Graph showing geographic variation in fecundity and maturity
dev.off()

# Script to save file as a PNG
png("figure4.png",         # File name
    width = 1150, height = 1400, # Width and height in inches
    bg = "white",          # Background color
    res = 200)             # Resolution
par(mfrow=c(2,1), oma=c(5,4,0,0)+0.5, mar=c(1,1,1,1)+0.5)

# Age of sexual maturity
plot(maturityF$MAT_C, maturityF$Estimate, lwd=2, axes=FALSE,
     cex.lab=1.4, cex.axis=2, ylab="Maturity age", xlab="",
     xlim=c(18,24))
females = lm(Estimate ~ MAT_C, data=maturityF)
wx = par("usr")[1:2]
new.x = seq(wx[1],wx[2],len=100)
pred = predict(females, new=data.frame(MAT_C=new.x), interval="conf")
lines(new.x,pred[,"fit"],lwd=2); lines(new.x,pred[,"lwr"],lty=3); lines(new.x,pred[,"upr"],lty=3)
axis(side=1, lwd=3, cex.axis=1.3)
axis(side=2, lwd=3, cex.axis=1.3)
mtext("Maturity age", side=2, line=3, cex=1.3, col="black")
text(23.5,19,"A", cex=3)
# For each deg C increase in MAT, female age of maturity (M) decreased
# by 1.41 years (0.19-2.62, 95% CI), which was significant (P=0.029).
# Intercept = 43.5
text(20,10,paste0("y = ",round(Mint,2)," + ",round(Mslope,2),"*MAT"), cex=1.3)
text(20,9,"P = 0.029", cex=1.3)

# Fecundity
plot(fecundity$MAT_C, fecundity$Estimate, lwd=2, axes=FALSE,
     cex.lab=1.4, cex.axis=2, ylab="Fecundity", xlab="",
     xlim=c(18,24), ylim=c(4,10))
fec = lm(Estimate ~ MAT_C, data=fecundity)
wx = par("usr")[1:2]
new.x = seq(wx[1],wx[2],len=100)
pred = predict(fec, new=data.frame(MAT_C=new.x), interval="conf")
lines(new.x,pred[,"fit"],lwd=2); lines(new.x,pred[,"lwr"],lty=3); lines(new.x,pred[,"upr"],lty=3)
axis(side=1, lwd=3, cex.axis=1.3)
axis(side=2, lwd=3, cex.axis=1.3)
mtext("Fecundity", side=2, line=3, cex=1.3, col="black")
mtext("Mean annual temperature (deg C)", side=1, line=3.5, cex=1.3, col="black")
text(23.5,9.6, "B", cex=3)
# For each degree increase in MAT, population-level fecundity
# increases by 0.52 eggs/clutch; intercept = -4.38; P < 0.001
text(20,9.5,paste0("y = ",round(Fint,2)," + ",round(Fslope,2),"*MAT"), cex=1.3)
text(20,9,"P < 0.001", cex=1.3)

# Closing the graphical device
dev.off()

# WorldClim automatically saves a folder called "Wc" that is pretty big (80 MB).
# Delete it.
unlink("wc10", recursive = TRUE)


#######################################################################
########## Section 3)                                        ##########
########## Construct algebraic population projection model   ##########
#######################################################################

### Specify datafile of tortoise populations for projection, tidy it all up,
### and specify other details for population projection
pops = droplevels(subset(datum, datum$Pop_Estimate_SSA != "NA")) # Remove NAs
pops = droplevels(subset(pops, pops$MAT_C != "NA"))
head(pops)
str(pops)
length(pops[,1]) #626, one population was lost
summary(pops$Pop_Estimate_SSA)

# Rename columns for tidyness
#head(pops)
pops = plyr::rename(pops, c("Site_Name"="SiteName", "Survey_Methodology"="SurveyMethod",
                      "Burrow_Number_SSA"="BurrowNumberSSA","Scope_Tortosie"="NumberTortoisesScoped",
                      "Conversion_Factor"="ConversionFactor","Pop_Estimate_SSA"="PopEst",
                      "Density_SSA"="Density","Pop_Est_LCL_SSA"="PopEstLCL","Pop_Est_UCL_SSA"="PopEstUCL",
                      "GT_Hab_ha_SSA"="Habitat_ha","Pop_Categ_SSA"="PopCateg","Notes_SSA"="Notes",
                      "Ac_600m"="Ac600m","HSI_Ac_600m"="HsiAc600m","Units"="GeneticUnits",
                      "Lat"="Latitude","Long"="Longitude","Landscape_pop_ID"="LandscapePopID"))
head(pops)


# Population estimates for many populations lack estimates of error; only a few
# have estimates of confidence intervals (95% CI). Use the information from observations
# with 95% CI to estimate the relationship between abundance and confidence limits
# use library(plotrix)
dev.off()
sub=subset(pops, pops$PopEstLCL != "NA")
plotCI(sub$PopEst, sub$PopEst, ui=sub$PopEstUCL, li=sub$PopEstLCL)
#as PopEst increases, 95% CI appears to increase as well

CI = sub$PopEstUCL-sub$PopEstLCL
plot(sub$PopEst, CI, xlab="Population estimate (N)", ylab="Confidence interval")

# A rule of thumb is that the breadth of 95% CI can be approximated as 4*SD
# Using this rule, first estimate SD for sites that have confidence limits
# by dividing the difference between UCL and LCL by 4
PopEstSD = (pops$PopEstUCL-pops$PopEstLCL)/4

# Next, for sites lacking confidence intervals, estimate the relationship
# between population size and SD and then predict SD, given population estimates
SDmodel = lm(CI/4 ~ sub$PopEst)
summary(SDmodel); confint(SDmodel)
#For each unit increase in estimated population size, SD increases by 0.178 units
SDint = summary(SDmodel)$coefficients[1,1]
SDslope = summary(SDmodel)$coefficients[2,1]
SDslopeSE = summary(SDmodel)$coefficients[2,2]

for (i in 1:length(PopEstSD)){
  if(is.na(PopEstSD[i]) == TRUE){ #if PopEstSD is in fact NA ("TRUE"), then
    PopEstSD[i] = SDint+SDslope*pops$PopEst[i]
  }
}

# Bind the 'PopEstSD' object to the populations object 'pops'
pops = cbind(pops, PopEstSD)

# Some populations have very small numbers of individuals are already functionally
#   extinct. Assuming a population with fewer than 3 adult females is functionally
#   extinct and that populations have 1:1 sex ratios and 3:1 adult to juvenile ratios,
#   then any population with fewer than 8 individuals would have less than 3 adult females
#   on average. Remove those from the dataset before performing population projections.
pops = droplevels(subset(pops, pops$PopEst >= 8))

# In the model, we will assume that if populations drop below sea level that they
#   go extinct. We will assess this by tracking how elevation (m asl) changes for populations
#   through time, after the sea level rises. Remove any populations that currently are
#   below sea level (< 0 m asl).
pops = droplevels(subset(pops, pops$Elev_m > 0))


### Summarize the number of populations and metapopulations by state and analysis unit
###   to report in the methods section
dim(pops)[1];length(table(pops$LandscapePopID)) # total number of pops & metapopulations

dim(subset(pops, pops$GeneticUnits == 1))[1]; length(table(subset(pops, pops$GeneticUnits==1)$LandscapePopID))#Genetic Unit 1
dim(subset(pops, pops$GeneticUnits == 2))[1]; length(table(subset(pops, pops$GeneticUnits==2)$LandscapePopID))#Genetic unit 2
dim(subset(pops, pops$GeneticUnits == 3))[1]; length(table(subset(pops, pops$GeneticUnits==3)$LandscapePopID))#Genetic unit 3
dim(subset(pops, pops$GeneticUnits == 4))[1]; length(table(subset(pops, pops$GeneticUnits==4)$LandscapePopID))#Genetic unit 4
dim(subset(pops, pops$GeneticUnits == 5))[1]; length(table(subset(pops, pops$GeneticUnits==5)$LandscapePopID))#Genetic unit 5

dim(subset(pops, pops$State == "FL"))[1]; length(table(subset(pops, pops$State=="FL")$LandscapePopID))#FL
dim(subset(pops, pops$State == "GA"))[1]; length(table(subset(pops, pops$State=="GA")$LandscapePopID))#GA
dim(subset(pops, pops$State == "MS"))[1]; length(table(subset(pops, pops$State=="MS")$LandscapePopID))#MS
dim(subset(pops, pops$State == "AL"))[1]; length(table(subset(pops, pops$State=="AL")$LandscapePopID))#AL
dim(subset(pops, pops$State == "SC"))[1]; length(table(subset(pops, pops$State=="SC")$LandscapePopID))#SC
dim(subset(pops, pops$State == "LA"))[1]; length(table(subset(pops, pops$State=="LA")$LandscapePopID))#LA




#######################################################################
################# Start a timer and run the for-loops #################
#######################################################################

ptm <- proc.time()

#### Run the model for a projection interval(s)
####    e.g., 80 years into the future

intervals = c(80)

for (u in 1:length(intervals)){

  # Define number of replicate simulations, time-period for population project,
  # and load packages for simulations
  r = 50             # Number of simulation replicates
  t = intervals[u]      # Number of years to project simulation replicates;
                        # NOT TO EXCEED 80!

  # Count the number of sites with a population to be projected
  sites = length(pops$PopEst) #

  # Create arrays to save simulation results for each population
  N = matrix(0,r,t)     # Population size (N) across replicates and time
  lam = matrix(0,r,t)   # Population growth
  Pext = matrix(0,r,t)  # Probability of quasi-extinction
  density = matrix(0,r,t) # Density (ind/ha) across replicates and time
  Nimm = matrix(0,r,t)	# Number of immigrants through reps and time
  #KJL
  LAM<-matrix(1,r,t) # metapopulation growth rate for across reps and time (was previously not stored)
  Ni = matrix(0, sites)     # Mean initial population size
  Njt = matrix(0,sites)     # Final juvenile population size
  Nat = matrix(0,sites)     # Final adult population size
  Nt = matrix(0,sites)      # Final total population size at time t
  NtLCL = matrix(0,sites)   # LCL of final population size at time t
  NtUCL = matrix(0,sites)   # UCL of final population size at time t
  MedDensity = matrix(0,sites) # Median density at time t
  meanlam = matrix(0,sites) # Average among yearly averages of lambda
  meanF = matrix(0,sites)   # Mean fecundity
  meanTja = matrix(0,sites) # Mean age of sexual maturity
  PE = matrix(0,sites)      # Mean probability of extinction

  Nmed = matrix(0,sites,t)  # Median population size among replicates
  Nlcl = matrix(0,sites,t)  # LCL population size among replicates
  Nucl = matrix(0,sites,t)  # UCL population size among replicates


  ##### Future conditions are characterized by uncertainty in demographic rates and threats
  ##### Let's consider different future scenarios to predict how tortoises might respond,
  ##### while accounting for uncertain future conditions in demography and threats.

  # 1) Let's first consider a baseline population scenario ('status quo') with adult
  #   female survival (0.96), no future climate warming, no future sea level rise,
  #   no future urbanization, the same habitat management regime as is currently used,
  #   immigration rate of 0.01, and density dependent limits on recruitment at or
  #   above 2 females/ha.

  # Previous demographic studies of gopher tortoises have found that populations are
  #   most sensitive to adult female survival, but there is uncertainty in what true
  #   survival rates are within and among populations.
  # Let's consider three additional scenarios, with adult female survival as
  # 1) 0.98,
  # 2) 0.94, and
  # 3) 0.92.

  # Some areas (particularly in peninsular Florida) have remarkably high densities of
  #   gopher tortoises in populations. Consider two scenarios where density-dependent
  #   limits on recruitment occur at:
  # 1) higher (4 females/ha) and
  # 2) lower densities (1 females/ha)

  # Climate change is driving increasing mean annual temperatures around the world.
  #   Tortoise maturation and reproductive rates currently increase along gradients
  #   of mean annual temperatures throughout the species range, and warming climates
  #   cause predictable shifts in those demographic rates.
  # Let's consider three warming scenarios:
  # 1) 1 deg C increase,
  # 2) 1.5 deg C increase,
  # 3) 2.0 deg C increase

  # Climate change is melting polar ice caps and driving sea-level rise around
  # the world. Tortoise populations in low-lying coastal areas may be at risk
  # of coastal inundation due to sea-level rise (SLR). Model three SLR scenarios from
  # the NOAA predictions:
  # 1) intermediate-high SLR,
  # 2) high SLR, and
  # 3) extreme SLR

  # While many tortoise populations occur on properties that are protected
  #   from future land-use change (i.e., no development), lands adjacent to tortoise
  #   populations may be subject to changing land-use and/or urbanization.
  #   Urbanization might disrupt dispersal dynamics among local populations within
  #   landscape populations (i.e., metapopulations) of tortoises.
  # Model three scenarios of predicted urbanization effects on landscape populations
  #   using the SLEUTH model:
  # 1) Land is urbanized in the future only if it has 95% probability or greater of being urbanized (least urbanization)
  # 2) Land is urbanized with 50% probability or greater (moderate urbanization)
  # 3)Land is urbanized with 20% probability or greater (most urbanization)

  # Habitat managers anticipate increased difficulty in managing land for high-
  #   quality tortoise habitat with prescribed fire due to climate warming.
  # Model four scenario involving habitat management with fire:
  # 1) decreased fire use predicted by RCP 4.5 ('less management')
  # 2) decreased fire use given RCP 8.5 ('much less management')
  # 3) increased fire use (or alternatives) predicted by the opposite of RCP4.5 ('more management')

  # Telemetry studies of tortoises have found that older individuals will occasionally
  #   undergo long-distance dispersal events at distances that could connect
  #   separate local populations (>2 km). However, the frequency at which these movements
  #   occur and the frequency which the individual finds another population to live at is
  #   unknown. Eubanks et al. (2003) found 2% of adults to emigrate from populations.
  # With a baseline immigration rate of 0.01 (1%) of adults that leave a population annually
  #   and successfully make it to another population, given no barriers to dispersal in the habitat,
  #   consider three other scenarios:
  # 1) no successful immigration (immigration rate = 0),
  # 2) higher successful immigration (0.02), and
  # 3) very high immigration (0.04).

  # While the aforementioned scenarios are useful for evaluating the sensitivity of the model
  #   to specific demographic rates and threats, future conditions for gopher tortoises
  #   are likely to be a product of multiple synergystic threats, and we still have uncertainty
  #   in two important demographic rates, survival and immigration.
  # To this end, consider 11 additional scenarios that have synergystic threats and uncertainty demographic rates:
  # (1) low threats + same fire, (2) medium threats + same fire, (3) high threats + same fire,
  # (4) more fire + medium threats, (5) less fire + medium threats, (6) much less fire + medium threats,
  # (7) high survival + medium threats, (8) low survival + medium threats,
  # (9) very high immigration + medium threats, (10) high immigration + medium threats, (11) no immigration + medium threats

  # These altogether yield 32 scenarios to be simulated as either sensitivity analysis or for
  #   future prediction
  nscens = 32
  Scenarios = matrix(NA, nscens, 8, dimnames=list(1:nscens, c("Scenario","AdultSurvival","ClimateWarming",
                                        "SeaLevelRise","Urbanization","Management","Immigration",
                                        "DensityDependence")))
  Scenarios[1,] = c("Status quo", 0.96, 0, "None", "None", "Same", 0.01, 2)
  Scenarios[2,] = c("Survival (high)", 0.98, 0, "None", "None", "Same", 0.01, 2)
  Scenarios[3,] = c("Survival (low)", 0.94, 0, "None", "None", "Same", 0.01, 2)
  Scenarios[4,] = c("Survival (very low)", 0.92, 0, "None", "None", "Same", 0.01, 2)
  Scenarios[5,] = c("Max density (high)", 0.96, 0, "None", "None", "Same", 0.01, 4)
  Scenarios[6,] = c("Max density (low)", 0.96, 0, "None", "None", "Same", 0.01, 1)
  Scenarios[7,] = c("Immigration (very high)", 0.96, 0, "None", "None", "Same", 0.04, 2)
  Scenarios[8,] = c("Immigration (high)", 0.96, 0, "None", "None", "Same", 0.02, 2)
  Scenarios[9,] = c("Immigration (zero)", 0.96, 0, "None", "None", "Same", 0, 2)
  Scenarios[10,] = c("Climate warming (low)", 0.96, 1, "None", "None", "Same", 0.01, 2)
  Scenarios[11,] = c("Climate warming (medium)", 0.96, 1.5, "None", "None", "Same", 0.01, 2)
  Scenarios[12,] = c("Climate warming (high)", 0.96, 2, "None", "None", "Same", 0.01, 2)
  Scenarios[13,] = c("Sea-level rise (low)", 0.96, 0, "IntHigh", "None", "Same", 0.01, 2)
  Scenarios[14,] = c("Sea-level rise (medium)", 0.96, 0, "High", "None", "Same", 0.01, 2)
  Scenarios[15,] = c("Sea-level rise (high)", 0.96, 0, "Extreme", "None", "Same", 0.01, 2)
  Scenarios[16,] = c("Urbanization (low)", 0.96, 0, "None", "Low", "Same", 0.01, 2)
  Scenarios[17,] = c("Urbanization (medium)", 0.96, 0, "None", "Medium", "Same", 0.01, 2)
  Scenarios[18,] = c("Urbanization (high)", 0.96, 0, "None", "High", "Same", 0.01, 2)
  Scenarios[19,] = c("Management (high)", 0.96, 0, "None", "None", "More", 0.01, 2)
  Scenarios[20,] = c("Management (low)", 0.96, 0, "None", "None", "Less", 0.01, 2)
  Scenarios[21,] = c("Management (very low)", 0.96, 0, "None", "None", "MuchLess", 0.01, 2)
  Scenarios[22,] = c("Low threats", 0.96, 1, "IntHigh", "Low", "Same", 0.01, 2)
  Scenarios[23,] = c("Medium threats", 0.96, 1.5, "High", "Medium", "Same", 0.01, 2)
  Scenarios[24,] = c("High threats", 0.96, 2, "Extreme", "High", "Same", 0.01, 2)
  Scenarios[25,] = c("Management (high) + medium threats", 0.96, 1.5, "High", "Medium", "More", 0.01, 2)
  Scenarios[26,] = c("Management (low) + medium threats", 0.96, 1.5, "High", "Medium", "Less", 0.01, 2)
  Scenarios[27,] = c("Management (very low) + medium threats", 0.96, 1.5, "High", "Medium", "MuchLess", 0.01, 2)
  Scenarios[28,] = c("Survival (high) + medium threats", 0.98, 1.5, "High", "Medium", "Same", 0.01, 2)
  Scenarios[29,] = c("Survival (low) + medium threats", 0.94, 1.5, "High", "Medium", "Same", 0.01, 2)
  Scenarios[30,] = c("Immigration (very high) + medium threats", 0.96, 1.5, "High", "Medium", "Same", 0.04, 2)
  Scenarios[31,] = c("Immigration (high) + medium threats", 0.96, 1.5, "High", "Medium", "Same", 0.02, 2)
  Scenarios[32,] = c("Immigration (zero) + medium threats", 0.96, 1.5, "High", "Medium", "Same", 0, 2)
  Scenarios = as.data.frame(Scenarios)

  ### Create lists to save summaries describing results from each scenario
  NumbScenarios = dim(Scenarios)[1]
  ScenarioResults = vector(mode="list", length=NumbScenarios)

  # Create arrays and lists to save results summarizing the scenarios
  variables = 23
  ScenSumm = matrix(0, length(Scenarios[,1]), variables)
  colnames(ScenSumm) = c("PopsTi","PopsTt","PopsTtLCL","PopsTtUCL","PopsTDelta",
                        "NMetaI","NMetaT","NMetaTLCL","NMetaTUCL","NMetaTDelta",
                        "NToti", "NTott","NTottLCL","NTottUCL","NTottDelta",
                        "Extant","ExtantPerc","VLikelyExtant","VLikelyExtantPerc",
                        "LikelyExtant","LikelyExtantPerc","UnlikelyExtant","UnlikelyExtantPerc")
  # Extant = Pe > 0.05, VLikelyExtant = Pe 0.05-0.2, LikelyExtant = Pe 0.2-0.5, UnlikelyExtant = Pe > 0.5

  # Create additional arrays/lists to save results for populations, units and states
  ScenSummStates = list(ScenSumm,ScenSumm,ScenSumm,ScenSumm,ScenSumm,ScenSumm) #Six states
  ScenSummUnits = list(ScenSumm,ScenSumm,ScenSumm,ScenSumm,ScenSumm) #Five genetic analysis units

  N_median = vector(mode="list", length=length(Scenarios[,1])) # Lists to save median (CI) for abundance at each site within each scenario
  N_lcl = vector(mode="list", length=length(Scenarios[,1]))
  N_ucl = vector(mode="list", length=length(Scenarios[,1]))

  PopVars = vector(mode = "list", length = sites)
  ScenarioPopVars = vector(mode="list", length=length(Scenarios[,1]))

  # Establish a minimum number of females for a population to stay extant
  MinFemales = 3

  # Track the number of extant populations through time during scenario projections
  ExtinctionTracker = matrix(NA, length(pops[,1]), t)
  ExtantTracker = ExtinctionTracker
  ScenExtantTracker = list(matrix(NA, length(Scenarios[,1]), t),
                           matrix(NA, length(Scenarios[,1]), t),
                           matrix(NA, length(Scenarios[,1]), t))
  names(ScenExtantTracker) = c("ExtantMean","ExtantLCL","ExtantUCL")


  ### A) Future scenario loop
  for (g in 1:length(Scenarios[,1])){
  #KJL make an empty list to fill with results for saving
    outlist<-list()
    ### Create objects with starting demographic rates (and associated error)
    ### for study populations that account for uncertainty by subject
    ### demographic rates to potential future scenarios
    params = matrix(NA, nrow=sites, ncol=12)
    colnames(params) = c("propB","F","phiN","propH","propF","phiH",
                         "phiJ","phiA","gJF","nJ","nA","Elevation")

    params[,1] <- 0.974    # Proportion females breeding
    params[,2] <- Fint + Fslope*pops$MAT_C # Mean fecundity
    params[,3] <- 0.349   # Mean nest survival
    params[,4] <- 0.85    # Proportion of eggs that hatch
    params[,5] <- 0.5     # Probability of female eggs (sex ratio)
    params[,6] <- 0.13    # Hatchling survival
    params[,7] <- 0.75    # Juvenile survival
    params[,8] <- as.numeric(Scenarios$AdultSurvival[g])    # Female survival
    params[,9] <- Mint + Mslope*pops$MAT_C # Age of sexual maturity
    params[,10] <- round(pops$PopEst*0.25*0.5) # Initial abundance of juvenile females
    params[,11] <- round(pops$PopEst*0.75*0.5) # Initial abundance of adult females
    params[,12] <- pops$Elev_m # Initial elevation of the population

    paramsSD = matrix(NA, nrow=sites, ncol=12)
    colnames(paramsSD) = c("propB_SD","F_SD","phiN_SD","propH_SD","propF_SD","phiH_SD",
                           "phiJ_SD","phiA_SD","gJF_SD","nJ_SD","nA_SD","Elev_SD")

    paramsSD[,1] <- 0.02    # SD of prop females breeding
    paramsSD[,2] <- FslopeSD # SD of mean fecundity
    paramsSD[,3] <- 0.103   # SE of mean nest survival
    paramsSD[,4] <- 0.105   # SD of proportion of eggs that hatch
    paramsSD[,5] <- 0.04    # SD of proportion female (sex ratio)
    paramsSD[,6] <- 0.05    # SD of hatchling survival
    paramsSD[,7] <- 0.06    # SD of juvenile survival
    paramsSD[,8] <- 0.03    # SD of female survival
    paramsSD[,9] <- MslopeSE # SE of age of sexual maturity
    paramsSD[,10] <- pops$PopEstSD*0.25*0.5 # SD of initial abundance of juvenile females
    paramsSD[,11] <- pops$PopEstSD*0.75*0.5 # SD of initial abundance of adult females
    paramsSD[,12] <- NA

    # Tortoise habitat often relies on prescribed fire to create high-quality habitat
    # for tortoises. This might require habitat to be burned every 2-5 years.
    # Specify a burn interval (BurnInt) and the probability of a population being burned (BurnProb)
    BurnInt = 4 # e.g., 4 year burn interval
    BurnProb = 1/BurnInt

    # Hunter & Rostal in press: tortoises probability of staying decreases 0.027
    # each consecutive year after a prescribed fire
    FireEffectOnSurvival = -0.027
    FireEffectSD = 0.0025

    ### Specify threat levels for each scenario
    warming = as.numeric(Scenarios$ClimateWarming[g])
    slr = Scenarios$SeaLevelRise[g]
    urbanization = Scenarios$Urbanization[g]
    fire = Scenarios$Management[g]
    densityDependence = Scenarios$DensityDependence[g]

    ### Specify immigration rate
    immigration = as.numeric(Scenarios$Immigration[g])

    # Specify rates of habitat loss due to urbanization specific to
    # projection period (40, 60, 80 years) and urbanization level (low, medium, high).
    # There are two urbanization effects in a list:
    # the first list item is the change in the proportion of landscape habitat available for dispersal,
    # and the second list item is the change in distance to urban area
    if (t <= 40) {if (urbanization=="None"){URB = list(rep(0,length(pops$dPDH_2060_950)),rep(0,length(pops$dDUA_2060_950)))} else {
                  if (urbanization=="Low"){URB = list(pops$dPDH_2060_950,pops$dDUA_2060_950)} else {
                  if (urbanization=="Medium"){URB = list(pops$dPDH_2060_500,pops$dDUA_2060_500)} else {
                    URB = list(pops$dPDH_2060_200,pops$dDUA_2060_200)}}}} else {
    if (t <= 60) {if (urbanization=="None"){URB = list(rep(0,length(pops$dPDH_2080_950)),rep(0,length(pops$dDUA_2080_950)))} else {
                  if (urbanization=="Low"){URB = list(pops$dPDH_2080_950,pops$dDUA_2080_950)} else {
                  if (urbanization=="Medium"){URB = list(pops$dPDH_2080_500,pops$dDUA_2080_500)} else {
                    URB = list(pops$dPDH_2080_200,pops$dDUA_2080_200)}}}} else {
    if (t <= 80) {if (urbanization=="None"){URB = list(rep(0,length(pops$dPDH_2100_950)),rep(0,length(pops$dDUA_2100_950)))} else {
                  if (urbanization=="Low"){URB = list(pops$dPDH_2100_950,pops$dDUA_2100_950)} else {
                  if (urbanization=="Medium"){URB = list(pops$dPDH_2100_500,pops$dDUA_2100_500)} else {
                    URB = list(pops$dPDH_2100_200,pops$dDUA_2100_200)}}}}}}

    # Specify rates of habitat loss due to sea-level rise specific to
    # projection period (40, 60, 80 years) and SLR level (IntHigh, High, Ext)
    # for both local and landscape populations
    if (t <= 40) {if (slr=="None"){SLR=list(0,c(rep(0,length(pops$dSLR_2060_IH))))} else {
                  if (slr=="IntHigh"){SLR=list(IH_40,pops$dSLR_2060_IH)} else {
                  if (slr=="High"){SLR = list(H_40,pops$dSLR_2060_H)} else {
                    SLR=list(EXT_40,pops$dSLR_2060_EXT)}}}} else {
    if (t <= 60) {if (slr=="None"){SLR=list(0,c(rep(0,length(pops$dSLR_2080_IH))))} else {
                  if (slr=="IntHigh"){SLR=list(IH_60,pops$dSLR_2080_IH)} else {
                  if (slr=="High"){SLR = list(H_60,pops$dSLR_2080_H)} else {
                    SLR=list(EXT_60,pops$dSLR_2080_EXT)}}}} else {
    if (t <= 80) {if (slr=="None"){SLR=list(0,c(rep(0,length(pops$dSLR_2100_IH))))} else {
                  if (slr=="IntHigh"){SLR=list(IH_80,pops$dSLR_2100_IH)} else {
                  if (slr=="High"){SLR = list(H_80,pops$dSLR_2100_H)} else {
                    SLR=list(EXT_80,pops$dSLR_2100_EXT)}}}}}}
      # Object 'SLR' is an indexed estimate of SLR for local populations (absolute change over interval) and
      # landscape populations (as rate of change over the projection interval)

    # Specify changes in fire management due to different scenarios
    # predicted by Kupfer et al. 2020:
    # same (no changes in fire use), less (less fire predicted by RCP 4.5),
    # much less (predicted by RCP 8.5), and more (more fire, the opposite of effect predicted by RCP 4.5)
    if (fire=="Same"){FIRE=0} else {
      if (fire=="Less"){FIRE=4.5} else {
        if (fire=="MuchLess"){FIRE=8.5} else {
          FIRE=-4.5}}}


    #### B) Population-level loop
    for (h in c(1:sites)){    # Subset to research site, h

      #### Demographic parameters at sites

      ## Pfb - proportion of females breeding
      mPfb = params[h,1]
      varPfb = paramsSD[h,1]
      aPfb = mPfb*((mPfb*(1-mPfb)/varPfb^2)-1)
      bPfb = (1-mPfb)*((mPfb*(1-mPfb)/varPfb^2)-1)
      Pfb = matrix(rbeta(r*t,aPfb,bPfb),r,1)
      #hist(Pfb)
      APfbi = matrix(0,r,1)
      BPfbi = matrix(0,r,1)
      Pfbt = matrix(0,r,t)

      ## F - fecundity of adult females (clutch size)
      muF = params[h,2]
      sdF = paramsSD[h,2]
      #mean(matrix(round(rlnorm(5000,getParmsLognormForMoments(muF,sdF^2)[1],getParmsLognormForMoments(muF,sdF^2)[2])),r,1))
      #hist(round(rlnorm(5000,getParmsLognormForMoments(muF,sdF^2)[1],getParmsLognormForMoments(muF,sdF^2)[2])),r,1)
      Fi = matrix(0,r,1)
      Fit = matrix(0,r,t)
      eggs = matrix(0,r,t)

      ## Pns - probability of nest surviving from predators
      mPns = params[h,3]
      varPns = paramsSD[h,3]
      aPns = mPns*((mPns*(1-mPns)/varPns^2)-1)
      bPns = (1-mPns)*((mPns*(1-mPns)/varPns^2)-1)
      Pns = matrix(rbeta(r*t,aPns,bPns),r,1)
      #hist(Pns)
      APnsi = matrix(0,r,1)
      BPnsi = matrix(0,r,1)
      Pnst = matrix(0,r,t)

      ## Ph - probability of eggs hatching (hatching success)
      mPh = params[h,4]
      varPh = paramsSD[h,4]
      aPh = mPh*((mPh*(1-mPh)/varPh^2)-1)
      bPh = (1-mPh)*((mPh*(1-mPh)/varPh^2)-1)
      Ph = matrix(rbeta(r*t,aPh,bPh),r,1)
      #hist(Ph)
      APhi = matrix(0,r,1)
      BPhi = matrix(0,r,1)
      Pht = matrix(0,r,t)

      ## PropF - proportion of eggs that are female (i.e., sex ratio)
      mPropF = params[h,5]
      varPropF = paramsSD[h,5]
      aPropF = mPropF*((mPropF*(1-mPropF)/varPropF^2)-1)
      bPropF = (1-mPropF)*((mPropF*(1-mPropF)/varPropF^2)-1)
      PropF = matrix(rbeta(r*t,aPropF,bPropF),r,1)
      #hist(PropF)
      APropFi = matrix(0,r,1)
      BPropFi = matrix(0,r,1)
      PropFt = matrix(0,r,t)

      ## Sh - survival of hatchlings
      mSh = params[h,6]                 # Mean hatchling survival
      varSh = paramsSD[h,6]							# Variance of mean
      aSh = mSh*((mSh*(1-mSh)/(varSh^2))-1)
      bSh = (1-mSh)*((mSh*(1-mSh)/(varSh^2))-1)
      Shi = matrix(rbeta(r,aSh,bSh),r,1)	# Parametric uncertainty
      #hist(Shi)
      #SDmShi = matrix(rinvgauss(r,varSh^2,1),r,1)
      AShi = matrix(0,r,1)				# beta distribution shape parameters
      BShi = matrix(0,r,1)				# beta distribution shape parameters
      Sht = matrix(0,r,t)					# annual variation survival

      ## Sj - survival of juveniles
      mSj = params[h,7]
      varSj = paramsSD[h,7]
      aSj = mSj*((mSj*(1-mSj)/(varSj^2))-1)
      bSj = (1-mSj)*((mSj*(1-mSj)/(varSj^2))-1)
      Sji = matrix(rbeta(r,aSj,bSj),r,1)
      #hist(Sji)
      #SDmSji = matrix(rinvgauss(r,varSj^2,1),r,1)
      ASji = matrix(0,r,1)
      BSji = matrix(0,r,1)
      Sjt = matrix(0,r,t)

      ## Sa - survival of adults
      mSa = params[h,8]
      varSa = paramsSD[h,8]
      aSa = mSa*((mSa*(1-mSa)/(varSa^2))-1)
      bSa = (1-mSa)*((mSa*(1-mSa)/(varSa^2))-1)
      Sai = matrix(rbeta(r,aSa,bSa),r,1)
      #hist(Sai)
      #SDmSai = matrix(rinvgauss(r,varSa^2,1),r,1)
      ASai = matrix(0,r,1)
      BSai = matrix(0,r,1)
      Sat = matrix(0,r,t)

      ## gJF - transition from juvenile to adult
      muTja = params[h,9]
      sdTja = paramsSD[h,9]
      #mean(matrix(round(rlnorm(5000,getParmsLognormForMoments(muTja,sdTja^2)[1],getParmsLognormForMoments(muTja,sdTja^2)[2])),r,1))
      #hist(round(rlnorm(5000,getParmsLognormForMoments(muTja,sdTja^2)[1],getParmsLognormForMoments(muTja,sdTja^2)[2])),r,1)
      Tjai = matrix(0,r,1)
      Tjait = matrix(0,r,t)
      age = matrix(0,r,t)

      ## Nj - initial abundance of juvenile females
      muNj = params[h,10]
      sdNj = paramsSD[h,10]
      #mean(matrix(round(rlnorm(5000,getParmsLognormForMoments(muNj,sdNj^2)[1],getParmsLognormForMoments(muNj,sdNj^2)[2])),r,1))
      #hist(round(rlnorm(5000,getParmsLognormForMoments(muNj,sdNj^2)[1],getParmsLognormForMoments(muNj,sdNj^2)[2])),r,1)
      Nji = matrix(0,r,1)
      Nj = matrix(0,r,t)

      ## Na - initial abundance of adult females
      muNa = params[h,11]
      sdNa = paramsSD[h,11]
      #mean(matrix(round(rlnorm(5000,getParmsLognormForMoments(muNa,sdNa^2)[1],getParmsLognormForMoments(muNa,sdNa^2)[2])),r,1))
      #hist(round(rlnorm(5000,getParmsLognormForMoments(muNa,sdNa^2)[1],getParmsLognormForMoments(muNa,sdNa^2)[2])),r,1)
      Nai = matrix(0,r,1)
      Na = matrix(0,r,t)

      ## Nm - initial abundance of metapopulation (females only)
      ## Nm is the sum of all other local populations in the metapopulation;
      ## the focal local population is not included
      muNm = round(sum(subset(pops,pops$LandscapePopID == pops[h,]$LandscapePopID)$PopEst))/2
      sdNm = 0.1*muNm
      Nmi = matrix(0,r,1)
      Nm = matrix(0,r,t)

      ## Elevation - meters of elevation above sea level
      elev = params[h,12]
      elevi = matrix(0,r,1)
      elevit = matrix(0,r,t)

      ## Immigration
      mImm = immigration
      varImm = 0.001
      if (mImm > 0){  # Immigration = 0 can't be modeled with error
        aImm = mImm*((mImm*(1-mImm)/varImm^2)-1)
        bImm = (1-mImm)*((mImm*(1-mImm)/varImm^2)-1)
        Immi = matrix(rbeta(r*t,aImm,bImm),r,1)
        #hist(Imm)
        AImmi = matrix(0,r,1)
        BImmi = matrix(0,r,1)
        Immt = matrix(0,r,t)}
        else {Immt = matrix(0,r,t)
      }

      ## Prescribed fire
      mDaysW = DaysW*0.766
      varDaysW = mDaysW*0.05
      #round(rlnorm(1,getParmsLognormForMoments(mDaysW,varDaysW^2)[1],getParmsLognormForMoments(mDaysW,varDaysW^2)[2]))
      DaysWi = matrix(0,r,1)
      DaysWit = matrix(0,r,t)

      mDaysSp = DaysSp*0.8
      varDaysSp = mDaysSp*0.05
      DaysSpi = matrix(0,r,1)
      DaysSpit = matrix(0,r,t)

      mDaysSu = DaysSu*0.645
      varDaysSu = mDaysSu*0.0
      DaysSui = matrix(0,r,1)
      DaysSuit = matrix(0,r,t)

      # Total number of burn days through time
      TotalBurnDays = matrix(0,r,t)

      # Change in total number of burn days through time
      ChangeBurnDays = matrix(0,r,t)

      if (FIRE == 4.5){
        # Effect of climate on burn window in winter
        # Estimates from the Kupfer et al. (2020) paper
        mBurnW = BurnW_45
        varBurnW = BurnW_45_SD
        aBurnW = mBurnW*((mBurnW*(1-mBurnW)/varBurnW^2)-1)
        bBurnW = (1-mBurnW)*((mBurnW*(1-mBurnW)/varBurnW^2)-1)
        BurnW = matrix(rbeta(r*t,aBurnW,bBurnW),r,1)
        BurnWi = matrix(0,r,1)
        ABurnWi = matrix(0,r,1)
        BBurnWi = matrix(0,r,1)
        BurnWit = matrix(0,r,t)

        #Spring
        mBurnSp = BurnSp_45*-1 #multiple by -1 to make values positive
        varBurnSp = BurnSp_45_SD
        aBurnSp = mBurnSp*((mBurnSp*(1-mBurnSp)/varBurnSp^2)-1)
        bBurnSp = (1-mBurnSp)*((mBurnSp*(1-mBurnSp)/varBurnSp^2)-1)
        BurnSp = matrix(rbeta(r*t,aBurnSp,bBurnSp),r,1)
        BurnSpi = matrix(0,r,1)
        ABurnSpi = matrix(0,r,1)
        BBurnSpi = matrix(0,r,1)
        BurnSpit = matrix(0,r,t)

        #Summer
        mBurnSu = BurnSu_45*-1 #multiple by -1 to make values positive
        varBurnSu = BurnSu_45_SD
        aBurnSu = mBurnSu*((mBurnSu*(1-mBurnSu)/varBurnSu^2)-1)
        bBurnSu = (1-mBurnSu)*((mBurnSu*(1-mBurnSu)/varBurnSu^2)-1)
        BurnSu = matrix(rbeta(r*t,aBurnSu,bBurnSu),r,1)
        BurnSui = matrix(0,r,1)
        ABurnSui = matrix(0,r,1)
        BBurnSui = matrix(0,r,1)
        BurnSuit = matrix(0,r,t)
      }

      if (FIRE == 8.5){
        mBurnW = BurnW_85
        varBurnW = BurnW_85_SD
        aBurnW = mBurnW*((mBurnW*(1-mBurnW)/varBurnW^2)-1)
        bBurnW = (1-mBurnW)*((mBurnW*(1-mBurnW)/varBurnW^2)-1)
        BurnW = matrix(rbeta(r*t,aBurnW,bBurnW),r,1)
        BurnWi = matrix(0,r,1)
        ABurnWi = matrix(0,r,1)
        BBurnWi = matrix(0,r,1)
        BurnWit = matrix(0,r,t)

        mBurnSp = BurnSp_85*-1 #multiple by -1 to make values positive
        varBurnSp = BurnSp_85_SD
        aBurnSp = mBurnSp*((mBurnSp*(1-mBurnSp)/varBurnSp^2)-1)
        bBurnSp = (1-mBurnSp)*((mBurnSp*(1-mBurnSp)/varBurnSp^2)-1)
        BurnSp = matrix(rbeta(r*t,aBurnSp,bBurnSp),r,1)
        BurnSpi = matrix(0,r,1)
        ABurnSpi = matrix(0,r,1)
        BBurnSpi = matrix(0,r,1)
        BurnSpit = matrix(0,r,t)

        mBurnSu = BurnSu_85*-1 #multiple by -1 to make values positive
        varBurnSu = BurnSu_85_SD
        aBurnSu = mBurnSu*((mBurnSu*(1-mBurnSu)/varBurnSu^2)-1)
        bBurnSu = (1-mBurnSu)*((mBurnSu*(1-mBurnSu)/varBurnSu^2)-1)
        BurnSu = matrix(rbeta(r*t,aBurnSu,bBurnSu),r,1)
        BurnSui = matrix(0,r,1)
        ABurnSui = matrix(0,r,1)
        BBurnSui = matrix(0,r,1)
        BurnSuit = matrix(0,r,t)
      }

      if (FIRE == 0){
        mBurnW = 1 # no effect of climate on burn days in winter
        BurnW = matrix(1,r,1)
        BurnWi = matrix(0,r,1)
        ABurnWi = matrix(0,r,1)
        BBurnWi = matrix(0,r,1)
        BurnWit = matrix(0,r,t)

        mBurnSp = 1 # no effect of climate on burn days in spring
        BurnSp = matrix(1,r,1)
        BurnSpi = matrix(0,r,1)
        ABurnSpi = matrix(0,r,1)
        BBurnSpi = matrix(0,r,1)
        BurnSpit = matrix(0,r,t)

        mBurnSu = 1 # no effect of climate on fire burn days in summer
        BurnSu = matrix(1,r,1)
        BurnSui = matrix(0,r,1)
        ABurnSui = matrix(0,r,1)
        BBurnSui = matrix(0,r,1)
        BurnSuit = matrix(0,r,t)
      }

      if (FIRE == -4.5){
        # Effect of climate on burn window in winter
        mBurnW = BurnW_45
        varBurnW = BurnW_45_SD
        aBurnW = mBurnW*((mBurnW*(1-mBurnW)/varBurnW^2)-1)
        bBurnW = (1-mBurnW)*((mBurnW*(1-mBurnW)/varBurnW^2)-1)
        BurnW = matrix(rbeta(r*t,aBurnW,bBurnW),r,1)
        BurnWi = matrix(0,r,1)
        ABurnWi = matrix(0,r,1)
        BBurnWi = matrix(0,r,1)
        BurnWit = matrix(0,r,t)

        #Spring
        mBurnSp = BurnSp_45*-1 #multiple by -1 to make values positive
        varBurnSp = BurnSp_45_SD
        aBurnSp = mBurnSp*((mBurnSp*(1-mBurnSp)/varBurnSp^2)-1)
        bBurnSp = (1-mBurnSp)*((mBurnSp*(1-mBurnSp)/varBurnSp^2)-1)
        BurnSp = matrix(rbeta(r*t,aBurnSp,bBurnSp),r,1)
        BurnSpi = matrix(0,r,1)
        ABurnSpi = matrix(0,r,1)
        BBurnSpi = matrix(0,r,1)
        BurnSpit = matrix(0,r,t)

        #Summer
        mBurnSu = BurnSu_45*-1 #multiple by -1 to make values positive
        varBurnSu = BurnSu_45_SD
        aBurnSu = mBurnSu*((mBurnSu*(1-mBurnSu)/varBurnSu^2)-1)
        bBurnSu = (1-mBurnSu)*((mBurnSu*(1-mBurnSu)/varBurnSu^2)-1)
        BurnSu = matrix(rbeta(r*t,aBurnSu,bBurnSu),r,1)
        BurnSui = matrix(0,r,1)
        ABurnSui = matrix(0,r,1)
        BBurnSui = matrix(0,r,1)
        BurnSuit = matrix(0,r,t)
      }

      # Prescribed fire probability at sites through time
      mFireProb = BurnProb
      varFireProb = 0.015
      aFireProb = mFireProb*((mFireProb*(1-mFireProb)/varFireProb^2)-1)
      bFireProb = (1-mFireProb)*((mFireProb*(1-mFireProb)/varFireProb^2)-1)
      FireProb = matrix(rbeta(r*t,aFireProb,bFireProb),r,1)
      FireProbi = matrix(0,r,1)
      AFireProbi = matrix(0,r,1)
      BFireProbi = matrix(0,r,1)
      FireProbit = matrix(0,r,t)
      FireProbChangeit = matrix(0,r,t)
      FireProbUrbit = matrix(0,r,t)

      # Prescribed fire at sites (yes/no)
      Fire = matrix(0,r,t)

      # Fire effect on survival
      mFireEffect = FireEffectOnSurvival*-1 #make positive for simulation
      varFireEffect = FireEffectSD
      aFireEffect = mFireEffect*((mFireEffect*(1-mFireEffect)/varFireEffect^2)-1)
      bFireEffect = (1-mFireEffect)*((mFireEffect*(1-mFireEffect)/varFireEffect^2)-1)
      FireEffect = matrix(rbeta(r*t,aFireEffect,bFireEffect),r,1)
      FireEffecti = matrix(0,r,1)
      AFireEffecti = matrix(0,r,1)
      BFireEffecti = matrix(0,r,1)
      FireEffectij = matrix(0,r,t)


      #### C) Replication loop; draws replicate-level means for each demographic parameter'
      ####     to model parametric uncertainty of parameters
      for (i in 1:r){

        #### Site-level demographic rates

        # Probability of females breeding
        APfbi[i] = 100*Pfb[i]
        BPfbi[i] = 100*(1-Pfb[i])

        # Fecundity of adult females
        Fi[i] = rlnorm(1,getParmsLognormForMoments(muF,sdF^2)[1],
                       getParmsLognormForMoments(muF,sdF^2)[2])

        # Probability of nest survival from predation
        APnsi[i] = 100*Pns[i]
        BPnsi[i] = 100*(1-Pns[i])

        # Probability of eggs hatching (hatching success)
        APhi[i] = 100*Ph[i]
        BPhi[i] = 100*(1-Ph[i])

        # Survival of hatchling females
        AShi[i] = 100*Shi[i]
        BShi[i] = 100*(1-Shi[i])

        # Survival of juvenile females
        ASji[i] = 100*Sji[i]
        BSji[i] = 100*(1-Sji[i])

        # Survival of adult females
        ASai[i] = 100*Sai[i]
        BSai[i] = 100*(1-Sai[i])

        # Age of maturing from juvenile to adult female
        Tjai[i] = rlnorm(1,getParmsLognormForMoments(muTja,sdTja^2)[1],
                         getParmsLognormForMoments(muTja,sdTja^2)[2])

        # Initial juvenile abundance
        if(muNj > 0) Nj[i,1] = round(rlnorm(1,getParmsLognormForMoments(muNj,sdNj^2)[1],
                        getParmsLognormForMoments(muNj,sdNj^2)[2])) else Nj[i,1] = 0

        # Initial female abundance
        if (muNa > 0) Na[i,1] = round(rlnorm(1,getParmsLognormForMoments(muNa,sdNa^2)[1],
                        getParmsLognormForMoments(muNa,sdNa^2)[2])) else Na[i,1] = 0

        # Initial metapopulation abundance
        if (muNm > 0) Nm[i,1] = round(rlnorm(1,getParmsLognormForMoments(muNm,sdNm^2)[1],
                          getParmsLognormForMoments(muNm,sdNm^2)[2])) else Nm[i,1] = 0

        # Probability of immigration
        if (mImm > 0){
          AImmi[i] = 100*Immi[i]
          BImmi[i] = 100*(1-Immi[i])
        }

        ### Site-level abiotic features

        # Elevation
        elevi[i] = elev

        ### Prescribed fire features

        # Initial days available for winter burn
        DaysWit[i,1] = round(rlnorm(1,getParmsLognormForMoments(mDaysW,varDaysW^2)[1],
                                  getParmsLognormForMoments(mDaysW,varDaysW^2)[2]))
        # Initial days available for spring burn
        DaysSpit[i,1] = round(rlnorm(1,getParmsLognormForMoments(mDaysSp,varDaysSp^2)[1],
                                  getParmsLognormForMoments(mDaysSp,varDaysSp^2)[2]))
        # Initial days available for summer burn
        DaysSuit[i,1] = round(rlnorm(1,getParmsLognormForMoments(mDaysSu,varDaysSu^2)[1],
                                  getParmsLognormForMoments(mDaysSu,varDaysSu^2)[2]))

        # Climate-change effects on prescribed fire
        ABurnWi[i] = 100*BurnW[i]	#Winter
        BBurnWi[i] = 100*(1-BurnW[i])

        ABurnSpi[i] = 100*BurnSp[i]	#Spring
        BBurnSpi[i] = 100*(1-BurnSp[i])

        ABurnSui[i] = 100*BurnSu[i]	#Summer
        BBurnSui[i] = 100*(1-BurnSu[i])

        # Prescribed fire probabilities through time
        AFireProbi[i] = 100*FireProb[i]
        BFireProbi[i] = 100*(1-FireProb[i])

        # Fire effect on survival
        AFireEffecti[i] = 100*FireEffect[i]
        BFireEffecti[i] = 100*(1-FireEffect[i])


        #### Projection loop; drawing annual demographic rates in each year from 1:t
        ####  to simulate temporal stochasticity
        for(j in 1:t){

          ## Abiotic parameters
          elevit[i,j] = elevi[i]-SLR[[1]]*j/t

          ## Management parameters
          BurnWit[i,j] = rbeta(1,ABurnWi[i],BBurnWi[i])     #Effect of climate change on winter burn days
          BurnSpit[i,j] = rbeta(1,ABurnSpi[i],BBurnSpi[i])  #Effect on spring burn days
          BurnSuit[i,j] = rbeta(1,ABurnSui[i],BBurnSui[i])  #Effect on summer burn days

          if (j >= 2 & FIRE == 4.5){ # Fire management scenario of DECREASING fire
            DaysWit[i,j] = round(DaysWit[i,1])+round(DaysWit[i,1])*((BurnWit[i,j-1]*j)/t) # bit more fire in winter
            DaysSpit[i,j] = round(DaysSpit[i,1])+round(DaysSpit[i,1])*((-1*BurnSpit[i,j-1]*j)/t) #multiple by -1 to make negative effect in spring
            DaysSuit[i,j] = round(DaysSuit[i,1])+round(DaysSuit[i,1])*((-1*BurnSuit[i,j-1]*j)/t) #multiple by -1 to make negative effect in summer
          } else {
            if (j >= 2 & FIRE == 8.5){ # Fire management scenario of VERY DECREASING fire
              DaysWit[i,j] = round(DaysWit[i,1])+round(DaysWit[i,1])*((BurnWit[i,j-1]*j)/t) # more fire in winter
              DaysSpit[i,j] = round(DaysSpit[i,1])+round(DaysSpit[i,1])*((-1*BurnSpit[i,j-1]*j)/t) #multiple by -1 to make negative effect in spring
              DaysSuit[i,j] = round(DaysSuit[i,1])+round(DaysSuit[i,1])*((-1*BurnSuit[i,j-1]*j)/t) #multiple by -1 to make negative effect in summer
            } else {
              if (j >= 2 & FIRE == 0){ # Fire management scenario of NO CHANGE in fire
                DaysWit[i,j] = round(DaysWit[i,1]) #+((BurnWit[i,j-1]*j)/t))
                DaysSpit[i,j] = round(DaysSpit[i,1]) #+((-1*BurnSpit[i,j-1]*j)/t)) #multiple by -1 to make burn rates negative again
                DaysSuit[i,j] = round(DaysSuit[i,1]) #+((-1*BurnSuit[i,j-1]*j)/t)) #multiple by -1 to make burn rates negative again
              } else {
                if (j >= 2 & FIRE == -4.5) { # Fire management scenario of INCREASING fire
                  DaysWit[i,j] = round(DaysWit[i,1])+round(DaysWit[i,1])*((-1*BurnWit[i,j-1]*j)/t) # opposite of RCP4.5 effect
                  DaysSpit[i,j] = round(DaysSpit[i,1])+round(DaysSpit[i,1])*((BurnSpit[i,j-1]*j)/t) # opposite of RCP4.5 effect
                  DaysSuit[i,j] = round(DaysSuit[i,1])+round(DaysSuit[i,1])*((BurnSuit[i,j-1]*j)/t) #NOT multiplied -1 to keep burn effect positive
                }
              }
            }
          }

          # Total burn days
          TotalBurnDays[i,j] = round(DaysWit[i,j]) + round(DaysSpit[i,j]) + round(DaysSuit[i,j])

          # Change in total burn days
          ChangeBurnDays[i,j] = TotalBurnDays[i,j]/TotalBurnDays[i,1]

          # Probability of prescribed fire
          FireProbit[i,j] = rbeta(1,AFireProbi[i],BFireProbi[i])

          # Probability of prescribed fire, with climate-driven changes in burn window
          FireProbChangeit[i,j] = FireProbit[i,j]*ChangeBurnDays[i,j]

          # Probability of prescribed fire, with urbanization
          DistToUrban = pops$Urban_SLEUTH_Current_Distance_km[h]+(URB[[2]][h]*(j/t)) #
          if(DistToUrban <= 0){UrbanEffect=0} else{ #Urban effect on fire probability
            if(DistToUrban > 0 & DistToUrban < 3.2){UrbanEffect=DistToUrban/3.2} else {
              UrbanEffect=1
            }
          }
          FireProbUrbit[i,j] = FireProbChangeit[i,j]*UrbanEffect
          if(FireProbUrbit[i,j] >= 1) next
          # look at this

          # Prescribed fire (1 = yes; 0 = no)
          Fire[i,j] = rbinom(n=1, size=1, FireProbUrbit[i,j])

          # Fire effect on survival
          FireEffectij[i,j] = rbeta(1,AFireEffecti[i],BFireEffecti[i])

          ### Demographic parameters
          # Probability of breeding
          Pfbt[i,j] = rbeta(1,APfbi[i],BPfbi[i])

          # Fecundity; modeled as a function of climate change
          Fit[i,j] = rlnorm(1,getParmsLognormForMoments(Fi[i]+Fslope*warming*j/t,sdF^2)[1],
                      getParmsLognormForMoments(Fi[i]+Fslope*warming*j/t,sdF^2)[2])

          # Prob of nests surviving predation
          Pnst[i,j] = rbeta(1,APnsi[i],BPnsi[i])

          # Prob of egg hatching
          Pht[i,j] = rbeta(1,APhi[i],BPhi[i])

          # Proportion of eggs that are female
          PropFt[i,j] = rbeta(1,aPropF,bPropF)

          # Survival of hatchlings
          Sht[i,j] = rbeta(1,AShi[i],BShi[i])

          # Survival of juveniles
          Sjt[i,j] = rbeta(1,ASji[i],BSji[i])

          # Maturation rates
          Tjait[i,j] = 1/(rlnorm(1,getParmsLognormForMoments(Tjai[i]+Mslope*warming*j/t,sdTja^2)[1],
                      getParmsLognormForMoments(Tjai[i]+Mslope*warming*j/t,sdTja^2)[2])-1)
              # calculated by 1 / (age of maturity - 1)
              # subtract one from denominator to remove hatchling year from juv stage

          # Adult survival for first three years
          if (j <= 3){Sat[i,j] = rbeta(1,ASai[i],BSai[i])} else {

          # Adult survival for years 4 and beyond; modeled as a function of year since last burn
          if (Fire[i,j] == 1){Sat[i,j] = rbeta(1,ASai[i],BSai[i])} else {
            if (Fire[i,j-1] == 1){Sat[i,j] = rbeta(1,ASai[i],BSai[i]) - FireEffectij[i,j]} else {
              if (Fire[i,j-2] == 1){Sat[i,j] = rbeta(1,ASai[i],BSai[i]) - 2*FireEffectij[i,j]} else {
                 Sat[i,j] = rbeta(1,ASai[i],BSai[i]) - 3*FireEffectij[i,j]}}}}

          # Probability of immigration as a function of habitat
          if (mImm > 0) {
            Immt[i,j] = rbeta(1,AImmi[i],BImmi[i])*(PDHi[h]+(URB[[1]][h]*(j/t))+(SLR[[2]][h]*(j/t)))
            if (Immt[i,j] < 0){Immt[i,j] = 0}
          } else {
            Immt[i,j] = 0
          }
          # Apparent survival + immigration cannot exceed 1.0.
          #   If it does, make immigration probability = 1 - survival.
          #   This would assume not mortality, which is unlikely, but is a more realistic.
          if (Sat[i,j] + Immt[i,j] > 1) Immt[i,j] = 1 - Sat[i,j]

          # Effect of sea-level rise or low abundance; if elevation is <1 m asl
          #   force abundance and reproduction to zero (i.e., local extinction)
          if (elevit[i,j] < 0){density[i,j] = 0; Na[i,j] = 0; eggs[i,j]=0; Nj[i,j] = 0; Nimm[i,j] = 0} else {
            #if (Na[i,j] < MinFemales){density[i,j] = 0; Na[i,j] = 0; eggs[i,j]=0; Nj[i,j] = 0; Nimm[i,j] = 0} else {

              ### Demographic population projection
              if (j>1) density[i,j] = (Na[i,j-1]+Nj[i,j-1])/(pops$Ac600m[h]*0.4047) #density = torts/ha = abun / ac*0.4047

              # Immigrants into the local population
              if (j>2) {
                if(pops$PopsPerMetapop[h] > 1){
                  #Nimm[i,j] = round(((Nm[i,j-1]-(Na[i,j-1]+Nj[i,j-1]))/(pops$PopsPerMetapop[h]-1))*0.75*Immt[i,j])
                  # No. immigrants = [metapopulation size minus size of local population)/number local pops minus the local population] times
                  #   the proportion of adults in the local pop (assumed to be 0.75) times the immigration rate
                  immPool = Nm[i,j-1]-(Na[i,j-1]+Nj[i,j-1])
                  if (immPool < 0){immPool = 0}
                  Nimm[i,j] <- rbinom(n=1, size=round((immPool/(pops$PopsPerMetapop[h]-1))*0.75), Immt[i,j])

                } else {
                  #Nimm[i,j] = round(Nm[i,j-1]*0.75*Immt[i,j])
                  Nimm[i,j] <- rbinom(n=1, size=round(Nm[i,j-1]*0.75), Immt[i,j])

                }
              }

              # Adults
              #KJL - fix density dependence issue for 3percent and dd versions
              if(version !="original"){densityDependence<-as.numeric(densityDependence)}
                if (j>1){ if (density[i,j-1] < densityDependence){
                  Na[i,j] = round(Na[i,j-1]*Sat[i,j-1])+round(Nj[i,j-1]*Sjt[i,j-1]*Tjait[i,j-1])+round(Nimm[i,j-1])} else {
                  Na[i,j] = round(Na[i,j-1]*Sat[i,j-1])+round(Nj[i,j-1]*Sjt[i,j-1]*Tjait[i,j-1]*0)+round(Nimm[i,j-1]*0)}}
                  # Recruitment to adult state is density-dependent; if density exceeds a scenario-defined limit in females/ha, then
                  # recruitment from the juvenile state and immigration becomes zero.

              # If number of females <3, force to extinction
              #if (j>1){if(Na[i,j-1] < 3){Na[i,j] = 0}}


              # Meta-population size
              #KJL replace LAM with LAM[i,j] so we keep track of these through time
              if (j>1) LAM[i,j] = Na[i,j]/Na[i,j-1]
              #KJL cap growth at 3% annually for metapopulation if version = 3percent
              if(version == "3percent"){
                if (j>1) {if(!is.na(LAM[i,j]) & (LAM[i,j]>1.03)) {LAM[i,j]<-1.03} }
              }
              if (j>1){if(is.na(LAM[i,j]) == FALSE & is.infinite(LAM[i,j]) == FALSE){Nm[i,j] = round(Nm[i,j-1]*LAM[i,j])} else {
                if(is.na(LAM[i,j]) == TRUE){Nm[i,j] = round(Nm[i,j-1]*1)} else {Nm[i,j] = round(Nm[i,j-1]*1)}}}

              # No. of eggs produced per year, while accounting for sex ratio
              eggs[i,j] = round(round(round(round(round(Na[i,j]*Pfbt[i,j])*Pnst[i,j])*Fit[i,j])*Pht[i,j])*PropFt[i,j])
                # No. of eggs is a product of:
                # no. of adults*prob breeding*prob nest surviving*fecundity/clutch*prob hatching*proportion of female hatchlings

              # Juveniles
                if (j>1) Nj[i,j] = round(Nj[i,j-1]*Sjt[i,j-1]*(1-Tjait[i,j-1])) + round(eggs[i,j-1]*Sht[i,j-1])
                # No. of juveniles is a product of:
                # No. of juveniles that survive and don't mature + no. of eggs that hatched and survive for one year
            #} # End-if < minFemales
          } # End-if < 3m above sea level

          # Calculate population growth rate
          if(j>1) lam[i,j] = (Na[i,j]+Nj[i,j])/(Na[i,j-1]+Nj[i,j-1])

          # Calculate extinction risk, the proportion of replicates that end w/ <3 females
          if (Na[i,j] < MinFemales) Pext[i,j]=1 else Pext[i,j]=0

        } # Close projection loop

      } # Close replication loop

      ### KJL For each population * scenario, save individual replicate level data
      out<-data.frame(
        i=rep(1:r,80),
        j=floor(0:(80*r-1)/r)+1,
        scen=g,
        pop=h,
        Nj=as.vector(Nj),
        Na=as.vector(Na),
        Sa=as.vector(Sat),
        Sj=as.vector(Sjt),
        LAM=as.vector(LAM), #KJL: note, i changed this to be a matrix so each rep * timestep value is recorded
        Nm=as.vector(Nm),
        Nimm=as.vector(Nimm),
        density=as.vector(density),
        densitythresh=rep(densityDependence,length(as.vector(density)))
      )
    #KJL - make sure it's in ascending order of time and reps
      out<-out[with(out, order(i, j)),]
      outlist[[h]]<-out #save it to the outlist (should be 457 long at end, one for each pop)
      #KJL now save the output


      ### For each population, save pertinent results in lists
      Variables = list(Nj, Na, Nj+Na, Nimm, Nm, density,
                        Sht, Sjt, Sat, Tjai, Fit, eggs,
                        elevit, Fire, FireProbit, FireProbChangeit, FireProbUrbit,
                        ChangeBurnDays, TotalBurnDays, DaysWit, DaysSpit, DaysSuit,
                        DistToUrban)
      names(Variables) = c("Nj", "Na", "N", "Nimm", "Nm", "Density",
                           "Sh", "Sj", "Sa", "Tja", "Fec", "Eggs",
                           "Elev", "Fire", "FireProb", "FireProbChange", "FireProbUrb",
                           "ChangeBurnDays", "TotalBurnDays", "BurnDaysWi", "BurnDaysSp", "BurnDaysSu",
                           "DistToUrban")
      PopVars[[h]] = Variables
        # This creates a massive list object with 626 entries (populations);
        # for each population, there are >20 sub-objects with 100x80 entries.
        # Saving this file will eat up multiple gigabites of storage. So, don't do this.
        # Examine this object at the end to make sure all objects are behaving as they should,
        # but don't save it with the other summarized results.

      ### Calculate important parameters, and save them in arrays

      # Median initial population size
      Ni[h] = median(Nj[,1]) + median(Na[,1])

      # Mean fecundity
      meanF[h] = mean(Fi)

      # Mean age of sexual maturity
      meanTja[h] = mean(Tjai)

      # Total population size
      N = Na + Nj

      # Final median juvenile population size
      Njt[h] = median(Nj[,t])

      # Final median adult population size
      Nat[h] = median(Na[,t])

      # Final total population size at time t
      Nt[h] = median(Nj[,t])+median(Na[,t])
      NtLCL[h] = round(apply(as.data.frame(N[,t]),2,quantile,probs=c(0.075)))
      NtUCL[h] = round(apply(as.data.frame(N[,t]),2,quantile,probs=c(0.925)))

      # Median density at time t
      MedDensity[h] = (median(Nj[,t])+median(Na[,t]))/(pops$Ac600m[h]*0.4047)

      # Mean (yearly) lambda
      mlamda = apply(lam, 2, mean, na.rm=TRUE)

      # Overall average lambda among all years t
      meanlam[h] = mean(mlamda)

      # Extinction risk
      pe = apply(Pext,2,sum)/r
      PEt = pe[t]
      PE[h] = mean(Pext[,t])

      # Track the extinction likelihood for each population across any given year
      ExtinctionTracker[h,] = apply(Pext,2,sum)/r

      ## Create objects summarizing abundance with confidence intervals through time
      # Median and CI of total population size through time
      Nmed[h,] = round(apply(N,2,median))
      Nlcl[h,] = round(apply(N,2,quantile,probs=c(0.075)))  #Lower 85%
      Nucl[h,] = round(apply(N,2,quantile,probs=c(0.925)))  #Upper 85%

    } # Close site-specific demographic loop
    dir.create(sprintf("%s/output",version),showWarnings = F)
    #dir.create("output_matrix_Nm",showWarnings = F)
    #dir.create(sprintf("%s/output_problems_Nm",version),showWarnings=F)
    saveRDS(outlist,file=paste0(version,"/output/scen-",g,".rds"))

    ## Create objects summarizing results from each SITE
    ## Bind projection results to original 'pops' object, and rename columns
    popsH = cbind(pops, meanF, meanTja, Ni, Njt, Nat, Nt, NtLCL, NtUCL, MedDensity, meanlam, PE)  # pops object saved for loop 'h'
    popsH = rename(popsH, c("meanF"="MeanFecundity", "meanTja"="MeanAgeMaturity",
                          "Ni"="AbunInitial", "Njt"="AbunJuvProjected",
                          "Nat"="AbunAdProjected", "Nt"="AbunProjected", "NtLCL"="AbunProjectedLCL","NtUCL"="AbunProjectedUCL",
                          "MedDensity"="DensityMedian", "meanlam"="MeanLamda", "PE"="ProbExtinction"))

    # Save results of each scenario to ScenarioResults and ScenarioPopVars
    ScenarioResults[[g]] = popsH
    ScenarioPopVars[[g]] = PopVars

    ## Save median abundance and CI in a list indexed by SCENARIO
    N_median[[g]] = Nmed
    N_lcl[[g]] = Nlcl
    N_ucl[[g]] = Nucl

    ## Summarize results for each SCENARIO, bind them to 'scenParams' object,
    ## and save in 'ScenSumm' object by indexing row 'g' and columns for variables
    reps = 100 # number of replicate for binomial simulations below
    MetaPopsExtant = matrix(NA,summary(pops$LandscapePopID)[6],100)
    MedianMetaPopsExtant = matrix(NA,summary(pops$LandscapePopID)[6],3) #median, LCL, UCL

    for (w in 1:t){
      ScenExtantTracker$ExtantMean[g,w] = median(replicate(reps, sum(rbinom(length(popsH[,1]), size=1, prob=1-ExtinctionTracker[,w])))) #Median of 100 simulations randomly drawing whether populations survive (1-ProbExtinction) and summing how many survive
      ScenExtantTracker$ExtantLCL[g,w] = apply(data.frame(replicate(reps,sum(rbinom(length(popsH[,1]), size=1, prob=1-ExtinctionTracker[,w])))),2,quantile,probs=c(0.075)) #LCL of 1000 simulations randomly drawing whether populations survive
      ScenExtantTracker$ExtantUCL[g,w] = apply(data.frame(replicate(reps,sum(rbinom(length(popsH[,1]), size=1, prob=1-ExtinctionTracker[,w])))),2,quantile,probs=c(0.925)) #LCL of 1000 simulations randomly drawing whether populations survive
    }


    # Overall results
    ScenSumm[g,"PopsTi"] = dim(popsH)[1] #Total number of populations modeled
    ScenSumm[g,"PopsTt"] = median(replicate(reps, sum(rbinom(length(popsH[,1]), size=1, prob=1-popsH$ProbExtinction)))) #Median of 100 simulations randomly drawing whether populations survive (1-ProbExtinction) and summing how many survive
    ScenSumm[g,"PopsTtLCL"] = apply(data.frame(replicate(reps,sum(rbinom(length(popsH[,1]), size=1, prob=1-popsH$ProbExtinction)))),2,quantile,probs=c(0.075)) #LCL of 1000 simulations randomly drawing whether populations survive (1-ProbExtinction) and summing how many survive
    ScenSumm[g,"PopsTtUCL"] = apply(data.frame(replicate(reps,sum(rbinom(length(popsH[,1]), size=1, prob=1-popsH$ProbExtinction)))),2,quantile,probs=c(0.925)) #UCL of 1000 simulations randomly drawing whether populations survive (1-ProbExtinction) and summing how many survive
    ScenSumm[g,"PopsTDelta"] = ((ScenSumm[g,"PopsTt"]-ScenSumm[g,"PopsTi"])/ScenSumm[g,"PopsTi"])*100 # Percent change in number of populations during projection interval
    ScenSumm[g,"NToti"] = round(sum(popsH$AbunInitial)) #Total number of tortoises modeled
    ScenSumm[g,"NTott"] = round(sum(popsH$AbunProjected)) #Total number of tortoise projected after interval
    ScenSumm[g,"NTottLCL"] = round(sum(popsH$AbunProjectedLCL)) #LCL of total number of tortoises projected
    ScenSumm[g,"NTottUCL"] = round(sum(popsH$AbunProjectedUCL)) #UCL of total number of tortoises projected
    ScenSumm[g,"NTottDelta"] = ((ScenSumm[g,"NTott"]-ScenSumm[g,"NToti"])/ScenSumm[g,"NToti"])*100 # Percent change in number of tortoise over projection interval
    ScenSumm[g,"NMetaI"] = length(table(popsH$LandscapePopID)) #Total number of metapopulations modeled
    for (m in 1:summary(popsH$LandscapePopID)[6]){
      mpop = subset(popsH, popsH$LandscapePopID == m)
      MetaPopsExtant[m,] = replicate(reps, sum(rbinom(length(mpop[,1]), size=1, prob=1-mpop$ProbExtinction)))}
    MetaPopsExtant[MetaPopsExtant > 0] = 1
    for (m in 1:summary(popsH$LandscapePopID)[6]){
      MedianMetaPopsExtant[m,1] = median(MetaPopsExtant[m,]) #Median
      MedianMetaPopsExtant[m,2] = summary(MetaPopsExtant[m,])[2] #1Q
      MedianMetaPopsExtant[m,3] = summary(MetaPopsExtant[m,])[5]} #3Q
    ScenSumm[g,"NMetaT"] = sum(MedianMetaPopsExtant[,1]) #Median total number of metapopulations surviving
    ScenSumm[g,"NMetaTLCL"] = round(sum(MedianMetaPopsExtant[,2])) #1Q total number of metapopulations surviving
    ScenSumm[g,"NMetaTUCL"] = round(sum(MedianMetaPopsExtant[,3])) #3Q total number of metapopulations surviving
    ScenSumm[g,"NMetaTDelta"] = ((ScenSumm[g,"NMetaT"]-ScenSumm[g,"NMetaI"])/ScenSumm[g,"NMetaI"])*100 # Percent change in number of metapopulations
    ScenSumm[g,"Extant"] = length(subset(popsH, 1-popsH$ProbExtinction >= 0.95)[,1])
    ScenSumm[g,"VLikelyExtant"] = length(subset(popsH, 1-popsH$ProbExtinction < 0.949 & 1-popsH$ProbExtinction >= 0.80)[,1])
    ScenSumm[g,"LikelyExtant"] = length(subset(popsH, 1-popsH$ProbExtinction < 0.799 & 1-popsH$ProbExtinction >= 0.5)[,1])
    ScenSumm[g,"UnlikelyExtant"] = length(subset(popsH, 1-popsH$ProbExtinction < 0.5)[,1])
    ScenSumm[g,"ExtantPerc"] = ScenSumm[g,"Extant"]/ScenSumm[g,"PopsTi"]
    ScenSumm[g,"VLikelyExtantPerc"] = ScenSumm[g,"VLikelyExtant"]/ScenSumm[g,"PopsTi"]
    ScenSumm[g,"LikelyExtantPerc"] = ScenSumm[g,"LikelyExtant"]/ScenSumm[g,"PopsTi"]
    ScenSumm[g,"UnlikelyExtantPerc"] = ScenSumm[g,"UnlikelyExtant"]/ScenSumm[g,"PopsTi"]

    # Perform the same calculations as above, but for each genetic representation unit (i.e., Analysis Units)
    units = max(pops$GeneticUnits)
    for (u in 1:units){
      unit = subset(popsH, popsH$GeneticUnits == u)
      MetaPopsExtantUnit = matrix(NA,summary(unit$LandscapePopID)[6],100)
      MedianMetaPopsExtantUnit = matrix(NA,summary(unit$LandscapePopID)[6],3) #median, LCL, UCL

      ScenSummUnits[[u]][g,"PopsTi"] = dim(unit)[1]
      ScenSummUnits[[u]][g,"PopsTt"] = round(sum(1-unit$ProbExtinction))
      ScenSummUnits[[u]][g,"PopsTtLCL"] = apply(data.frame(replicate(1000,sum(rbinom(length(unit[,1]), size=1, prob=1-unit$ProbExtinction)))),2,quantile,probs=c(0.075)) #LCL of 1000 simulations randomly drawing whether populations survive (1-ProbExtinction) and summing how many survive
      ScenSummUnits[[u]][g,"PopsTtUCL"] = apply(data.frame(replicate(1000,sum(rbinom(length(unit[,1]), size=1, prob=1-unit$ProbExtinction)))),2,quantile,probs=c(0.925)) #UCL of 1000 simulations randomly drawing whether populations survive (1-ProbExtinction) and summing how many survive
      ScenSummUnits[[u]][g,"PopsTDelta"] = ((ScenSummUnits[[u]][g,"PopsTt"]-ScenSummUnits[[u]][g,"PopsTi"])/ScenSummUnits[[u]][g,"PopsTi"])*100
      ScenSummUnits[[u]][g,"NMetaI"] = length(table(unit$LandscapePopID)) #Total number of metapopulations modeled
      for (m in 1:summary(unit$LandscapePopID)[6]){
        mpop = subset(unit, unit$LandscapePopID == m)
        MetaPopsExtantUnit[m,] = replicate(reps, sum(rbinom(length(mpop[,1]), size=1, prob=1-mpop$ProbExtinction)))}
      MetaPopsExtantUnit[MetaPopsExtantUnit > 0] = 1
      for (m in 1:summary(unit$LandscapePopID)[6]){
        MedianMetaPopsExtantUnit[m,1] = median(MetaPopsExtantUnit[m,]) #Median
        MedianMetaPopsExtantUnit[m,2] = summary(MetaPopsExtantUnit[m,])[2] #1Q
        MedianMetaPopsExtantUnit[m,3] = summary(MetaPopsExtantUnit[m,])[5]} #3Q
      ScenSummUnits[[u]][g,"NMetaT"] = sum(MedianMetaPopsExtantUnit[,1]) #Median total number of metapopulations surviving
      ScenSummUnits[[u]][g,"NMetaTLCL"] = sum(MedianMetaPopsExtantUnit[,2]) #1Q total number of metapopulations surviving
      ScenSummUnits[[u]][g,"NMetaTUCL"] = sum(MedianMetaPopsExtantUnit[,3]) #3Q total number of metapopulations surviving
      ScenSummUnits[[u]][g,"NMetaTDelta"] = ((ScenSummUnits[[u]][g,"NMetaT"]-ScenSummUnits[[u]][g,"NMetaI"])/ScenSummUnits[[u]][g,"NMetaI"])*100 #change in number of metapopulations
      ScenSummUnits[[u]][g,"NToti"] = round(sum(unit$AbunInitial))
      ScenSummUnits[[u]][g,"NTott"] = round(sum(unit$AbunProjected))
      ScenSummUnits[[u]][g,"NTottLCL"] = round(sum(unit$AbunProjectedLCL)) #LCL of total number of tortoises projected
      ScenSummUnits[[u]][g,"NTottUCL"] = round(sum(unit$AbunProjectedUCL)) #UCL of total number of tortoises projected
      ScenSummUnits[[u]][g,"NTottDelta"] = ((ScenSummUnits[[u]][g,"NTott"]-ScenSummUnits[[u]][g,"NToti"])/ScenSummUnits[[u]][g,"NToti"])*100
      ScenSummUnits[[u]][g,"Extant"] = length(subset(unit, 1-unit$ProbExtinction > 0.95)[,1])
      ScenSummUnits[[u]][g,"VLikelyExtant"] = length(subset(unit, 1-unit$ProbExtinction < 0.949 & 1-unit$ProbExtinction >= 0.80)[,1])
      ScenSummUnits[[u]][g,"LikelyExtant"] = length(subset(unit, 1-unit$ProbExtinction < 0.799 & 1-unit$ProbExtinction >= 0.5)[,1])
      ScenSummUnits[[u]][g,"UnlikelyExtant"] = length(subset(unit, 1-unit$ProbExtinction < 0.5)[,1])
      ScenSummUnits[[u]][g,"ExtantPerc"] = ScenSummUnits[[u]][g,"Extant"]/ScenSummUnits[[u]][g,"PopsTi"]
      ScenSummUnits[[u]][g,"VLikelyExtantPerc"] = ScenSummUnits[[u]][g,"VLikelyExtant"]/ScenSummUnits[[u]][g,"PopsTi"]
      ScenSummUnits[[u]][g,"LikelyExtantPerc"] = ScenSummUnits[[u]][g,"LikelyExtant"]/ScenSummUnits[[u]][g,"PopsTi"]
      ScenSummUnits[[u]][g,"UnlikelyExtantPerc"] = ScenSummUnits[[u]][g,"UnlikelyExtant"]/ScenSummUnits[[u]][g,"PopsTi"]
    }

    # Track progress looping through scenarios for
    if (g %in% 1:length(Scenarios[,1])){print((g/length(Scenarios[,1]))*100); print("%")}

  } # Close the scenario loop

# Save the results as a list
ScenSummTotal = cbind(Scenarios, ScenSumm)
row.names(ScenExtantTracker$ExtantMean) = c(ScenSummTotal[,1])
row.names(ScenExtantTracker$ExtantLCL) = c(ScenSummTotal[,1])
row.names(ScenExtantTracker$ExtantUCL) = c(ScenSummTotal[,1])
#RESULTS = list(ScenSummTotal, ScenSummStates, ScenSummUnits, ScenarioResults, ScenExtantTracker)
RESULTS = list(ScenSummTotal, ScenSummUnits, ScenExtantTracker, ScenarioResults)
names(RESULTS) = c("ScenSummTotal", "ScenSummUnits", "ScenPopsExtantTracker", "ScenarioResults")
names(RESULTS$ScenarioResults) = RESULTS$ScenSummTotal[,1]
#KJL save results to version folder
saveRDS(RESULTS, paste0(version,"/pva-output-results-",t,"yr.rds"))

#### Runtimes:
#### for t = 80, takes about ~2.1 hours

} # End interval
Tjait
# Stop the clock
proc.time() - ptm


####
### Examine summaries for all scenarios
# KJL read results from version folder
results = readRDS(paste0(version,"/pva-output-results-80yr.rds"))
ScenSummTotal = results[[1]]
ScenSummUnits = results[[2]]
ScenExtantTracker= results[[3]]
ScenarioResults= results[[4]]

## Examine summaries by genetic analysis unit
head(ScenSummUnits[[1]]) # 1 - LA+MS+AL
head(ScenSummUnits[[2]]) # 2 - AL/FL
head(ScenSummUnits[[3]]) # 3 - west GA
head(ScenSummUnits[[4]]) # 4 - east GA/SC
head(ScenSummUnits[[5]]) # 5 - peninsular FL

## Examine individual scenarios and their predictions for each population
## e.g, scenario 1 (low threats)
head(ScenarioResults[[1]],10)





#######################################################################
########## Section 4)                                        ##########
########## Examine the results of the population projection  ##########
########## using tables, figures, & appendices               ##########
#######################################################################

# Read in results of the population projections
# KJL read results from version folder
RESULTS80 = readRDS(paste0(version,"/pva-output-results-",80,"yr.rds"))
names(RESULTS80$ScenarioResults) = RESULTS80$ScenSummTotal[,1]

head(RESULTS80)
head(RESULTS80$ScenSummTotal) #Scenarios
head(RESULTS80$ScenSummUnits) #Genetic units
head(RESULTS80$ScenPopsExtantTracker) #Extinction tracker by scenario
head(RESULTS80$ScenarioResults) #Results for each population in each scenario

ScenSummTotal = RESULTS80[[1]]
ScenSummUnits = RESULTS80[[2]]
ScenExtantTracker= RESULTS80[[3]]
ScenarioResults= RESULTS80[[4]]


###############
### Table 2 ###
###############

### Create a table that lists the scenarios analyzed by the model
table2 = RESULTS80$ScenSummTotal[,1]

# Copy and paste table into Excel for more formal formatting
write_clip(table2)
#KJL save as csv file
write.csv(table2,file=paste0(version,"/table2.csv"),row.names=F)


###############
### Table 3 ###
###############

### Create a table that summarizes results of model predictions for each scenario
###   at 80 years into the future
table3 = RESULTS80$ScenSummTotal[,c("Scenario", "NToti","NTott","NTottLCL","NTottUCL","NTottDelta",
                           "PopsTi","PopsTt","PopsTtLCL","PopsTtUCL","PopsTDelta",
                           "NMetaI","NMetaT","NMetaTLCL","NMetaTUCL","NMetaTDelta")]
write_clip(table3)
#KJL save as csv file
write.csv(table3,file=paste0(version,"/table3.csv"),row.names=F)

###############
### Table 4 ###
###############

### Create a table that summarizes results of extinction risk
### for each population in each scenario
table4 = RESULTS80$ScenSummTotal[c(1:length(table2)),][,c("Scenario", "Extant","ExtantPerc","VLikelyExtant","VLikelyExtantPerc",
                        "LikelyExtant","LikelyExtantPerc","UnlikelyExtant","UnlikelyExtantPerc")]
write_clip(table4)
#KJL save as csv file
write.csv(table4,file=paste0(version,"/table4.csv"),row.names=F)


##################
### Appendix 2 ###
##################

### Summarize projection results for each genetic unit
### after an 80 YEAR PROJECTION ONLY
app2 = rbind(RESULTS80$ScenSummUnits[[1]][,c("NToti","NTott","NTottLCL","NTottUCL","NTottDelta",
                                    "PopsTi","PopsTt","PopsTtLCL","PopsTtUCL","PopsTDelta",
                                    "NMetaI","NMetaT","NMetaTLCL","NMetaTUCL","NMetaTDelta")],
             RESULTS80$ScenSummUnits[[2]][,c("NToti","NTott","NTottLCL","NTottUCL","NTottDelta",
                                    "PopsTi","PopsTt","PopsTtLCL","PopsTtUCL","PopsTDelta",
                                    "NMetaI","NMetaT","NMetaTLCL","NMetaTUCL","NMetaTDelta")],
             RESULTS80$ScenSummUnits[[3]][,c("NToti","NTott","NTottLCL","NTottUCL","NTottDelta",
                                    "PopsTi","PopsTt","PopsTtLCL","PopsTtUCL","PopsTDelta",
                                    "NMetaI","NMetaT","NMetaTLCL","NMetaTUCL","NMetaTDelta")],
             RESULTS80$ScenSummUnits[[4]][,c("NToti","NTott","NTottLCL","NTottUCL","NTottDelta",
                                    "PopsTi","PopsTt","PopsTtLCL","PopsTtUCL","PopsTDelta",
                                    "NMetaI","NMetaT","NMetaTLCL","NMetaTUCL","NMetaTDelta")],
             RESULTS80$ScenSummUnits[[5]][,c("NToti","NTott","NTottLCL","NTottUCL","NTottDelta",
                                    "PopsTi","PopsTt","PopsTtLCL","PopsTtUCL","PopsTDelta",
                                    "NMetaI","NMetaT","NMetaTLCL","NMetaTUCL","NMetaTDelta")])
app2 = data.frame(cbind(rep(table2,5),app2))
colnames(app2)[1] = "Scenarios"
app2 = subset(app2, app2$Scenarios %in% c("Low threats","Medium threats","High threats",
                                          "Management (high) + medium threats","Management (low) + medium threats",
                                          "Management (very low) + medium threats"))
write_clip(app2)
#KJL save as csv file
write.csv(app2,paste0(version,"/app2.csv"),row.names=F)


#################
### Figures 1 ###
#################

# Made manually in Microsoft PowerPoint


################
### Figure 2 ###
################

# Map illustrating initial population size of populations,
#   and sites with demographic data

### A map of each population modeled within each genetic analysis units
### and scaled by initial population size

# Load packages for mapping data, create simple maps, and then
# plot the tortoise sites onto the map
state_coords = read.csv("state-coords.csv", header=TRUE)
state_coords[4,c(2,3)] = c(-80.8,26.6)
colnames(state_coords) = c("State","Longitude","Latitude")

require(dplyr); require(RColorBrewer); require(ggplot2)
require(mapdata); require(maptools)

usa = map_data("usa")

(SEmap = ggplot() +
    geom_polygon(data = usa, aes(x=long, y = lat, group = group),
                 fill = "grey", color="black", lwd=1.1) +
    coord_fixed(xlim = c(-91, -80), ylim = c(25.5, 34), ratio = 1.2) +
    theme(line = element_blank(),
          text = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.border=element_rect(size=2, fill = NA)))

# Add gopher tortoise distribution polygon to map
#range <- data.frame(getKMLcoordinates("TortoiseRange.kml", ignoreAltitude=TRUE))
unit1 <- data.frame(getKMLcoordinates("Unit1.kml", ignoreAltitude=TRUE))
unit2 <- data.frame(getKMLcoordinates("Unit2.kml", ignoreAltitude=TRUE))
unit3 <- data.frame(getKMLcoordinates("Unit3.kml", ignoreAltitude=TRUE))
unit4 <- data.frame(getKMLcoordinates("Unit4.kml", ignoreAltitude=TRUE))
unit5 <- data.frame(getKMLcoordinates("Unit5.kml", ignoreAltitude=TRUE))
colnames(unit1) = c("Long","Lat")
colnames(unit2) = c("Long","Lat")
colnames(unit3) = c("Long","Lat")
colnames(unit4) = c("Long","Lat")
colnames(unit5) = c("Long","Lat")

torts = SEmap +
  geom_polygon(data=unit1, aes(x = Long, y = Lat, group = NA),
               fill = "lightgreen", lty=5, size=1.5, alpha=1/2) +
  geom_polygon(data=unit2, aes(x = Long, y = Lat, group = NA),
               fill = "lightgreen", lty=5, size=1.5, alpha=1/2) +
  geom_polygon(data=unit2, aes(x = Long, y = Lat, group = NA),
               fill = "skyblue", lty=5, size=1.5, alpha=1/2) +
  geom_polygon(data=unit3, aes(x = Long, y = Lat, group = NA),
               fill = "cyan", lty=5, size=1.5, alpha=1/2) +
  geom_polygon(data=unit4, aes(x = Long, y = Lat, group = NA),
               fill = "pink", lty=5, size=1.5, alpha=1/2) +
  geom_polygon(data=unit5, aes(x = Long, y = Lat, group = NA),
               fill = "orange", lty=5, size=1.5, alpha=1/2)
plot(torts)

# Plot the states where tortoises occur
states <- map_data("state")
tortoise_states <- subset(states, region %in% c("alabama","georgia","florida","louisiana","mississippi"))
(state_map <- torts + geom_polygon(data = tortoise_states, aes(x=long, y = lat, group=group), fill = NA, color="black", lwd=1.1))

# Read in functions to create a scale-bar, northing, etc.
source("functions.R")

# Using  above North America map as a base & create a map for the Gopher Tortoise
i = -90; j = 27 # Specify the coordinate position (i,j) of the legend
a = -83; b = 26; l = 0.75 # Specify  coordinate positions (a,b) and length (l) of northing bar
#a = -80; b = 31; l = 0.75 # Specify  coordinate positions (a,b) and length (l) of northing bar
z = 3 # Specify the line type for the custom line in the map

## Create a custom dataframe for the legends
lat = c(26.7,26.7)
long = c(-89.5,-88)
leg1 = data.frame(cbind(lat,long))

lat = c(27.2,27.2)
long = c(-89.5,-88)
leg2 = data.frame(cbind(lat,long))

# Add points for each population & color by genetic analysis unit
head(pops)
fig2 <- state_map +  geom_point(data=pops, aes(x=Longitude, y=Latitude),
                                colour="black", fill=pops$GeneticUnits+2,
                                shape=21, size=log(pops$PopEst))
plot(fig2)

# Add a legend, northing, etc. to the make the map more informative
fig2a <- fig2 +
  ## Label the states
  annotate(x=state_coords$Longitude, y=state_coords$Latitude,
           label=state_coords$State, colour="black", geom="text", size=6)	+
  ## Add a northing
  geom_segment(arrow=arrow(length=unit(5,"mm")),
               aes(x=a, xend=a, y=b, yend=b+l), lwd=2, colour="black") +
  annotate(x=a, y=b-0.4, label="N", colour="black", geom="text", size=10) +
  ## Add a scale bar
  ggsn::scalebar(x.min = -92, x.max = -84,
                 y.min = 25.5, y.max = 32,
                 dist = 100, dist_unit = "km",
                 st.bottom = FALSE, st.color = "black", st.dist	= 0.03,
                 transform = TRUE, model = "WGS84") +
  ## Add legend for points
  geom_point(aes(x=i-0.8, y=j+1.5), colour = "black", fill=3, shape=21, size=6) +
  geom_point(aes(x=i-0.8, y=j+1.0), colour = "black",  fill=4, shape=21, size=6) +
  geom_point(aes(x=i-0.8, y=j+0.5), colour = "black", fill=5, shape=21, size=6) +
  geom_point(aes(x=i-0.8, y=j+0), colour = "black", fill=6, shape=21, size=6) +
  geom_point(aes(x=i-0.8, y=j-0.5), colour = "black", fill=7, shape=21, size=6) +
  annotate(x=i, y=j+1.5, label="Western Genetic Unit", geom="text", size=6, hjust=0) +
  annotate(x=i, y=j+1, label="Central Genetic Unit", geom="text", size=6, hjust=0) +
  annotate(x=i, y=j+0.5, label="West Georgia Genetic Unit", geom="text", size=6, hjust=0) +
  annotate(x=i, y=j+0, label="East Georgia Genetic Unit", geom="text", size=6, hjust=0) +
  annotate(x=i, y=j-0.5, label="Florida Genetic Unit", geom="text", size=6, hjust=0)
plot(fig2a)

# Add points showing locations for fecundity and maturity age estimates
fecundity = subset(demorates, demorates$Parameter == "Fecundity")
matage = subset(demorates, demorates$Parameter == "Maturity_F")
(fig2b = fig2a +
    geom_point(aes(x=fecundity$Long, y=fecundity$Lat), colour="black",
               fill="grey", shape=22, size=4) +
    geom_point(aes(x=matage$Long, y=matage$Lat), colour="black",
               fill="lightgoldenrod1", shape=23, size=3) +
    geom_point(aes(x=i-0.8, y=j-1), colour="black", fill="grey", shape=22, size=4) +
    geom_point(aes(x=i-0.8, y=j-1.5), colour="black", fill="lightgoldenrod1", shape=23, size=3) +
    annotate(x=i, y=j-1, label="Fecundity estimate", geom="text", size=6, hjust=0) +
    annotate(x=i, y=j-1.5, label="Maturity age estimate", geom="text", size=6, hjust=0))

## Add a picture of Gopher Tortoise to the second map
## using library(png) and library(grid)
img <- readPNG("tortoise.png")
g <- rasterGrob(img, interpolate=TRUE)
fig2c <- fig2b + annotation_custom(g, xmin=-85.7, xmax=-83.7, ymax=27.5, ymin=26)

## Plot the two maps side by side
plot(fig2c)
#KJL Save in version folder
ggsave(paste0(version,"/figure2.png"), width=8, height=8)

################
### Figure 3 ###
################

# Made in Powerpoint


################
### Figure 4 ###
################

# Made above


################
### Figure 5 ###
################

t = 80 # number of years during projection
ExtantPops = melt(ScenExtantTracker$ExtantMean)
LCL = melt(ScenExtantTracker$ExtantLCL)
UCL = melt(ScenExtantTracker$ExtantUCL)
ExtantPops = cbind(ExtantPops,LCL[,3],UCL[,3])
colnames(ExtantPops) = c("Scenario","Year","PopsExtant","LCL","UCL")
ScenarioNumber = rep(1:length(table2), t)
ExtantPops = cbind(ScenarioNumber,ExtantPops)

ExtantPops$Label <- NA
ExtantPops$Label[which(ExtantPops$Year == max(ExtantPops$Year))] <- 1:length(table2)
ExtantPops$Label2 <- NA
ExtantPops$Label2 <-  rep(1:32, 80)
ExtantPops$Scenario2 <- paste(ExtantPops$Label2," = ", ExtantPops$Scenario)
ExtantPops$Scenario2 <- factor(ExtantPops$Scenario2,
                               levels=c('1  =  Status quo','2  =  Survival (high)','3  =  Survival (low)','4  =  Survival (very low)',
                                        '5  =  Max density (high)','6  =  Max density (low)','7  =  Immigration (very high)','8  =  Immigration (high)',
                                        '9  =  Immigration (zero)','10  =  Climate warming (low)','11  =  Climate warming (medium)','12  =  Climate warming (high)',
                                        '13  =  Sea-level rise (low)','14  =  Sea-level rise (medium)','15  =  Sea-level rise (high)','16  =  Urbanization (low)',
                                        '17  =  Urbanization (medium)','18  =  Urbanization (high)','19  =  Management (high)','20  =  Management (low)',
                                        '21  =  Management (very low)','22  =  Low threats','23  =  Medium threats','24  =  High threats','25  =  Management (high) + medium threats',
                                        '26  =  Management (low) + medium threats','27  =  Management (very low) + medium threats','28  =  Survival (high) + medium threats',
                                        '29  =  Survival (low) + medium threats','30  =  Immigration (very high) + medium threats',
                                        '31  =  Immigration (high) + medium threats','32  =  Immigration (zero) + medium threats'))
colnames(ExtantPops) = c("ScenarioNumber","Scenario1","Year","PopsExtant","LCL","UCL","Label","Label2","Scenario")

# Panel A) Sensitivity analysis of single-factor scenarios
nA = 21 # number of scenarios for panel A
resultsA = droplevels(subset(ExtantPops, ExtantPops$ScenarioNumber %in% c(1:nA)))

(fig5a = ggplot(data=resultsA, aes(x=Year, y=PopsExtant, group=Scenario)) +
    geom_line(aes(color=Scenario), lwd=1.3) +
    geom_line(aes(y=UCL, color=Scenario), lwd=0.8, linetype="dashed") +
    geom_line(aes(y=LCL, color=Scenario), lwd=0.8, linetype="dashed") +
    theme(legend.position=c(0.15, 0.33),
          legend.background = element_rect(fill = "white"),
          axis.text=element_text(size=18, color = "black"),
          axis.title=element_text(size=20),
          axis.line.x = element_line(colour = 'black', size=1.4, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=1.4, linetype='solid')) +
    labs(x="", y="Number of extant populations") +
    xlim(0,85) + ylim(0,475) +
    geom_label_repel(aes(label = Label), max.overlaps = 15,
                     nudge_x = 25, nudge_y = 25, na.rm = TRUE) +
    theme(legend.position = "right",
          legend.text = element_text(size=15),
          legend.title = element_text(size=18),
          legend.margin=margin(t = 0, unit='cm')) +
    guides(col = guide_legend(ncol = 1)) +
    annotate(geom="text", x=80, y=430, label="A",
             color="black", cex=15)
)

# Panel B) Future prediction analysis for multi-factor scenarios
nB = 11 # number of scenarios for panel A
resultsB = droplevels(subset(ExtantPops, ExtantPops$ScenarioNumber %in% c((nA+1):(nA+nB))))
(fig5b = ggplot(data=resultsB, aes(x=Year, y=PopsExtant, group=Scenario)) +
    geom_line(aes(color=Scenario), lwd=1.3) +
    geom_line(aes(y=UCL, color=Scenario), lwd=0.8, linetype="dashed") +
    geom_line(aes(y=LCL, color=Scenario), lwd=0.8, linetype="dashed") +
    theme(legend.position=c(0.15, 0.33),
          legend.background = element_rect(fill = "white"),
          axis.text=element_text(size=18, color = "black"),
          axis.title=element_text(size=20),
          axis.line.x = element_line(colour = 'black', size=1.4, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=1.4, linetype='solid')) +
    labs(x="Years", y="Number of extant populations") +
    xlim(0,85) + ylim(0,475) +
    geom_label_repel(aes(label = Label), max.overlaps = 15,
                     nudge_x = 25, nudge_y = 25, na.rm = TRUE) +
    theme(legend.position = "right",
          legend.text = element_text(size=15),
          legend.title = element_text(size=18)) +
    guides(col = guide_legend(ncol = 1)) +
    annotate(geom="text", x=80, y=430, label="B",
             color="black", cex=15)
)

## Plot the two maps side by side
figure5ab = fig5a + fig5b + plot_layout(ncol=1, nrow=2)
plot(figure5ab)
#KJL Save in version folder
ggsave(paste0(version,"/figure5.png"), width=12, height=12)
ggsave(paste0(version,"figure5_tall.png"), width=10, height=16)

# Open file with Microsoft Paint, and manually adjust panel A legend
#   to be more flush with the graph, like in Panel B





################
### Figure 6 ###
################

### A map of persistence probabilities for each local pop and metapopulation

# (A) local populations first
head(RESULTS80$ScenarioResults[[26]]) # Results for 'Management (less) + medium threats' scenario
prediction = RESULTS80$ScenarioResults[[26]] # Results for medium threats + less management
prediction <- prediction[, !duplicated(colnames(prediction))]
PP = 1-prediction$ProbExtinction ## PP = Persistence Probability

# Categorize persistence probabilities by categories:
# Extremely likely extant (>95%), v. likely extant (80-94.9%),
# likely extant (50-79.9%), and unlikely extant (<50%).
# Then, bind persistence data to data frame
PPCat=matrix(NA,length(PP),1)
for (q in 1:length(PP)){
  if(PP[q] >= 0.95){PPCat[q] = "blue"} else{
    if(PP[q] >= 0.8 & PP[q] < 0.95){PPCat[q] = "green"} else{
      if(PP[q] >= 0.5 & PP[q] < 0.799){PPCat[q] = "yellow"} else{
        PPCat[q] = "orange"}
    }
  }
}
prediction = cbind(prediction, PP, PPCat)

# Plot persistence probability categories by population
i = -90.2; j = 26.3 # Specify the coordinate position (i,j) of the legend
fig6a <- state_map +
  geom_point(data=prediction, aes(x=Longitude, y=Latitude, shape=PPCat),
             colour="black", fill=prediction$PPCat, shape=as.numeric(factor(prediction$PPCat))+20, size=2) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  annotate(x=state_coords$Longitude, y=state_coords$Latitude,
           label=state_coords$State, colour="black", geom="text",size=6)	+
  geom_segment(arrow=arrow(length=unit(5,"mm")),
               aes(x=a, xend=a, y=b, yend=b+l), lwd=2, colour="black") +
  annotate(x=a, y=b-0.4, label="N", colour="black", geom="text", size=10) +
  scale_bar(lon=-91.2, lat=25.5, distance_lon=200,
            distance_lat=20, distance_legend=-25, dist_unit="km", orientation=FALSE) +
  geom_point(aes(x=i-0.8, y=j+1.5), colour = "black", fill="blue", shape=21, size=5) +
  geom_point(aes(x=i-0.8, y=j+1.0), colour = "black",  fill="green", shape=22, size=5) +
  geom_point(aes(x=i-0.8, y=j+0.5), colour = "black", fill="yellow", shape=23, size=5) +
  geom_point(aes(x=i-0.8, y=j+0), colour = "black", fill="orange", shape=24, size=5) +
  annotate(x=i, y=j+1.5, label="Extremely Likely to Persist", geom="text", size=5, hjust=0) +
  annotate(x=i, y=j+1, label="Very Likely to Persist", geom="text", size=5, hjust=0) +
  annotate(x=i, y=j+0.5, label="More Likely Than Not to Persist", geom="text", size=5, hjust=0) +
  annotate(x=i, y=j+0, label="Unlikely to Persist", geom="text", size=5, hjust=0)
plot(fig6a) # blue = PP>0.95, green = 0.95>PP>=0.80, yellow = 0.80>PP>=0.5, orange = 0.50>PP

table(prediction$PPCat, prediction$GeneticUnits)
table(prediction$GeneticUnits, prediction$PPCat)


# Figure 6b
# Adapt the objects made immediately above to map
#   persistence probability for each metapopulation
# Create a matrix to subset and save results for each metapopulation
# Also save metapopulation coordinates as an average of all coordinates
#   for each constituent local population
MetaNumbers = dimnames(table(prediction$LandscapePopID))[[1]]
MetaRows = matrix(NA, length(MetaNumbers), 1)
LatitudeMeta = matrix(NA, length(MetaNumbers), 1)
LongitudeMeta = matrix(NA, length(MetaNumbers), 1)
for (w in 1:length(MetaNumbers)){
  meta = subset(prediction, prediction$LandscapePopID == MetaNumbers[w])
  meta = meta[order(meta$PP, decreasing=TRUE),]
  MetaRows[w,1] = row.names(meta[1,])
  LatitudeMeta[w,1] = mean(meta$Latitude)
  LongitudeMeta[w,1] = mean(meta$Longitude)
}
MetaRows
# For each metapopulation in the 'prediction' object, this identifies the row number
# of the local population with the highest persistence probability


# Subset the 'prediction' object to isolate the best representative local population
# for each metapopulation
MetaPops = prediction[MetaRows,]
MetaPops = cbind(MetaPops, LatitudeMeta, LongitudeMeta)

# Plot persistence probability categories by metapopulation
fig6b <- state_map +
  geom_point(data=MetaPops, aes(x=LongitudeMeta, y=LatitudeMeta, shape=PPCat),
             colour="black", fill=MetaPops$PPCat, shape=as.numeric(factor(MetaPops$PPCat))+20, size=2) +
  annotate(x=state_coords$Longitude, y=state_coords$Latitude,
           label=state_coords$State, colour="black", geom="text",size=6)
plot(fig6b) # blue = PP>0.95, green = 0.95>PP>=0.80, yellow = 0.80>PP>=0.5, orange = 0.50>PP


## Plot the two maps side by side
figure6 = fig6a + fig6b + plot_layout(byrow=TRUE)
plot(figure6)
#KJL Save in version folder
ggsave(paste0(version,"/figure6.png"), width=14, height=6.5)
