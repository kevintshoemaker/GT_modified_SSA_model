Code README file for:
Inflated predictions from a flawed model influenced the decision to deny federal protection for the gopher tortoise
A response to Folt et al. 2022 Global Ecology & Conservation
Loope, Ackakaya & Shoemaker

Created by: KJ Loope
12/7/23

Scripts
tortoise-pva_response.R
This file contains a modified version of the original tortoise-pva.R script.  It is modified only as described in the Supplementary Methods, and each modification is indicated by a comment containing the modifier's initials (KJL).  At the top of the document, the user can choose between three versions of the model by specifying "this.version".  this.version=1 runs the original model and outputs results into a folder called "original".  Specifying this.version=2 runs the model with a correction for the density dependence threshold typo, and outputs results into a folder called "dd".  Specifying this.version=3 runs the model with a correction for the density dependence typo and implements a 3% annual growth cap on metapopulations, and outputs results into a folder called "3percent".  It requires the same input files as the original code:
demographic-rates.csvfunctions.Rstate-coords.csvtortoise.pngUnit4.kmlUnit5.kmlUnit1.kmlUnit2.kmlUnit3.kml
tortoise-surveys.csv [available by writing to USFWS Jacksonville Office; see main text]
[See https://code.usgs.gov/cooperativeresearchunits/modeling-gopher-tortoises-population-viability]

tortoise-pva_response_ages.R
This file implements the original model with corrected density dependence typo, a 3% cap on annual metapopulation growth, and a matrix model that models all juvenile age classes explicitly.  It requires the same input files as the original code.  Each modification from the original code is indicated by the modifier's initials (KJL).  The code is set up to run in parallel using a computing cluster (set option cluster=T to run with 32 cores; runs in ~11 minutes).  If not using a computing cluster, set cluster = F; this will run using 2 cores on a personal computer (runs in ~3h).   

response_figs.R
This script generates all figures and analyses reported in the text.  It imports replicate-level output data from the above original version of the model (Figures 2 & 3).  This requires the 44mb file Code/original/output/scen-26.rds to create Figure 2 and S2.  A summary of the relevant data required to make Figure 3 is provided in Code/original/original_Nm_summary.csv; this can also be re-created from the full output of all scenarios if the model script is re-run. These output files are available upon request (~1.27gb) if you prefer not to re-run the model script.  To create figures reporting results across model versions (Figs 5 and S3-4), the script imports table3.csv from each of the four model version output folders.  These are included but will be re-created if model scripts are run.  The script also requires the file pop_metapop_key.csv to link populations to their metapopulations to create Figure S1. 

Data Files from Folt et al. 2022:
demographic-rates.csv [from Folt et al. 2022]functions.R [from Folt et al. 2022]state-coords.csv [from Folt et al. 2022]tortoise.png [from Folt et al. 2022]Unit4.kml [from Folt et al. 2022]Unit5.kml [from Folt et al. 2022]Unit1.kml [from Folt et al. 2022]Unit2.kml [from Folt et al. 2022]Unit3.kml [from Folt et al. 2022]
tortoise-surveys.csv [from Folt et al. 2022; available by writing to USFWS Jacksonville Office; see main text]

Data Files from this letter to the editor:
original/original_Nm_summary.csv [summary of outputs from original run used to create Figure 3.]
pop_metapop_key.csv [used to link pops to metapops to create Figure S1]



Note on packages:  the original code relies on early versions of elevatR and maptools that have since been updated to reflect the retirement of rgdal, maptools and sp packages.  The original code does not run with the updated versions of these packages.  To run as is, one must use an archived version of R (=<4.2.3) and install archived versions of these packages; for exact versions, see below.
  
We executed the scripts using the following platforms:

R version 4.1.0 (2021-05-18)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux

Attached packages and versions:
 [1] legendMap_1.0   png_0.1-7       patchwork_1.1.1 ggrepel_0.9.1  
 [5] reshape_0.8.8   mapdata_2.3.1   maps_3.4.1      ggsn_0.5.0     
 [9] dplyr_1.0.8     maptools_1.1-6  ggmap_3.0.2     ggplot2_3.4.2  
[13] clipr_0.8.0     plotrix_3.8-1   plyr_1.8.6      lognorm_0.1.10 
[17] elevatr_0.4.2   raster_3.6-3    sp_1.6-0        installr_0.23.2


The same code also runs on a Intel MacBook Air with the following packages and versions:

R version 4.2.3 (2023-03-15)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.7

 [1] legendMap_1.0   png_0.1-8       patchwork_1.1.2 ggrepel_0.9.3  
 [5] reshape_0.8.9   mapdata_2.3.1   maps_3.4.1      ggsn_0.5.0     
 [9] dplyr_1.1.1     maptools_1.1-6  ggmap_3.0.2     ggplot2_3.4.2  
[13] clipr_0.8.0     plotrix_3.8-2   plyr_1.8.8      lognorm_0.1.10 
[17] elevatr_0.4.2   raster_3.5-15   sp_1.6-0        installr_0.23.4

