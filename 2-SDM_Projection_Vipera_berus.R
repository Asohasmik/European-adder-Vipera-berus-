
# Title: Projecting the distribution models of common adder (Vipera berus)
# Author: Sajad Noori
# Date: 16.10.2025

library(sp)
library(terra)
library(dismo)
library(raster)
library(biomod2)
library(foreach)
library(eeptools)
library(tidyverse)
library(doParallel)


# Select the given taxon in the species list
taxon <- c("Vipera berus")
project <- "V_berus"


# set working directory
wdir <- paste0("C:/...", project)
setwd(wdir)



################################################################################
# Environmental variables for species
################################################################################
# Loading the environmental variables for current and futuer 
# List of the selected environmental variables for species
clim_list <- list.files(paste(wdir, "/climate", sep=""), ".asc", full.names=TRUE) # Environmental variables
clim <- stack(clim_list)


load(paste0(wdir, "/SDM_results/myBiomodData.RData"))
load(paste0(wdir, "/SDM_results/myBiomodModelOut.RData"))
load(paste0(wdir, "/SDM_results/myBiomodEM.RData"))


# selected models
test <- read.table(paste0(wdir, "/SDM_results/", project, "_test_single_SDMs.csv"), header = TRUE, sep = ",") 
head(test)

all.models <- get_built_models(myBiomodModelOut)
TSS <- test[test$metric.eval == "AUCroc",]
TSS_clean <- TSS[!is.na(TSS$validation),]

# Select the best models
selected.models <- c()
for (s in 1:length(TSS_clean$validation)){
  if (TSS_clean$validation [s] > 0.7) {
    selected.models <- c(selected.models,all.models[s])
  }
}


###############################################################################
# PROJECTING ENVIRONMNETAL VARIABLES 
###############################################################################

senarios <- c("ssp126", "ssp585")
projects <- c("current","2041-2070", "2071-2100")  


options(scipen = 10)
# Name of the selected variables
clim_names = names(clim)

for(p in 1:length(projects)){
  for(s in 1:length(senarios)){


if (p == 1){
  list <- list.files(paste0(wdir, "/proj/", projects [p], "/"), 
                     ".tif", full.names=TRUE) 
  st <- stack(list)
  vars <- subset(st, clim_names)
  names(vars) <- trimws(clim_names)
}else{
  list <- list.files(paste0(wdir, "/proj/", projects [p], "/", senarios [s], "/"), 
                     ".tif", full.names=TRUE) 
  st <- stack(list)
  vars <- subset(st, clim_names)
  names(vars) <- trimws(clim_names)
  }

   # plot(vars[[1]])
###############################################################################
# Projecting the distribution of the species 
###############################################################################
  

  proj <- projects [p]
  new_env <- vars
  crs(new_env) <- crs(clim)
  
  myBiomodProj <- BIOMOD_Projection(
    bm.mod = myBiomodModelOut,
    new.env = new_env,
    proj.name = proj,
    xy.new.env = NULL,
    models.chosen = selected.models,
    binary.meth = NULL,
    filtered.meth = NULL,
    compress = 'xz',
    clamping.mask = TRUE)  
  
  
  ###############################################################################
  # Ensemble the resulting forecasting 
  ###############################################################################
  
  myBiomodEF <- BIOMOD_EnsembleForecasting(
    bm.em = myBiomodEM,  
    bm.proj = myBiomodProj,
    models.chosen = "all",
    total.consensus = TRUE,   
    binary.meth = NULL,
    filtered.meth = NULL,
    build.clamping.mask = TRUE,
    parallel = TRUE)    
  

  
  spplot<- terra::unwrap(myBiomodEF@proj.out@val) [[2]]
  #plot(spplot)
  
  
  ###############################################################################
  # Re-scaling SDM results
  ###############################################################################
  
  df <- as.data.frame(spplot, xy=TRUE)
  prediction <- df [,3]
  element <- cutoff(prediction, .1) #return minimum number of elements to account 10 percent of total
  cutoff <- prediction [element]
  
  Ensemble <- spplot
  Ensemble[Ensemble < cutoff] <- NA
  cell_values <-values(Ensemble)
  r.min <- min(cell_values, na.rm = TRUE)
  r.max <- max(cell_values, na.rm = TRUE)
  
  r.scale <- ((Ensemble - r.min) / (r.max - r.min)) 
  names(r.scale) <- "habitat suitability"
   plot(r.scale)
  
  if (p == 1){
    dir = paste0(wdir,"/SDM_results/", project, "_",projects [p])
  }else{
    dir = dir = paste0(wdir,"/SDM_results/", project, "_",projects [p], "_", senarios [s])
  }
  
  output_path <- file.path(paste0(dir, ".tif"))
  terra::writeRaster(r.scale, filename = output_path, overwrite = TRUE)
  # plot(r.scale)
  ###############################################################################
  # Generating the MESS layer
  ###############################################################################
  
  # Multivariate environmental similarity surfaces: Elith et al., (2010)
  bg <- cbind(myBiomodData@data.species, myBiomodData@data.env.var)
  bg <- bg[!is.na(bg[,1]),]
  mess <- mess(vars , bg[,2:ncol(bg)], full=F)
  
  mess[values(mess>=0)] <- NA
  mess[values(mess<0)] <- 1
  # plot(mess)
  writeRaster(mess, paste0(dir, "_MESS.tif"), overwrite=TRUE)  
  
}}







