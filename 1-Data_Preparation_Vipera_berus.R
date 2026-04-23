# Title: Selection of the climate variables for common adder (Vipera berus)
# Author: Sajad Noori
# Date: 16.10.2025


library(sf)
library(usdm)
library(ade4) # dudi.pca
library(raster)
library(spThin)
library(ecospat)
library(corrplot)
library(tidyverse)
library(rnaturalearth)



# set working directory
wdir <- "C:/..."
setwd(wdir)


# Select the given taxon in the species list
taxon <- c("Vipera berus")
project <- "V_berus"

################################################################################
# Cleaning and preparing the occurrence dataset
################################################################################

# Import the occurrence data as csv file
dataset <- read.csv2(paste0(wdir,"/csv/All_data.csv"), header = TRUE, sep = ",") %>%     
  filter(year >= 1980) %>% 
  select(species,x, y) %>% 
  distinct()

head(dataset)
dim(dataset)


# DefinTRUE# Define the cleaning and thinning function
clean_thin_occurrences <- function(data, taxon, wdir, project, thin_km , save) {

  n_original <- nrow(data)
  
  # Select and rename columns
  data <- dataset 
  
  
  # Filter by species
  myspecies <- data %>% filter(species == taxon)
  
  # Remove NAs and duplicates
  myspecies <- myspecies %>%
    filter(!is.na(x) & !is.na(y) & !is.na(species)) %>%
    distinct(species, x, y, .keep_all = TRUE)
  
  myspecies$x <- as.numeric(myspecies$x)
  myspecies$y <- as.numeric(myspecies$y)
  
  n_cleaned <- nrow(myspecies)
  
  # Thin the dataset spatially
  temp_dir <- tempfile(pattern = "thin_output_")
  dir.create(temp_dir) # use a temporary output directory
  suppressMessages(
    suppressWarnings(
      capture.output({
    thinned <- thin(loc.data = myspecies,
                    lat.col = "y",
                    long.col = "x",
                    spec.col = "species",
                    thin.par = thin_km, # in km
                    reps = 50,        
                    write.files = TRUE,
                    out.dir = temp_dir,
                    verbose = FALSE)
    # Read back the best replicate
    thinned_files <- list.files(temp_dir, full.names = TRUE)
    best_thin <- read.csv(thinned_files[1])
    n_thinned <- nrow(best_thin)
    unlink(temp_dir, recursive = TRUE)  # Clean up temp files
      }, file = NULL)
  ))
  
  # Show summary
  cat("=== Record Summary ===\n")
  cat("Original records:      ", n_original, "\n")
  cat("After cleaning (NA + duplicate removal):", n_cleaned, "\n")
  cat("After thinning:        ", n_thinned, "\n")

  # Save result
  if (save) {
    dir.create(file.path(wdir, project), showWarnings = FALSE, recursive = TRUE)
    write.csv(best_thin, file = file.path(wdir, project, "sp_occ.csv"), row.names = FALSE)
  }
  
  return(best_thin)
}

sp_occ <- clean_thin_occurrences(dataset, taxon, wdir, project, 1, TRUE)
head(sp_occ)
dim(sp_occ)

sp_occ <- read.csv2(paste0(wdir,"/", project, "/sp_occ.csv"), header = TRUE, sep = ",")
sp_occ$x <- as.numeric(sp_occ$x)
sp_occ$y <- as.numeric(sp_occ$y)


################################################################################
# Variables selection
################################################################################

# Importing the environmental variables
var_list <- list.files(path =paste0("C:/..."), 
                         pattern=".tif$", all.files=TRUE, full.names=TRUE) # Environmental variables



var_selection <- function(var_list){
  
  variables <- stack(var_list)
  # Creating a directory for variables selection
  dir.create(paste0(wdir,"/",project, "/var_selection/")) 
  
  # Extract the environmental variables for each grid 
  spatial_points <- SpatialPoints(coords = sp_occ[c("x", "y")], 
                                  proj4string = CRS("+proj=longlat +datum=WGS84"))
  sp_Valuaes <- raster::extract(variables, spatial_points)
  #rasterVAlues <- getValues(variables)
  
  M <- cor(sp_Valuaes,use = "pairwise.complete.obs")
  ## for better interpretation we use a threshold
  tmp = M # Copy matrix
  tmp[ tmp > -0.75 & tmp < 0.75 ] = 0
  png(paste0(wdir,"/",project,"/var_selection/corrplot.png"),
      width = 20, height = 20, units = "cm", res = 300)
  corrplot(tmp,
           method = "color",
           addCoef.col="white",
           order = "AOE",
           number.cex=0.75,
           tl.cex = 1,
           addgrid.col = "grey")
  dev.off()
  
  # A VIF greater than 10 is a signal that the model has a collinearity problem
  vars <- data.frame(sp_Valuaes)  # Convert to a data frame
  myvar <- vifstep(vars, th = 5)
  write.table(myvar@results, paste0(wdir,"/",project,"/var_selection/VIFcor.csv"), 
              sep = ",", row.names = FALSE)
  name <- myvar@results$Variables
  # return(name)
  # PCA
  env <- na.omit(as.data.frame(sp_Valuaes))
  
  pca.clim <- dudi.pca(env, center = TRUE,
                       scale = TRUE, scannf = FALSE, nf = 2)
  png(paste0(wdir,"/",project,"/var_selection/PCA.png"),
      width = 20, height = 20, units = "cm", res = 300)
  contrib <- ecospat.plot.contrib(pca.clim$co , eigen = pca.clim$eig)
  dev.off()
  return(myvar)
}

myvariable <- var_selection(var_list)
myvariable

# Here I checked the results of variable selection part and came up with the following variables
var_names <- c("bio02","bio06","bio07","bio08","bio09", "bio11", "bio16","bio19")




# Define the study area
country <- ne_states("Germany" , returnclass = "sf") %>%
  st_transform(4326)

# Select southern states
area <- country %>%
  filter(name %in% c( "Baden-Württemberg", "Bayern"))


# Cropping environmental variables for species ranges
var_cropping <- function(area, rasterList, var_names){
  options(scipen = 10)
  
  
  clim <- stack (rasterList) # Make stack from all variables
  names(clim)
  subset_clim <- subset(clim, var_names) # A subset of all variables for species 
  
  
  dir.create(paste0(wdir,"/", project,"/climate"))
  dir <- paste0(wdir,"/", project,"/climate")
  
  
  # Step 1: Crop each raster to the buffer
  cropped_rasters <- lapply(rasterList, function(f) {
    r <- raster(f)
    crop(r, area)
  })
  
  cropped_rasters <- crop(subset_clim, area)
  
  # Step 4: Write to disk
  file_path <- file.path(dir, paste0(var_names, ".asc"))
  
  writeRaster(cropped_rasters, file_path,
              names(subset_clim), bylayer = TRUE,
              driver = 'HFA', format = "ascii", overwrite = TRUE)
  
  return(subset_clim)
  
}




variables <- var_cropping(area, rasterList, var_names)


