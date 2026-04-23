# Title:    Assessing the niche of the Vipera berus across climate scenarios
# Author:   Sajad Noori
# Date:     24.09.2025

library(sf)
library(usdm)
library(ade4) # dudi.pca
library(httr)
library(vegan)
library(terra)
library(dismo)
library(raster)
library(cowplot)
library(ecospat)
library(ggExtra)
library(tidyverse)
library(hypervolume)
library(rnaturalearth)


# Set working directory
wdir <- "C:/..."
setwd(wdir)



################################################################################
# Uploading required data and datasets
################################################################################

# occurrnce dataset
tbl <- read.table(paste0(wdir, "/csv/All_data.csv"), sep = ",", dec = ".", header = TRUE) %>% 
  filter(!is.na(population))

head(tbl)
dim(tbl)
unique(tbl$population)
summary(tbl$year)


# Climate data

clim_dir <- "C:/.../"

clim <- raster::stack(
  list.files(clim_dir, pattern = "\\.tif$", full.names = TRUE)
)



# Load Germany states
germany_states <- ne_states(country = "Germany", returnclass = "sf") %>%
  st_transform(4326)

# Select southern states
S_germany <- germany_states %>%
  filter(name %in% c( "Baden-Württemberg", "Bayern"))
#plot(st_geometry(S_germany))


Tcover <- raster(paste0(wdir, "/rasters/Tcover_S_Germany.tif"))
#plot(Tcover)

################################################################################
# Climate space of southern Germany usign PCA
################################################################################
# Get random points across the study area to extract the climatic variable 
n_points <- 10000
random_pts <- st_sample(S_germany, size = n_points, type = "random")
random_pts_sp <- as(random_pts, "Spatial") # Convert to Spatial

clim_values <- raster::extract(clim , random_pts_sp) # Extract climate values at those points
tcover_val <- raster::extract(Tcover , random_pts_sp)
clim_df <- cbind(st_coordinates(random_pts), clim_values, tcover_val) %>% as.data.frame()

head(clim_df)
dim(clim_df)
# Find the most relavent environmental variables
# PCA
env <- na.omit(as.data.frame(clim_df [, 3:24]))

pca.clim <- dudi.pca(env, center = TRUE,
                     scale = TRUE, scannf = FALSE, nf = 2)
png(paste0(wdir,"/graphs/PCA.png"),
    width = 20, height = 20, units = "cm", res = 300)
contrib <- ecospat.plot.contrib(pca.clim$co , eigen = pca.clim$eig)
# text(pca.clim$co[,1], pca.clim$co[,2], labels = rownames(pca.clim$co), pos = 3, cex = 0.7)
dev.off()

vars <- data.frame(clim_df [, 3:24])  # Convert to a data frame
myvar <- vifstep(vars, th = 5)
myvar


selected_vars <- c("bio_18","bio_16", "bio_14", "bio_15", "bio_9", "bio_8", "bio_6",  "bio_5", "bio_7")

bioclim <- subset(clim, selected_vars) 



################################################################################
# Species morphs across PCA climate space
################################################################################

# Bioclimatic rasters (subset of selected variables)
selected_vars <- c("bio_18","bio_16", "bio_14", "bio_15",
                   "bio_9", "bio_8", "bio_6",  "bio_5", "bio_7")

bioclim <- subset(clim, selected_vars) 


# Extract climate values for the study area
# Get random points across the study area to extract the climatic variable 
n_points <- 10000
random_pts <- st_sample(S_germany, size = n_points, type = "random")
random_pts_sp <- as(random_pts, "Spatial") # Convert to Spatial

clim_values <- raster::extract(bioclim , random_pts_sp) # Extract climate values at those points
clim_df <- cbind(st_coordinates(random_pts), clim_values) %>% 
  as.data.frame() %>% 
  rename(x = X, y = Y) %>% 
  mutate(data = "BG") 

head(clim_df)
dim(clim_df)



################################################################################
# First unsupervised clustering
################################################################################
set.seed(42)
library(factoextra)

df <- tbl %>% 
  select(x, y)

gsp <- SpatialPoints(df[c("x", "y")], 
                     proj4string = CRS("+proj=longlat +datum=WGS84"))

ext <- terra::extract(bioclim, gsp) %>% 
  bind_cols(df) %>% 
  na.omit()


# Perform PCA on climate variables
pca <- prcomp(ext [,1:9], scale. = TRUE)
pca_var_perc <- summary(pca)

pca_scores <- as.data.frame(pca$x[, 1:3]) %>%  # for example, use the first 3 PCs
  select_if(is.numeric)

png(paste0(wdir,"/graphs/fviz_nbclust.png"),
    width = 20, height = 20, units = "cm", res = 300)
# Try different cluster numbers
fviz_nbclust(pca_scores, kmeans, method = "wss")  # within-cluster sum of squares (Elbow method)
dev.off()

pca$rotation
summary(pca)
# Run k-means 
k_clust <- kmeans(pca_scores, centers = 5, nstart = 25)

# Add cluster labels to your data
pca_clustered <- cbind(pca_scores, cluster = as.factor(k_clust$cluster))
head(pca_clustered)


# Add cluster assignments as a new column in your occurrence dataset
occurrence_with_clusters <- cbind(ext, cluster = as.factor(k_clust$cluster))



pallet <- c(
  "1" = "#D2691E",
  "2" = "#2E8B57",
  "3" = "#483D8B",
  "4" =  "#6B8E23",
  "5" = "#DAA520"
)

pl <- ggplot() +
  geom_sf(data = S_germany, fill = NA, color = "black", linewidth = 0.6) +
  geom_point(data = occurrence_with_clusters, aes(x = x, y = y, color = cluster)) +
  scale_color_manual(values = pallet)+
  # theme_minimal(base_size = 14) +
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = c(0.12, 0.90),
        # margin = unit(0, "cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"))
pl

ggsave(paste0(wdir,"/graphs/Unsuperviszed_PCA_clusters.png"), pl, width = 12, height = 8, dpi = 300)
