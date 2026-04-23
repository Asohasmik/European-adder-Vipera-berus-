# Title: Current vs. future between geographical distribution and niche shifts for common adder (Vipera berus)
# Author: Sajad Noori
# Date: 16.10.2025

library(sf)
library(sp)
library(terra)
library(dismo)
library(raster)
library(ggpubr)
library(ggExtra)
library(ecospat)
library(rstatix)
library(tidyverse)
library(rnaturalearth)

# Select the given taxon in the species list
taxon <- c("Vipera berus")
project <- "V_berus"

# set working directory
wdir <- paste0("C:/.../", project)
setwd(wdir)

# Define the study area
country <- ne_states("Germany" , returnclass = "sf") %>%
  st_transform(4326)

# Select southern states
area <- country %>%
  filter(name %in% c( "Baden-Württemberg", "Bayern"))
area_sp <- as(area, "Spatial")
################################################################################
# Resulting SDM 
################################################################################
current <- raster(paste0(wdir, "/SDM_results/V_berus_current.tif"))
ssp126_70 <- raster(paste0(wdir, "/SDM_results/V_berus_2041-2070_ssp126.tif"))
ssp585_70 <- raster(paste0(wdir, "/SDM_results/V_berus_2041-2070_ssp585.tif"))
ssp126_100 <- raster(paste0(wdir, "/SDM_results/V_berus_2071-2100_ssp126.tif"))
ssp585_100 <- raster(paste0(wdir, "/SDM_results/V_berus_2071-2100_ssp585.tif"))



# Define your raster files
r_list <- list(
  current ,
  ssp126_70 ,ssp585_70 ,
  ssp126_100, ssp585_100
)

masked_list <- lapply(r_list, function(r) {
  mask(r, area_sp)
})


project_to_equal_area <- function(r) {
  projectRaster(r, crs = CRS("+proj=moll"), method = "bilinear")
}

rasters <- lapply(masked_list, project_to_equal_area)
names(rasters) <- c("current", "ssp126_70", "ssp585_70", "ssp126_100", "ssp585_100")
r_names <- c("current", "ssp126_70", "ssp585_70", "ssp126_100", "ssp585_100")
################################################################################
# Changing in the habitat suitability
################################################################################

th <- c(0.25, 0.5) # thresholds

# Function to calculate area above a threshold

calculate_area_km2 <- function(r, threshold) {
  suitable <- r > threshold
  pixel_area_km2 <- prod(res(r)) / 1e6  # pixel area in km²
  sum(getValues(suitable), na.rm = TRUE) * pixel_area_km2
}


# Create results table
results <- expand.grid(
  scenario = names(rasters),
  threshold = th
) %>%
  arrange(threshold) %>%
  mutate(area_km2 = NA)

# Fill in results
for (i in 1:nrow(results)) {
  r <- rasters[[results$scenario[i]]]
  t <- results$threshold[i]
  results$area_km2[i] <- calculate_area_km2(r, t)
}

# Pivot to wide format to compute change
results_wide <- results %>%
  pivot_wider(names_from = scenario, values_from = area_km2) %>%
  mutate(
    # 2041-70
    change_ssp126_70 = ssp126_70 - current,
    percent_change_ssp126_70 = 100 * (ssp126_70 - current) / current,
    change_ssp585_70 = ssp585_70 - current,
    percent_change_ssp585_70 = 100 * (ssp585_70 - current) / current,
    
    # 2071-100
    change_ssp126_100 = ssp126_100 - current,
    percent_change_ssp126_100 = 100 * (ssp126_100 - current) / current,
    change_ssp585_70 = ssp585_70 - current,
    percent_change_ssp585_100 = 100 * (ssp585_100 - current) / current
  )

# View and save
print(results_wide)
write.csv(results_wide, "SDM_results/suitability_area_changes.csv", row.names = FALSE)

################################################################################
# Habitat suitability change 
################################################################################

base <- rasters[[1]] # Current distribution as baseline

for (r in 2:length(rasters)) {
for (k in 1:length(th)) {

  # Current threshold value
  thr <- th[k]
  
  change_map <- overlay(
    rasters[[r]], 
    base, 
    fun = function(future, curr) {
      ifelse(curr >= thr & future < thr, 1,
             ifelse(curr >= thr & future >= thr, 2,
                    ifelse(curr < thr & future >= thr, 3,
                           0)))
    }
  )
  

  # Define colors including the 0 category (unsuitable)
  colors <- c("gray", "orangered", "gold","forestgreen")
  
  png(filename = paste0(wdir,"/SDM_results/current_vs_", r_names[r],  "_", thr, ".png"), 
      width = 8, height = 6, units = "in", res = 300)
  
  # Plot with legend (use 'levels' to associate categories)
  plot(change_map, col = colors,
       legend = TRUE, 
       main = paste0("Change Map of ", r_names[r], " threshold = ", thr),
       legend.args=list(text='Habitat change categories', side=4, font=2, line=2.5, cex=0.8))
  
  legend("topright", legend=c("Unsuitable","Loss","Stable","Gain"), fill=colors, bty="n")
  dev.off()
  # Save raster
  writeRaster(change_map, 
              filename = paste0(wdir, "/SDM_results/current_vs_", r_names[r],  "_", thr, ".tif"), 
              overwrite = TRUE)
  

}}


################################################################################
# Altitudinal shifts in the habitat suitability
################################################################################
library(ggsignif)

elevation <- raster("C:/...//rasters/elev.tif")
crs(elevation) <- CRS("+proj=moll")

# Initialize empty data frame
elev_df <- data.frame()



thr <- 0.25
for (r in 2:length(r_names)) {

  shift <- raster(paste0(wdir, "/SDM_results/current_vs_", r_names[r],  "_", thr, ".tif"))
  shift_ll <- projectRaster(shift, crs = CRS("+proj=longlat +datum=WGS84"), method="ngb")
  
  # Get logical rasters for Gain and Loss
  gain_r <- shift_ll >= 2.5
  loss_r <- shift_ll >= 0.5 & shift_ll <= 1.5
  
  # Convert raster cells to points (coordinates + value)
  gain_points <- rasterToPoints(gain_r, fun = function(x) x == 1)
  loss_points <- rasterToPoints(loss_r, fun = function(x) x == 1)
  
  
  
  gain_sp <- SpatialPointsDataFrame(coords = gain_points[,1:2],
                                    data = data.frame(value = gain_points[,3]),
                                    proj4string = crs(elevation))
  
  loss_sp <- SpatialPointsDataFrame(coords = loss_points[,1:2],
                                    data = data.frame(value = loss_points[,3]),
                                    proj4string = crs(elevation))

  
  gain_elev <- raster::extract(elevation, gain_sp)
  loss_elev <- raster::extract(elevation, loss_sp)
  
  gain_elev <- as.numeric(gain_elev[!is.na(gain_elev)])
  loss_elev <- as.numeric(loss_elev[!is.na(loss_elev)])
  
  # Add to data frame
  if(length(gain_elev) > 0) {
    elev_df <- rbind(elev_df,
                     data.frame(
                       scenario = r_names[r],
                       category = "Gain",
                       elevation = gain_elev
                     ))
  }
  
  if(length(loss_elev) > 0) {
    elev_df <- rbind(elev_df,
                     data.frame(
                       scenario = r_names[r],
                       category = "Loss",
                       elevation = loss_elev
                     ))
  }
}

head(elev_df)

elev_df$scenario <- factor(
  elev_df$scenario,
  levels = c("ssp126_70",  "ssp585_70",
             "ssp126_100", "ssp585_100")
)


# For each scenario
stats_df <- elev_df %>%
  group_by(scenario) %>%
  summarise(
    median_gain  = median(elevation[category == "Gain"], na.rm = TRUE),
    median_loss  = median(elevation[category == "Loss"], na.rm = TRUE),
    diff_median  = median_gain - median_loss,
    p_value      = wilcox.test(
      elevation[category == "Gain"],
      elevation[category == "Loss"]
    )$p.value
  )

stats_df <- stats_df %>%
  mutate(
    label = paste0("Δmedian = ", round(diff_median,1), 
                   "\nWilcox p = ", signif(p_value,3))
  )




ggplot(elev_df, aes(x = category, y = elevation, fill = category)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~scenario) +
  scale_fill_manual(values = c("Loss"= "orangered", "Gain"="forestgreen")) +
  geom_signif(
    comparisons = list(c("Loss", "Gain")),
    map_signif_level = TRUE,
    vjust = 0.5,
    textsize = 4 
  )+
  geom_text(
    data = stats_df,
    aes(x = 1.5, y = max(elev_df$elevation, na.rm = TRUE)*0.95, label = label),
    inherit.aes = FALSE,
    hjust = 0.5,
    size = 4
  )+
  labs(
    x = "",
    y = "Elevation (m)"
  )+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),        # x/y tick labels
    axis.title = element_text(size = 14),       # x/y axis titles
    strip.text = element_text(size = 14)       # facet labels

  )

ggsave(paste0(wdir, "/SDM_results/Altitudinal_shifts.png"), width = 10, height = 8, dpi = 300)


################################################################################
# Niche dis-similarity: base line sv. future scenarios
################################################################################

# List of the selected environmental variables for species
clim <- stack(list.files(paste(wdir, "/climate", sep=""), ".asc", full.names=TRUE))# Environmental variables
scenarios <- c("current", "ssp126_70", "ssp585_70", "ssp126_100", "ssp585_100")


output_file <- paste0(wdir,"/SDM_results/climate_values.csv") 
  

process_scenarios <- function(r_list, scenarios, th, clim, output_file) {
  results_df <- data.frame()
  
  for (i in 1:length(scenarios)) {
    for (j in 1:length(th)) {
      sp_df <- as.data.frame(r_list[[i]], xy = TRUE, na.rm = TRUE) %>% 
        filter(habitat.suitability > th[j])
      
      if(nrow(sp_df) == 0) next  # Skip if no values passed the threshold
      
      point <- SpatialPoints(sp_df[c("x", "y")], 
                             proj4string = CRS("+proj=longlat +datum=WGS84"))
      
      ex_val <- as.data.frame(raster::extract(clim, point))
      ex_clim <- cbind(sp_df, ex_val)
      
      df_ex_clim <- ex_clim %>%
        mutate(
          scenario = scenarios[i],
          threshold = th[j],
        ) %>%
        drop_na()
      
      results_df <- rbind(results_df, df_ex_clim)
    }
  }
  
  #write.csv(results_df, output_file, row.names = FALSE)
  return(results_df)
}
results <- process_scenarios(r_list, scenarios, th, clim, "my_output.csv")


# Define number of samples per group
n_subsample <- 1000

subsampled_df <- results %>%
  group_by(scenario, threshold) %>%
  sample_n(size = min(n(), n_subsample)) %>%  # take up to 1000, or all if fewer
  ungroup()
  
head(subsampled_df)

pca_df <- subsampled_df %>%
  select(starts_with("bio")) %>%
  mutate(across(everything(), as.numeric))
  
# Perform PCA
pca <- prcomp(pca_df, scale. = TRUE)

# Attach PCA scores to original data
pca_df <- as.data.frame(pca$x) %>%
  bind_cols(subsampled_df %>% select(scenario, threshold))

head(pca_df)

# Calculate group means (centroids)
centroids_df <- pca_df %>%
  group_by(scenario, threshold) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")


#  Loop over thresholds and make plots
# set the base line 
baseline <- "current"

# Get the PCA contribution result
pca_var_perc <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)


cols <- c(
  "current" = "forestgreen",
  "ssp126_70" = "gold",
  "ssp585_70" = "orangered",
  "ssp126_100" = "#C71585",
  "ssp585_100" = "#8B008B"
)



for (p in 1:length(th)){
  t <- as.character(th[p])
  # Example subset
  sub <- pca_df %>% filter(threshold == t)
  
  scenarios <- setdiff(unique(sub$scenario), baseline)
  scenario_order <- c("current", "ssp126_70", "ssp585_70", "ssp126_100", "ssp585_100")
  sub$scenario <- factor(sub$scenario, levels = scenario_order)
 
  
  # Calculate p-values for PC1 and PC2 against baseline
  pc1_p <- sapply(scenarios, function(scen) {
    t.test(sub$PC1[sub$scenario == baseline], sub$PC1[sub$scenario == scen])$p.value
  })
  
  pc2_p <- sapply(scenarios, function(scen) {
    t.test(sub$PC2[sub$scenario == baseline], sub$PC2[sub$scenario == scen])$p.value
  })
  
  # Create text for annotation
  annot_text <- paste0(names(pc1_p), ": PC1 p=", signif(pc1_p, 3), ", PC2 p=", signif(pc2_p, 3))
  annot_text <- paste(annot_text, collapse = "\n")
  
  # PCA plot
  p <- ggplot(sub, aes(PC1, PC2, color = scenario, group = scenario)) +
    geom_point(alpha = 0.3) +
    stat_ellipse(aes(group = scenario), type = "norm", level = 0.95)+
    scale_color_manual(values = cols) +
    theme(legend.position = c(0.12, 0.9),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.key.size = unit(0.5, "cm"),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"))+
    # annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, 
    #          label = annot_text, size = 3.5, color = "black")+
    # annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
    #          label = paste0("Threshold: ", t), size = 5, fontface = "bold") +
    labs(x = paste0("PC1 (", pca_var_perc[1], "%)"), 
         y = paste0("PC2 (", pca_var_perc[2], "%)"))
  
  # Add marginal boxplots
  p_margin <- ggMarginal(p, type = "boxplot", groupColour = TRUE, groupFill = TRUE)
  p_margin
  
  ggsave(paste0(wdir, "/SDM_results/pca_", t , "_plot.png"), plot = p_margin, width = 8, height = 8, dpi = 300)
}


?wilcox.test

# Calculate Schoener's D

niche_test <- data.frame()

for (thresh in unique(pca_df$threshold)) {
  #thresh = 0.25

  sub <- pca_df %>% filter(threshold == thresh)
  
  baseline_data <- sub %>% filter(scenario == "current") %>% select(PC1, PC2)
  scenarios <- setdiff(unique(sub$scenario), baseline)
  
  for (scen in scenarios) {
    #scen = "ssp126_100"
    scenario_data <- sub %>% filter(scenario == scen) %>% select(PC1, PC2)
    
    # Prepare niche objects with ecospat::ecospat.grid.clim.dyn or similar
    # Here I assume you create niche objects from the PCA scores for each scenario
    niche1 <- ecospat.grid.clim.dyn(glob = pca_df %>% select(PC1, PC2),
                                    glob1 = baseline_data,
                                    sp = baseline_data,
                                    R = 100)
    
    niche2 <- ecospat.grid.clim.dyn(glob = pca_df %>% select(PC1, PC2),
                                    glob1 = scenario_data,
                                    sp = scenario_data,
                                    R = 100)
    
    # Compute niche overlap metrics (D and I)
    niche_overlap <- ecospat.niche.overlap(niche1, niche2, cor = TRUE)
    
    # Conduct t-test on PC1 (or on PC2 or combined)
    w_test <- wilcox.test(baseline_data$PC1, scenario_data$PC1)
    
    # Store results
    df <- data.frame(
      overlap = paste0(baseline, " vs ", scen),
      w = w_test$statistic,
      #df = w_test$parameter,
      p_value = case_when(
        w_test$p.value < 0.001 ~ 0.001,
        w_test$p.value < 0.01 ~ 0.01,
        w_test$p.value < 0.05 ~ 0.05,
        TRUE ~ w_test$p.value 
      ),
      D = round(niche_overlap$D, 2),
      I = round(niche_overlap$I, 2),
      threshold = thresh
    )
    
    niche_test <- rbind(niche_test, df)
  }
}
# Save the values for further analysis
write.csv(niche_test, paste0(wdir, "/SDM_results/Niche_qualification_test.csv") , row.names = FALSE)
citation ("stats")
################################################################################
# hypervolume
################################################################################
# Build hypervolumes
library(hypervolume)
citation("hypervolume")
# A hypervolume represents the multidimensional region that the morphs occupy in environmental space.

# Compute hypervolumes for all populations
# estimate the probabilistic shape of the population’s niche 
## in that space, using a Gaussian kernel density estimator (KDE)

scenarios <- c("current", "ssp126_100", "ssp126_70", "ssp585_100", "ssp585_70" )
hv_list <- list()

for (s in 1:length(scenarios)) {

  
  sen <- scenarios [s]
  
  scen <- pca_df %>%
    filter(scenario == sen & threshold == 0.25) %>%
    select(PC1, PC2, PC3)
  
  hv_list[[s]] <- hypervolume_gaussian(scen,
                                       name = paste0(sen),
                                       samples.per.point = 100,
                                       chunk.size = 1000)
}


names(hv_list)

# Compute pairwise overlaps between hypervolumes
overlap_results <- list()
counter <- 1

for (i in 1:(length(hv_list) - 1)) {
  for (j in (i + 1):length(hv_list)) {

    hv1 <- hv_list[[i]]
    hv2 <- hv_list[[j]]
    ?hypervolume_overlap_statistics
    hv_set <- hypervolume_set(hv1, hv2, check.memory = FALSE)
    ov_stats <- hypervolume_overlap_statistics(hv_set)
    vol <- get_volume(hv_set)
    
    
    
    overlap_results[[counter]] <- data.frame(
      scen1 = hv1@Name,
      scen2 = hv2@Name,
      
      sorensen = ov_stats[2],
      jaccard = ov_stats[1],
      
      # volume
      scen1_vol = vol [1],
      scen2_vol = vol [2],
      intersection_vol = vol [3])
    
    counter <- counter + 1
  }
}

overlap_df <- dplyr::bind_rows(overlap_results)%>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
  as.tibble() %>% 
  as.data.frame()
overlap_df

write.table(overlap_df, paste0(wdir, "/SDM_results/Hypervolum_test.csv"))

################################################################################
# Unsupervised niche clustering
################################################################################
library(factoextra)
citation("factoextra")
library(ggExtra)

# Importing the occurrence dataset
occ <- read.table(paste0(wdir,"/PA_dataset.csv"), header = TRUE, sep = ",") %>% 
  filter(PA == 1)
head(occ)

df <- occ %>% 
  select(x, y)
dim(df)

gsp <- SpatialPoints(df[c("x", "y")], 
                     proj4string = CRS("+proj=longlat +datum=WGS84"))

ext <- terra::extract(clim, gsp) %>% 
  bind_cols(df) %>% 
  na.omit()
dim(ext)

# Perform PCA on climate variables
pca <- prcomp(ext [,1:8], scale. = TRUE)
pca_var_perc <- summary(pca)

pca_scores <- as.data.frame(pca$x[, 1:3]) %>%  # for example, use the first 3 PCs
  select_if(is.numeric)

pca$rotation
# Try different cluster numbers
fviz_nbclust(pca_scores, kmeans, method = "wss")  # within-cluster sum of squares (Elbow method)
dev.off()

# Run k-means 
k_clust <- kmeans(pca_scores, centers = 5, nstart = 25)

# Add cluster labels to your data
pca_clustered <- cbind(pca_scores, cluster = as.factor(k_clust$cluster))
head(pca_clustered)
dim(pca_clustered)

# Add cluster assignments as a new column in your occurrence dataset
occurrence_with_clusters <- cbind(ext [c("x", "y")], cluster = as.factor(k_clust$cluster))

pallet <- c(
  "1" = "#D2691E",
  "2" = "#2E8B57",
  "3" = "#6B8E23",
  "4" = "#483D8B",
  "5" = "#DAA520"
)

pl <- ggplot() +
  geom_sf(data = area, fill = NA, color = "black", linewidth = 0.6)+ 
  geom_point(data = occurrence_with_clusters, aes(x = x, y = y, color = cluster), alpha = 0.4, size = 3)+ 
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  scale_color_manual(values = pallet)+
  # theme_minimal(base_size = 14) +
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = c(0.08, 0.90),
        # margin = unit(0, "cm"),
        legend.key.size = unit(0.5, "cm"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"))
pl

ggsave(paste0(wdir,"/SDM_results/Unsuperviszed_PCA_clusters.png"), pl, width = 12, height = 8, dpi = 300)

# Calculate group means (centroids)
centroids_df <- pca_clustered %>%
  group_by(cluster) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

# Get the PCA contribution result
pca_var_perc <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

# PCA plot
p <- ggplot(pca_clustered, aes(PC1, PC2, color = cluster, group = cluster)) +
  geom_point(alpha = 0.4, size = 3) +
  stat_ellipse(aes(group = cluster), type = "norm", level = 0.95)+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  scale_color_manual(values = pallet) +
  theme(legend.position = c(0.06, 0.9),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"))+
  # annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, 
  #          label = annot_text, size = 3.5, color = "black")+
  # annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1,
  #          label = paste0("Threshold: ", t), size = 5, fontface = "bold") +
  labs(x = paste0("PC1 (", pca_var_perc[1], "%)"), 
       y = paste0("PC2 (", pca_var_perc[2], "%)"))

# Add marginal boxplots
p_margin <- ggMarginal(p, type = "boxplot", groupColour = TRUE, groupFill = TRUE)
p_margin

ggsave(paste0(wdir, "/SDM_results/PCA_Unsuperviszed_groups.png"), plot = p_margin, width = 8, height = 8, dpi = 300)
