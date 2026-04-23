# Title: Distribution drivers: Common Europen viper across southern Germany 
# Author: Sajad Noori
# Date: 05.12.2025

library(sf)
library(sp)
library(pdp)
library(mgcv)
library(terra)
library(dismo)
library(raster)
library(tidyverse)
library(randomForest)
library(rnaturalearth)



# Select the given taxon in the species list
project <- "V_berus"


# set working directory
wdir <- paste0("C:/.../", project)
setwd(wdir)

################################################################################
# species occurrences
################################################################################

# Importing the occurrence dataset
occ <- read.table(paste0(wdir,"/sp_occ.csv"), header = TRUE, sep = ",") 
head(occ)

# Convert occurrences to spatial points
spatial_points <- SpatialPoints(coords = occ[c("x", "y")], 
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
occ_buffered <- raster::buffer(spatial_points, width = 200000) # 200 km



# Define the study area
DE <- ne_states("Germany" , returnclass = "sf") %>%
  filter(name %in% c( "Baden-Württemberg", "Bayern")) %>% 
  st_transform(4326)
  
################################################################################
# Resulting SDM 
################################################################################
current <- raster(paste0(wdir, "/SDM_results/", project, "_current.tif"))
m_current <- mask(crop(current, DE),DE)

sdm_raster <- m_current 
plot(sdm_raster)


################################################################################
# Environmental/anthropogenic drivers (raster files)
################################################################################
env_dir <- paste0("C:/Projects_SAJAD_NOORI_2025/Climate_data/Drivers")
drivers <- stack(list.files(paste0(env_dir), ".tif", full.names=TRUE))
drivers <- mask(crop(drivers, DE), DE)
ref <- rast(drivers[[1]])

plot(drivers[[1]])

# Wetland and Tree cover 
dist_wetland <- raster(paste0("C:/.../rasters/dist_to_wetland_de.tif"))
dist_wetland <- mask(crop(dist_wetland, DE), DE)
Tree <- raster(paste0("C:/.../rasters/Tcover_S_Germany.tif"))
Tree <- resample(Tree, drivers[[1]])

aspect <- terrain(drivers[[1]], v = "aspect", unit = "degrees")

my_drivers <- stack(drivers, dist_wetland, Tree, aspect)
names(my_drivers) <- c("Elevation", "Human_Footprint","Human_modification", 
                    "Precipitation", "Temperature", "Dist_wetland", "tree_cover", "aspect" )
################################################################################
# Extract the values of the drivers for species presence 
################################################################################


# Convert the area higher habitat suitability (threshold) to occurrence dataset
sp_df <- as.data.frame(sdm_raster, xy = TRUE, cells = FALSE, na.rm = TRUE) 

# Converting the occurrences dataset to spatial points
point <- SpatialPoints(sp_df[c("x", "y")], 
                       proj4string = CRS("+proj=longlat +datum=WGS84"))

# Extract the environmental values for these points 
ex_val <- as.data.frame(raster::extract (my_drivers, point))
ex_clim <- cbind(sp_df, ex_val)

set.seed(123)

# Define number of samples per group
n_subsample <- 10000

df_clean <- ex_clim %>%
  sample_n(size = min(n(), n_subsample)) %>%  # take up to 1000, or all if fewer
  drop_na() %>% 
  ungroup() %>% 
  select(- Human_Footprint)

head(df_clean)


################################################################################
# Fitting models to see the effect of drivers on Habitat suitability
################################################################################
dir <- "C:/.../"

# install.packages(c("sjPlot", "sjlabelled", "spdep", "adespatial", "MuMIn", "lme4", "glmmTMB"))
# packages for plots
library (sjPlot)
library (sjlabelled)

# Packages for spatial autocorrelation
library (spdep)
library (adespatial)

# packages for dredgeing models and averaging
library (MuMIn)

# Packages with extra modelling functions
library (lme4)
library (glmmTMB)


# Apply standard transformation (Smithson & Verkuilen 2006)
n <- nrow(df_clean)
df_clean$habitat_beta <- (df_clean$habitat.suitability * (n - 1) + 0.5) / n

# Scale predictors (essential for effect interpretation)
vars <- c(
  "Temperature", "Precipitation", "Elevation",
  "Human_modification", "Dist_wetland", "tree_cover",
   "aspect"
)

df_clean[vars] <- scale(df_clean[vars])

# Random effect: Quad must be a factor (spatial blocking)
df_clean$Quad <- interaction(
  floor(df_clean$x / 10),
  floor(df_clean$y / 10)
)


head(df_clean)

# Building the global models

My_model <- glmmTMB(
  habitat_beta ~
    Temperature +
    Precipitation +
    Elevation +
    Human_modification +
    Dist_wetland +
    tree_cover +
    slope +
    aspect +
    (1 | Quad),
  family = beta_family(link = "logit"),
  data = df_clean
)

summary(My_model)
capture.output(summary(My_model), file = paste0(dir, "results/gam_model_summary.txt"))

# Dredging and averaging
options(na.action = "na.fail")
Dredged_bee_S <- dredge(My_model, rank = "AIC")


# ModelAvg_bee_S <- model.avg(Dredged_bee_S,subset=delta<2, fit = TRUE)
best_model <- get.models(Dredged_bee_S, subset = delta == 0)[[1]]
summary(best_model)
capture.output(summary(best_model), file = paste0(dir, "results/gam_best_model_summary.txt"))
options(na.action = "na.omit")

install.packages("broom.mixed")
library(broom.mixed)

coef_df <- tidy(
  best_model,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
) %>%
  filter(term != "(Intercept)") %>%   # drop intercept if desired
  mutate(
    Coefficient = estimate,
    CI_low  = conf.low,
    CI_high = conf.high,
    
    col = ifelse(Coefficient > 0, "Positive", "Negative"),
    
    Parameter_label = case_when(
      term == "aspect" ~ "Aspect",
      term == "Dist_wetland" ~ "Distance to wetland",
      term == "Human_modification" ~ "Human modification",
      term == "Precipitation" ~ "Precipitation",
      term == "Temperature" ~ "Temperature",
      term == "tree_cover" ~ "Tree cover",
      TRUE ~ term
    )
  )

coef_df <- coef_df %>%
  mutate(
    Parameter_label = factor(
      Parameter_label,
      levels = Parameter_label[order(Coefficient)]
    )
  )


# plot it
p <- ggplot(coef_df,
            aes(x = Coefficient, y = Parameter_label)) +
  
  # zero reference line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  
  # Error bars
  geom_segment(
    aes(
      x = CI_low,
      xend = CI_high,
      y = Parameter_label,
      yend = Parameter_label,
      color = col
    ),
    size = 3,
    lineend = "round"
  ) +
  
  # coefficient points
  geom_point(
    aes(color = col),
    size = 5
  ) +
  
  # coefficient values
  geom_text(
    aes(label = round(Coefficient, 2)),
    nudge_x = 0,
    nudge_y = 0.3,
    size = 5, 
    fontface = "bold"
  ) +
  
  # colors
  scale_color_manual(
    values = c(
      "Negative" = "darkred",
      "Positive"   = "steelblue"
    )
  ) +
  
  labs(
    x = "Standardized effect size",
    y = NULL,
    # title = "Standardized effects of environmental predictors"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "noun",
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = "bold"),
    plot.title = element_text(face = "bold")
  )


p
ggsave(paste0(wdir, "/SDM_results/Drivers.png"), plot = p, width = 10, height = 5, dpi = 300)



################################################################################
# Extinction rate vs. distance to the wetlands
################################################################################

# Importing the occurrence dataset
ex <- read.table("C:/.../csv/Extinction_dataset_E1.csv"
                 , header = TRUE, sep = ",") 
head(ex)
dim(ex)

# Convert occurrences to spatial points
points <- SpatialPoints(coords = ex[c("x", "y")], 
                                proj4string = CRS("+proj=longlat +datum=WGS84"))

# Extract the distance to the wetlands values for these points 
val <- as.data.frame(raster::extract (dist_wetland, points))
colnames(val) <- "dist"
ex_clim <- cbind(ex, val) 
head(ex_clim)


set.seed(123)

# Define number of samples per group
n_subsample <- 10000

df_clean <- ex_clim %>%
  sample_n(size = min(n(), n_subsample)) %>%  # take up to 1000, or all if fewer
  drop_na() %>% 
  ungroup()

head(df_clean)

library(ggsignif)

stats_df <- df_clean %>%
  group_by(status) %>%   # just to get medians per group
  summarise(
    median_dist = median(dist, na.rm = TRUE),
    .groups = "drop"
  )

# Compute delta + test
diff_median <- diff(stats_df$median_dist[stats_df$status == c("extinct", "present")])

p_value <- wilcox.test(dist ~ status, data = df_clean)$p.value

label_df <- data.frame(
  x = 1.5,
  y = max(df_clean$dist, na.rm = TRUE) * 0.95,
  label = paste0(
    "Δ median = ", round(diff_median, 1), " m\n",
    "Wilcoxon p = ", signif(p_value, 3)
  )
)

# Plot it
ggplot(df_clean, aes(x = status, y = dist, fill = status)) +
  
  # Violin in the background
  geom_violin(
    trim = FALSE,
    alpha = 0.5,
    width = 0.9,
    color = NA
  ) +
  
  # Boxplot on top
  geom_boxplot(
    width = 0.2,
    alpha = 0.6,
    outlier.size = 2,
    outlier.alpha = 0.5
  ) +
  
  geom_signif(
    comparisons = list(c("extinct", "present")),
    test = "wilcox.test",
    map_signif_level = TRUE,
    y_position = max(df_clean$dist, na.rm = TRUE) * 1.05,
    textsize = 5
  ) +
  
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 4.5,
    fontface = "bold"
  ) +
 
  
  scale_fill_manual(
    values = c("extinct" = "palegoldenrod", "present" = "darkred")
  ) +
  scale_x_discrete(
    labels = c(
      "extinct" = "Extinct",
      "present" = "Extant"
    )
  )+
  
  labs(
    x = NULL,
    y = "Distance to wetlands (m)"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title.y  = element_text(size = 14, face = "bold")
  )

ggsave(paste0(wdir, "/SDM_results/Distance_to_wetlands.png"), width = 8, height = 6, dpi = 300)
