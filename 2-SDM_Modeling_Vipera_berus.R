# Title: Modeling the distribution of common adder (Vipera berus)
# Author: Sajad Noori
# Date: 16.10.2025

library(sp)
library(terra)
library(dismo)
library(raster)
library(biomod2)
library(foreach)
library(ggtext)
library(eeptools)
library(tidyverse)
library(doParallel)


citation("biomod2")
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
var_list <- list.files(paste(wdir, "/climate", sep=""), ".asc", full.names=TRUE) # Environmental variables
var <- stack(var_list)

################################################################################
# Data preparation for SDM
################################################################################
#Creat a folder to save all the SDM results
dir.create(paste0(wdir, "/SDM_results"))

# Importing the occurrence dataset
occ <- read.table(paste0(wdir,"/sp_occ.csv"), header = TRUE, sep = ",") %>% 
  mutate(PA = as.numeric(1))
head(occ)

# Save the SDM results ins folder

# Set up the requierd data for Biomod function
myResp = as.numeric(occ[,'PA'])
myExpl = var
myRespXY = occ [, c('x', 'y')]
myRespName <- project


# Format Data with pseudo-absences : random method
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 3,
                                     PA.nb.absences = 10000,
                                     PA.strategy = 'random')
plot(myBiomodData)
save(myBiomodData, file = paste0(wdir, "/SDM_results/myBiomodData.RData"))
###############################################################################
# Algorithms configurations
##############################################################################

# Cross validation
block_CV <- bm_CrossValidation(
  bm.format = myBiomodData,
  strategy = "kfold",         # or "random"
  k = 5,                   #  the number of partitions
  nb.rep = 2              # number of repetitions
)
dim(block_CV)

Selected_algos <- c("MAXENT", "RF", "GBM") # 
default_options <- bm_ModelingOptions(data.type = 'binary',
                                      models = Selected_algos,
                                      strategy = 'default',    
                                      bm.format = myBiomodData,
                                      calib.lines = block_CV)
# default_options@options

###############################################################################
# Modeling 
###############################################################################
?BIOMOD_Modeling
# Model single models
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    modeling.id = 'AllModels',
                                    models = Selected_algos,
                                    CV.strategy = 'user.defined',
                                    CV.user.table = block_CV,
                                    #nb.rep = 10,
                                    OPT.strategy = 'bigboss',
                                    metric.eval = c('TSS','ROC'),
                                    var.import = 2, # the number of permutations
                                    seed.val = 42)

#stopImplicitCluster()
save(myBiomodModelOut, file = paste0(wdir, "/SDM_results/myBiomodModelOut.RData"))
# load(paste0(wdir, "/SDM_results/myBiomodModelOut.RData"))
###############################################################################
# Modeling summary and select the best models
###############################################################################

# Get the model evaluation results
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
myfilename <- paste0(wdir, "/SDM_results/", project, "_test_single_SDMs.csv")
write.csv(myBiomodModelEval, myfilename)
test <- read.csv(myfilename, h=TRUE, row.names=1)

all.models <- get_built_models(myBiomodModelOut)
ROC <- test[test$metric.eval == "AUCroc",]
ROC_clean <- ROC[!is.na(ROC$validation),]

# Select the best models
selected.models <- c()
for (s in 1:length(ROC_clean$validation)){
  if (ROC_clean$validation [s] > 0.7) {
    selected.models <- c(selected.models,all.models[s])
  }
}

# number of selected models
length(selected.models)

###############################################################################
# Modeling variable importance
###############################################################################

# Get a resport on important variables 
myVarImport <- get_variables_importance(myBiomodModelOut)

VarNames <- as.data.frame(names(var))
myVarImport <- cbind(VarNames, myVarImport)
myfilename <- paste0(wdir, "/SDM_results/", project,"_VarImport_single_SDMs.csv")
write.csv(myVarImport, myfilename)

# plot 
ggplot(myVarImport, aes(x = expl.var, y = var.imp)) +
  geom_boxplot(fill = "gray",outlier.size = 0.5) +
  geom_point(alpha = 0.3)+
  labs(x = "Environmental Variable", y = "Variable Importance") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "noun"
  )
ggsave(paste0(wdir, "/SDM_results/", project,"_VarImport_summary.png"), height = 6, width = 10, dpi = 300)
###############################################################################
# Ensemble the resulting models
###############################################################################

myBiomodEM <- BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModelOut,
  models.chosen = 'all',
  em.by = 'all',
  em.algo = c('EMmean', 'EMca'),
  metric.select = 'AUCroc',
  metric.select.thresh = 0.70,
  metric.eval = c('TSS', 'ROC'),
  EMci.alpha = 0.05,
  EMwmean.decay = "proportional",
  var.import = 0
)
save(myBiomodEM, file = paste0(wdir, "/SDM_results/myBiomodEM.RData"))
#load(paste0(wdir, "/SDM_results/myBiomodEM.RData"))
# ensemble_modeling_outputs
# get evaluation scores
testEM <- get_evaluations(myBiomodEM)                                                 
myfilename <- paste0(wdir, "/SDM_results/", project, "_test_ensemble.csv")
write.csv(testEM, myfilename)

###############################################################################
# get response plot
###############################################################################

my_EM_models <- BIOMOD_LoadModels(myBiomodEM)


png(paste0(wdir, "/SDM_results/", project, "_Response_Plot_2D.png"), 
    width = 20, height = 20, res = 300, units = "cm")
response.plot2 <- bm_PlotResponseCurves(
  bm.out = myBiomodEM,
  models.chosen = "all",
  new.env = get_formal_data(myBiomodModelOut, "expl.var"),
  show.variables = get_formal_data(myBiomodModelOut, "expl.var.names"),
  fixed.var = "mean",
  do.bivariate = FALSE,
  do.plot = TRUE,
  do.progress = TRUE
)
dev.off()

