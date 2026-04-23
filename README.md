# Current and future distribution of Common European adder across southern Germany

This repository provides R code for modeling the current and future distribution of the European adder in southern Germany. It includes a complete workflow from data preparation to niche analysis and driver identification.

1. Data preparation
Occurrence records (post-1980) were spatially thinned (1 km). Environmental variables were selected using pairwise correlation (r > 0.75), Variance Inflation Factor (VIF), and Principal Component Analysis (PCA), then cropped to the study area.

2. Modeling
Species distribution models were built using three algorithms in the BIOMOD2 framework. The workflow includes pseudo-absence selection, cross-validation, model evaluation and selection, variable importance, ensemble modeling, and response curves.

3. Projection
The ensemble model was projected onto current and future climate scenarios. Outputs were rescaled (0–1), and Multivariate Environmental Similarity Surfaces (MESS) were generated to identify interpolation and extrapolation areas.

4. Post-modeling analysis
Changes in habitat suitability were quantified across scenarios using multiple thresholds. Altitudinal shifts were assessed for areas of habitat gain and loss. PCA-based analyses were used to estimate niche space, hypervolume, and overlap between scenarios.

5. Drivers of distribution
A Generalized Linear Mixed Model (GLMM) was used to identify key climatic and anthropogenic drivers of current distribution patterns, including analysis of distances to wetlands for extinct vs. extant populations.

6. Niche structure
Unsupervised clustering based on PCA was applied to identify niche substructures and ecological differentiation within the study area.
