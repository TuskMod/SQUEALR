# Spatial QUantitative Evaluation for ASFv eLimination (SQUEAL) model R package

## Background

## Goals
The goal of this project is to predict ASFv epidemic dynamics in feral pig populations in the contiguous United States following a hypothetical introduction under a range of feasible conditions. Specifically, we ask whether feral pig density, inter-sounder contact, and ASFv virulence interact with sounder movement patterns to affect ASFv establishment and spread rates.

## Model Overview
The simulations use feral pig habitat preference and sounder movement models to simulate sounder movement in a subset of real-work landscapes. These simulations include introduction of ASFv at a single point in the landscape. Transmission is modeled based on inter-sounder contact and contact with infected carcasses, a major source of persistent ASFv presence on landscapes.
Simulations include a range of feasible inter-sounder contact, ASFv virulence, and feral pig population density parameters to test sensitivity to these factors and address unknowns that might affect real-world introductions. Measured outputs are probability of establishment, maximum incidence, epidemic geographic coverage, epidemic wave speed, escape time (5+% of ASFv-affected cells outside of 10km radius from introudction point), and total disease burden.

### Habitat Preference, Sounder Movement, and Landscape Tile Selection
We model feral pig habitat preference
by connecting GPS collar data with habitat attributes in those
locations and extrapolate to locations for which no feral pig movement data is
available. 
We similarly use GPS collar data to associate movement patterns with landscape attributes such as ruggedness, tree cover, masting species presence, roads, etc. to obtain gamma density dispersion parameters for each landscape tile, aggregated from predictions for HUC12 watersheds.
We then join habitat preference with feral swine
movement for 430 100x100 km landscape tiles covering historic or current feral pig ranges in the
contiguous United States. Fifty of these tiles are selected based on a classified Latin hypercube algorithm to optimally represent the ranges of landscape attributes where feral pigs are or could be found. 

### Simulations
We simulate feral swine movement and population
dynamics, as well as ASFv transmission, on selected landscapes for 78 simulated weeks (1.5 years), with 100 replicates per combination of landscape tile, feral pig density, inter-sounder contact parameters, and ASFv virulence. The simulations include pig reproduction from healthy pigs; division of sounders if they reach a large size; movement between 0.5x0.5 km spatial grid cells; and ASFv transmission within sounders, between sounders, and from infected carcasses to living individuals, both within and between spatial grid cells.

Infection status is handled by an augmented SIR model (SEIRCZ -- susceptible, exposed, infected, recovered, dead-infected, and dead-uninfected) handled as counds of individual pigs in each sounder.

### Analysis
Model outputs are analyzed using generalized linear mixed models with appropriate forms for each measured output variable, e.g. establishment probability is modeled using a beta distribution. GLMM's are applied to the entire output dataset for inference, and for prediction we use leave-one-out cross validation with the same GLMM's applied to all but one landscape tile (with each simulation landscape tile having a round as the "left out" data subset) to determine the sensitivity of output variables to specific datasets. The predicted relationships are used to extrapolate beyond the simulated landscape tiles to locations across the range of feral pigs.


## Potential hanging threads:
input landscapes not transferred with git
may have issues in landscape tile movement gamma distribution definitions

## Notes
This model is adapted from the ASF simulation model from Pepin et al. 2022, Optimizing response to an introduction of African Swine Fever in wild pigs, converted from Matlab to R/C++.    

Forked from kchalkowski/ASF_optimal_radius and now detached for R package build
Full history is available up to the beginning of this repo in ESodja/ASF_optirad fork off of kchalkowski/ASF_optimal_radius



