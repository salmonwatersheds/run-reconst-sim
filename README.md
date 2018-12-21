# Evaluating the consequences of assumptions in run-reconstructions when assessing the biological status of salmon
## R code for model simulations
### Stephanie J. Peacock <stephanie.j.peacock@gmail.com>
### December 21, 2018

This repo contains R code related to the PSF project "Evaluating the consequences of assumptions in run-reconstructions when assessing the biological status of salmon"

## Model overview

This project will develop and apply a stochastic simulation model of salmon population dynamics that allows control over various biological and management factors that may influence the accuracy of status assessments.  This model will be based on previous studies by Carrie Holt and colleagues (e.g., [Holt and Folkes 2015](http://dx.doi.org/10.1016/j.fishres.2015.01.002), Holt et al. 2018) that developed a simulation model comprised of sub-models for salmon population dynamics, observation of spawners, assessment, harvest, and performance (Fig. 1).

![Fig. 1. Schematic of the simulation model, including sub-models for population dynamics, harvest, observation, assessment, and performance.  The entire process will be repeated for different autocorrelation in residuals among sub-populations, inter-annual variability in age-at-return, bias in harvest, and observation errors (see research questions, above). Adapted from [Holt et al. (2016)](https://www.psc.org/fund-project/adapting-benchmarks/).](model.png)

## Model equations

### Population sub-model
The population dynamics of multiple sub-populations, designated as indicator or non-indicator streams, are simulated within a single hypothetical CU.  The population dynamics follow a Ricker type stock-recruitment relationship, with parameters based on observations from central coast chum CUs ([Connors et al. 2018](https://salmonwatersheds.ca/library/lib_442/)).  Target harvest rates will be predetermined, but include year-to-year variability and potential biases in realized harvest.  