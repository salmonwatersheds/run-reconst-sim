library(mvtnorm)
library(here)
library(corpcor) # make.positive.definite function for Sigma
library(data.table) # for shift function
library(gsl) # for lambert_W0 function
library(doParallel) # for parallelizing function application over parameters or MCMC iterations


source("model/populationSubmodFns.R")
source("model/obsSubmodFns.R")
source("model/expansionFactors.R")
source("model/benchmarkFns.R")
source("model/performanceFns.R")
source("model/reconstrSimulator.R")
source("model/runSensitivity.R")
source("model/plottingFns.R")

# Load base parameter values
simPar <- read.csv(here("data/baseSimPar.csv"), stringsAsFactors = F)
simPar <- simPar[simPar$scenario == "base",]

# Determined that 4000 simulations is adequate
# See nSim.R for justification.

nSim <- 4000

###############################################################################
# Wrapper function for key params
###############################################################################
parSubset <- c(ageErr = 0.8, obs_bias = -0.4, correlPop = 0.46)

# a function returning a outputs a matrix having as columns the model outputs.
morrisSim <- function(parSubset){
	
	par <- simPar
	par[match(names(parSubset), names(simPar))] <- parSubset
	
	nSim <- 1000
	perfMeasures <- matrix(nrow = nSim, ncol = 7, dimnames = list(NULL, c("avgS", "Sgen1", "Smsy", "S25", "S50", "SRwrong", "HSwrong")))
	for(i in 1:nSim){
		out <- reconstrSim(par)
		perfMeasures[i, 1] <- out$performance$currentMPE
		perfMeasures[i, 2:3] <- out$performance$benchMPE[1, ]
		perfMeasures[i, 4:5] <- out$performance$benchMPE[2, ]
		perfMeasures[i, 6] <- as.numeric(out$performance$statusDiff[1] > 3)
		perfMeasures[i, 7] <- as.numeric(out$performance$statusDiff[2] > 3)
	}
	
	perfMean <- c(apply(perfMeasures[, 1:5], 2, mean, na.rm=TRUE), apply(perfMeasures[, 6:7], 2, sum)/nSim)
	
		
	
} # end morrisSim
