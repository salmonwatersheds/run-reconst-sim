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
source("runSims/runSensitivity.R")

# Load base parameter values
simPar_all <- read.csv(here("data/baseSimPar.csv"), stringsAsFactors = F)

# Determined that 4000 simulations is adequate
# See nSim.R for justification.

nSim <- 4000

# Three base cases:
# 1) high productivity, HCR -> true green
# 2) low productivity, target harvest 0.42
# 3) low productivity, target harvest 0.60

# For each of the sensitivity analyses, we want to produce results
# for each of these base cases.
# Loop though base case scenarios for main results and save output
# as RDS.

for(baseCaseNum in 1:3){
	
	if(baseCaseNum == 1) simPar <- simPar_all[simPar_all$scenario == "baseGreen",]
	if(baseCaseNum == 2) simPar <- simPar_all[simPar_all$scenario == "baseAmber",]
	if(baseCaseNum == 3) simPar <- simPar_all[simPar_all$scenario == "baseRed",]
	
	###############################################################################
	# Base case simulations
	###############################################################################
	#
	# Test expansion factors I and II by also simulating with 100% monitoring
	simPar_base100Mon <- list(simPar, simPar, simPar, simPar)

	# a) Base case

	# b) 100% monitoring of indicator
	simPar_base100Mon[[2]]['ppnSampled_ind'] <- 1
	simPar_base100Mon[[2]]['ppnChange_ind'] <- 0

	# b) 100% monitoring of non indicator
	simPar_base100Mon[[3]]['ppnSampled_nonInd'] <- 1
	simPar_base100Mon[[3]]['ppnChange_nonInd'] <- 0

	# c) 100% monitoring of indicator and non-indicator
	simPar_base100Mon[[4]] <- simPar_base100Mon[[2]]
	simPar_base100Mon[[4]]['ppnSampled_nonInd'] <- 1
	simPar_base100Mon[[4]]['ppnChange_nonInd'] <- 0

	out_base100Mon <- runSensitivity(parList  = simPar_base100Mon, nSim = nSim, nCores = 4)

	print(paste("time for 100Mon = ", out_base100Mon[[2]]))
	out_base100Mon2 <- delistSensitivity(out_base100Mon[[1]])

	if(baseCaseNum == 1) saveRDS(object = out_base100Mon2, file="workspaces/base100Mon_delisted_baseGreen.rds")
	if(baseCaseNum == 2) saveRDS(object = out_base100Mon2, file="workspaces/base100Mon_delisted_baseAmber.rds")
	if(baseCaseNum == 3) saveRDS(object = out_base100Mon2, file="workspaces/base100Mon_delisted_baseRed.rds")

}

	###############################################################################
	# Real monitoring scenarios x decline in capacity ** Fig 3 **
	###############################################################################
	# Monitoring coverage scenarios:
	# 1) No decline - keep at ppnSampled_ind = 0.762 and ppnSampled_nonInd = 0.719
	# 2) Observed decline since mid 1980s across all systems - ppnChange_ind = -0.047 and ppnChange_nonInd = -0.667
	# 3) Observed decline since mid 1980s over for chum
	# 4) Observed decline since 2014 across all systems

	for(baseCaseNum in 1:3){
		
		if(baseCaseNum == 1) simPar <- simPar_all[simPar_all$scenario == "baseGreen",]
		if(baseCaseNum == 2) simPar <- simPar_all[simPar_all$scenario == "baseAmber",]
		if(baseCaseNum == 3) simPar <- simPar_all[simPar_all$scenario == "baseRed",]
		simPar$correlPop <- 0.9
		
		simPar_mon <- list(); length(simPar_mon) <- 4
	
	
	# 1) No decline
	simPar_mon[[1]] <- simPar
	simPar_mon[[1]]['ppnChange_ind'] <- 0
	simPar_mon[[1]]['ppnChange_nonInd'] <- 0

	# 2) Base case
	simPar_mon[[2]] <- simPar

	# 3) Chum decline
	simPar_mon[[3]] <- simPar
	simPar_mon[[3]]['ppnChange_ind'] <- -0.2029197
	simPar_mon[[3]]['ppnChange_nonInd'] <- -0.2433794
	simPar_mon[[3]]['ppnSampled_ind'] <- 0.8604703
	simPar_mon[[3]]['ppnSampled_nonInd'] <- 0.3323478

	# 4) Recent decline for indicator streams
	simPar_mon[[4]] <- simPar
	simPar_mon[[4]]['ppnSampled_ind'] <- 0.7720889
	simPar_mon[[4]]['ppnChange_ind'] <- -0.208303
	simPar_mon[[4]]['samplingDeclStart_ind'] <- 46

	# Changes in capacity
	greenHab <- c(rev(seq(0, 100, 50)), 0)
	redHab <- c(seq(0, 50, 25), 100)
	amberHab <- c(seq(0, 50, 25), 0)

	nPar <- length(simPar_mon) * length(greenHab)

	combo <- cbind(cap = rep(1:length(greenHab), length(simPar_mon)), mon = rep(1:length(simPar_mon), each = length(greenHab)))

	# Have to make this list manually since multiple parameters change at once
	simPar_capacityXmon <- list(); length(simPar_capacityXmon) <- nPar
	for(i in 1:nPar){
		simPar_capacityXmon[[i]] <- simPar_mon[[combo[i,'mon']]]

		simPar_capacityXmon[[i]][which(names(simPar) == "greenHab")] <- greenHab[combo[i,'cap']]
		simPar_capacityXmon[[i]][which(names(simPar) == "amberHab")] <- amberHab[combo[i,'cap']]
		simPar_capacityXmon[[i]][which(names(simPar) == "redHab")] <- redHab[combo[i,'cap']]
	}

	out_capacityXmon <- runSensitivity(parList = simPar_capacityXmon, nSim = 4000, nCores = 8)

	print(paste("time for capacityXmon = ", out_capacityXmon[[2]]))

	out_capacityXmon2 <- delistSensitivity(out_capacityXmon[[1]])

	if(baseCaseNum == 1) saveRDS(object = out_capacityXmon2, file="workspaces/capacityXmon_delisted_baseGreen_correlPop09.rds")
	if(baseCaseNum == 2) saveRDS(object = out_capacityXmon2, file="workspaces/capacityXmon_delisted_baseAmber_correlPop09.rds")
	if(baseCaseNum == 3) saveRDS(object = out_capacityXmon2, file="workspaces/capacityXmon_delisted_baseRed_correlPop09.rds")
}

	###############################################################################
	# Number of indicator & non-indicator streams
	###############################################################################
	
		#Cases:
	# 1) Base: 15 ind + 20 non-ind = 35
	# 2) SmallLow: 3 in + 7 non-ind = 10
	# 3) SmallHigh: 8 ind + 2 non-ind = 10
	# 4) LargeLow: 42 ind + 58 non-ind = 140
	# 5) LargeHigh: 119 ind + 21 non-ind = 140
	
	nInd <- c(15, 3, 8, 42, 119)
	nNonInd <- c(20, 7, 2, 58, 21)
	
	simPar_nPop <- list(simPar, simPar, simPar, simPar, simPar)
	
	for(i in 1:length(simPar_nPop)){
		simPar_nPop[[i]]$nIndicator <- nInd[i]
		simPar_nPop[[i]]$nNonIndicator<- nNonInd[i]
		simPar_nPop[[i]]$nPop<- nInd[i] + nNonInd[i]
	}
	
	out_nPop <- runSensitivity(parList = simPar_nPop, nSim = 4000, nCores = 5)
	
	print(paste("time for nPop = ", out_nPop[[2]]))
	
	out_nPop2 <- delistSensitivity(out_nPop[[1]])
	
	if(baseCaseNum == 1) saveRDS(object = out_nPop2, file="workspaces/nPop_delisted_baseGreen.rds")
	if(baseCaseNum == 2) saveRDS(object = out_nPop2, file="workspaces/nPop_delisted_baseAmber.rds")
	if(baseCaseNum == 3) saveRDS(object = out_nPop2, file="workspaces/nPop_delisted_baseRed.rds")
	
	###############################################################################
	# Over increasing bias in observation
	###############################################################################
	
		obsBias <- seq(-1.6, 0, 0.2)
	simPar_obsBias <- makeParList(basePar = simPar, sensName = "obs_bias", sensValues = obsBias)
	
	out_obsBias <- runSensitivity(parList = simPar_obsBias, nSim = 4000, nCores = 10)
	
	print(paste("time for obsBias = ", out_obsBias[[2]]))
	out_obsBias2 <- delistSensitivity(out_obsBias[[1]])
	
	if(baseCaseNum == 1) saveRDS(object = out_obsBias2, file="workspaces/obsBias_delisted_baseGreen.rds")
	if(baseCaseNum == 2) saveRDS(object = out_obsBias2, file="workspaces/obsBias_delisted_baseAmber.rds")
	if(baseCaseNum == 3) saveRDS(object = out_obsBias2, file="workspaces/obsBias_delisted_baseRed.rds")
	
	###############################################################################
	# Over increasing bias in catch
	###############################################################################
	catchBias <- seq(-1, 1, 0.2)

	simPar_catchBias <- makeParList(basePar = simPar, sensName = "catch_bias", sensValues = catchBias)

	out_catchBias <- runSensitivity(parList = simPar_catchBias, nSim = 4000, nCores = 6)
	
	print(paste("time for catchBias = ", out_catchBias[[2]]))
	out_catchBias2 <- delistSensitivity(out_catchBias[[1]])
	
	if(baseCaseNum == 1) saveRDS(object = out_catchBias2, file="workspaces/catchBias_delisted_baseGreen.rds")
	if(baseCaseNum == 2) saveRDS(object = out_catchBias2, file="workspaces/catchBias_delisted_baseAmber.rds")
	if(baseCaseNum == 3) saveRDS(object = out_catchBias2, file="workspaces/catchBias_delisted_baseRed.rds")
	
	###############################################################################
	# Changes in observation bias part-way through time series 
	###############################################################################
	
		obs_bias2 <- c(-1.6, -0.7, -0.4,  0)

	simPar_obsBiasChange <- makeParList(basePar = simPar, sensName = "obs_bias2", sensValues = obs_bias2)

	out_obsBiasChange <- runSensitivity(parList = simPar_obsBiasChange, nSim = 4000, nCores = 4)

	print(paste("time for obsBiasChange = ", out_obsBiasChange[[2]]))
	out_obsBiasChange2 <- delistSensitivity(out_obsBiasChange[[1]])
	
	if(baseCaseNum == 1) saveRDS(object = out_obsBiasChange2, file="workspaces/obsBiasChange_delisted_baseGreen.rds")
	if(baseCaseNum == 2) saveRDS(object = out_obsBiasChange2, file="workspaces/obsBiasChange_delisted_baseAmber.rds")
	if(baseCaseNum == 3) saveRDS(object = out_obsBiasChange2, file="workspaces/obsBiasChange_delisted_baseRed.rds")
	
	###############################################################################
	# Over increasing interannual variability in age-at-maturity
	###############################################################################		
for(baseCaseNum in 1:3){
		
		if(baseCaseNum == 1) simPar <- simPar_all[simPar_all$scenario == "baseGreen",]
		if(baseCaseNum == 2) simPar <- simPar_all[simPar_all$scenario == "baseAmber",]
		if(baseCaseNum == 3) simPar <- simPar_all[simPar_all$scenario == "baseRed",]
		
	ageErr <- seq(0.2, 1.6, 0.2)
	
	simPar_ageErr <- makeParList(basePar = simPar, sensName = "ageErr", sensValues = ageErr)
	
	out_ageErr <- runSensitivity(parList = simPar_ageErr, nSim = 4000, nCores = 8)
	
	print(paste("time for ageErr = ", out_ageErr[[2]]))
	out_ageErr2 <- delistSensitivity(out_ageErr[[1]])
	
	if(baseCaseNum == 1) saveRDS(object = out_ageErr2, file="workspaces/ageErr_delisted_baseGreen.rds")
	if(baseCaseNum == 2) saveRDS(object = out_ageErr2, file="workspaces/ageErr_delisted_baseAmber.rds")
	if(baseCaseNum == 3) saveRDS(object = out_ageErr2, file="workspaces/ageErr_delisted_baseRed.rds")
	
} # end loop over baseCaseNum 1-3
	
	
	###############################################################################
	# Over decreasing correlation in recruitment deviates
	###############################################################################
	for(baseCaseNum in 2:3){

		if(baseCaseNum == 1) simPar <- simPar_all[simPar_all$scenario == "baseGreen",]
		if(baseCaseNum == 2) simPar <- simPar_all[simPar_all$scenario == "baseAmber",]
		if(baseCaseNum == 3) simPar <- simPar_all[simPar_all$scenario == "baseRed",]

		correlPop <- seq(0, 0.9, 0.1)

		simPar_correlPop <- makeParList(basePar = simPar, sensName = "correlPop", sensValues = correlPop)

		out_correlPop <- runSensitivity(parList = simPar_correlPop, nSim = 4000, nCores = 8)

		print(paste("time for correlPop = ", out_correlPop[[2]]))
		out_correlPop2 <- delistSensitivity(out_correlPop[[1]])

		if(baseCaseNum == 1) saveRDS(object = out_correlPop2, file="workspaces/correlPop_delisted_baseGreen.rds")
		if(baseCaseNum == 2) saveRDS(object = out_correlPop2, file="workspaces/correlPop_delisted_baseAmber.rds")
		if(baseCaseNum == 3) saveRDS(object = out_correlPop2, file="workspaces/correlPop_delisted_baseRed.rds")

	} # end loop over baseCaseNum 1-3