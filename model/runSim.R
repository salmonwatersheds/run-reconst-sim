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
# Base case simulations
###############################################################################

a_mean <- c(0.77, 1.4, 1.97) 
simPar_a <- makeParList(basePar = simPar, sensName = "a_mean", sensValues = a_mean)

out_base <- runSensitivity(parList  = simPar_a, nSim = nSim, nCores = length(a_mean))
# 3.45 mins for 3 par sets on 3 cores

time_base <- out_base[[2]]
out_base <- out_base[[1]]

# With 100% monitoring
simPar_base2 <- list(simPar, simPar, simPar, simPar)

# a) 100% monitoring of indicator
simPar_base2[[1]]['ppnSampled_ind'] <- 1
simPar_base2[[1]]['ppnChange_ind'] <- 0

# b) 100% monitoring of non indicator
simPar_base2[[2]]['ppnSampled_nonInd'] <- 1
simPar_base2[[2]]['ppnChange_nonInd'] <- 0

# c) 100% monitoring of indicator and non-indicator
simPar_base2[[3]] <- simPar_base2[[1]]
simPar_base2[[3]]['ppnSampled_nonInd'] <- 1
simPar_base2[[3]]['ppnChange_nonInd'] <- 0

# c) 100% monitoring and no Expansion Factor III

out_base2 <- runSensitivity(parList  = simPar_base2, nSim = nSim, nCores = 3)

time_base <- out_base[[2]]
out_base2 <- out_base2[[1]]

###############################################################################
# Monitoring coverage scenarios
###############################################################################

# Base case: decline in capacity for 50% of streams
# Extreme case: decline in capacity for 100% of streams
# Hypothesis: Effect of reduced monitoring will be exacerabted by declines in capacity

# Monitoring coverage scenarios:
# 1) No decline - keep at ppnSampled_ind = 0.762 and ppnSampled_nonInd = 0.719
# 2) Observed decline since mid 1980s across all systems - ppnChange_ind = -0.047 and ppnChange_nonInd = -0.667
# 3) Observed decline since mid 1980s over for chum
# 4) Observed decline since 2014 across all systems


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

out_mon <- runSensitivity(parList = simPar_mon, nSim = 4000, nCores = 4)

time_mon <- out_mon[[2]]
out_mon <- out_mon[[1]]

out_mon2 <- delistSensitivity(out_mon)

# saveRDS(object = out_mon2, file="workspaces/mon_delisted.rds")

# Second scenario: more extreme decline in capacity
for(i in 1:4){
	simPar_mon[[i]]['greenHab'] <- 0
	simPar_mon[[i]]['amberHab'] <- 0
	simPar_mon[[i]]['redHab'] <- 100
}

out_mon <- runSensitivity(parList = simPar_mon, nSim = 4000, nCores = 4)

time_mon <- out_mon[[2]]
out_mon <- out_mon[[1]]

out_monDecl <- delistSensitivity(out_mon)
# saveRDS(object = out_monDecl, file="workspaces/mon_delisted.rds")

#------------------------------------------------------------------------------
# "Extreme fake" monitoring scenarios...where do we get an effect!?
#------------------------------------------------------------------------------
ppnSampled <- rev(seq(0.2, 1, 0.2))
ppnChange <- -seq(0.2, 0.8, 0.2)
nPar <- length(ppnSampled) + length(ppnChange)

simPar_fakeMon <- list(); length(simPar_fakeMon) <- 5#nPar
for(i in 1:length(ppnSampled)){ #constant
	simPar_fakeMon[[i]] <- simPar
	simPar_fakeMon[[i]]$ppnChange_ind <- 0
	simPar_fakeMon[[i]]$ppnChange_nonInd <- 0
	simPar_fakeMon[[i]]$ppnSampled_ind <- ppnSampled[i]
	simPar_fakeMon[[i]]$ppnSampled_nonInd <- ppnSampled[i]
}
for(i in 1:length(ppnChange)){ #
	simPar_fakeMon[[i+length(ppnSampled)]] <- simPar_fakeMon[[1]]
	simPar_fakeMon[[i+length(ppnSampled)]]$ppnChange_ind <- ppnChange[i]
	simPar_fakeMon[[i+length(ppnSampled)]]$ppnChange_nonInd <- ppnChange[i]
}

out_fakeMon <- runSensitivity(parList = simPar_fakeMon, nSim = 100, nCores = 5)

time_fakeMon <- out_fakeMon[[2]]
out_fakeMon <- out_fakeMon[[1]]

out_fakeMon2 <- delistSensitivity(out_fakeMon)

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

time_nPop <- out_nPop[[2]] # 36 minutes!
out_nPop <- out_nPop[[1]]

out_nPop2 <- delistSensitivity(out_nPop)

# saveRDS(object = out_nPop2, file="workspaces/nPop_delisted.rds")

###############################################################################
# Over increasing bias in observation
###############################################################################
# obsBias <- log(1/seq(1, 5, 0.2))

obsBias <- seq(-1.6, 0, 0.2)
# obsBias <- seq(-3, 3, 0.5) # Wider range Apr 8

simPar_obsBias <- makeParList(basePar = simPar, sensName = "obs_bias", sensValues = obsBias)

# for(i in 1:length(obsBias)) simPar_obsBias[[i]]$sigma_obs <- 0.1

out_obsBias <- runSensitivity(parList = simPar_obsBias, nSim = 4000, nCores = 10)


time_obsBias <- out_obsBias[[2]]
out_obsBias <- out_obsBias[[1]]

out_obsBias2 <- delistSensitivity(out_obsBias)

# saveRDS(object = out_obsBias2, file="workspaces/obsBias_delisted.rds")

###############################################################################
# Over increasing bias in catch
###############################################################################
# catchBias <- log(seq(0.5, 1.5, 0.1))

catchBias <- seq(-1, 1, 0.2)
simPar_catchBias <- makeParList(basePar = simPar, sensName = "catch_bias", sensValues = catchBias)

out_catchBias <- runSensitivity(parList = simPar_catchBias, nSim = 4000, nCores = 10)

time_catchBias <- out_catchBias[[2]]
out_catchBias <- out_catchBias[[1]]

out_catchBias2 <- delistSensitivity(out_catchBias)

# saveRDS(object = out_catchBias2, file="workspaces/catchBias_delisted.rds")

# Underestimating catch (catchBias = -1) led to fewer misclassifications 
# than overestimating catch (catchBias = 1).
# Is this due to green being misclassified as amber or amber being misclassified
# as red?

greenASamber <- numeric(length(catchBias))
amberASred<- numeric(length(catchBias))
greenASred <- numeric(length(catchBias))
for (i in 1:length(catchBias)){
	greenASamber[i] <- sum(out_catchBias[[i]]$SR[, 3] == 7)/nSim
	amberASred[i] <- sum(out_catchBias[[i]]$SR[, 3] == 4)/nSim
	greenASred[i] <- sum(out_catchBias[[i]]$SR[, 3] == 5)/nSim
}

plot(catchBias, greenASamber, "l", col=statusCols['a'], lwd=2, ylim=c(0,0.4))
lines(catchBias, amberASred, col=statusCols['r'], lwd=2)
lines(catchBias, greenASred, col=statusCols['r'], lty=2)


###############################################################################
# Over increasing decline in capacity
###############################################################################
greenHab <- rev(seq(0, 100, 10))
redHab <- seq(0, 50, 5)
amberHab <- seq(0, 50, 5)

# Check
greenHab + amberHab + redHab

# Have to make this list manually since multiple parameters change at once
simPar_capacity <- list(); length(simPar_capacity) <- length(greenHab)
for(i in 1:length(greenHab)){
	simPar_capacity[[i]] <- simPar
	simPar_capacity[[i]][which(names(simPar) == "greenHab")] <- greenHab[i]
	simPar_capacity[[i]][which(names(simPar) == "amberHab")] <- amberHab[i]
	simPar_capacity[[i]][which(names(simPar) == "redHab")] <- redHab[i]
}

out_capacity <- runSensitivity(parList = simPar_capacity, nSim = 4000, nCores = 10)

time_capacity <- out_capacity[[2]]
out_capacity <- out_capacity[[1]]

out_capacity2 <- delistSensitivity(out_capacity)
# saveRDS(object = out_capacity2, file="workspaces/capacity_delisted.rds")

greenASamber <- numeric(length(catchBias))
amberASred<- numeric(length(catchBias))
greenASred <- numeric(length(catchBias))

for (i in 1:length(greenHab)){
	greenASamber[i] <- sum(out_capacity[[i]]$HS[, 3] == 7)/nSim
	amberASred[i] <- sum(out_capacity[[i]]$HS[, 3] == 4)/nSim
	greenASred[i] <- sum(out_capacity[[i]]$HS[, 3] == 5)/nSim
}

plot(amberHab + redHab, greenASamber, "l", col=statusCols['a'], lwd=2, ylim=c(0,0.4))
lines(amberHab + redHab, amberASred, col=statusCols['r'], lwd=2)
lines(amberHab + redHab, greenASred, col=statusCols['r'], lty=2)

###############################################################################
# Monitoring decline + decline in capacity
###############################################################################
greenHab <- rev(seq(0, 100, 20))
redHab <- seq(0, 50, 10)
amberHab <- seq(0, 50, 10)

ppnSampled <- c(0.8, 0.4, 0.8, 0.8)
ppnChange <- c(0, 0, -0.4, -0.4)
samplingDeclStart <- c(NA, NA, 1, 41)

nPar <- length(greenHab)*length(ppnSampled)
combo <- cbind(rep(1:length(greenHab), length(ppnSampled)), rep(1:length(ppnSampled), each = length(greenHab)))

# Have to make this list manually since multiple parameters change at once
simPar_capacityXmon <- list(); length(simPar_capacityXmon) <- nPar
for(i in 1:nPar){
	simPar_capacityXmon[[i]] <- simPar
	
	simPar_capacityXmon[[i]][which(names(simPar) == "greenHab")] <- greenHab[combo[i,1]]
	simPar_capacityXmon[[i]][which(names(simPar) == "amberHab")] <- amberHab[combo[i,1]]
	simPar_capacityXmon[[i]][which(names(simPar) == "redHab")] <- redHab[combo[i,1]]
	
	simPar_capacityXmon[[i]][which(names(simPar) == "ppnSampled_ind")] <- ppnSampled[combo[i,2]]
	simPar_capacityXmon[[i]][which(names(simPar) == "ppnSampled_nonInd")] <- ppnSampled[combo[i,2]]
	simPar_capacityXmon[[i]][which(names(simPar) == "ppnChange_ind")] <- ppnChange[combo[i,2]]
	simPar_capacityXmon[[i]][which(names(simPar) == "ppnChange_nonInd")] <- ppnChange[combo[i,2]]
	simPar_capacityXmon[[i]][which(names(simPar) == "samplingDeclStart_ind")] <- samplingDeclStart[combo[i,2]]
	simPar_capacityXmon[[i]][which(names(simPar) == "samplingDeclStart_nonInd")] <- samplingDeclStart[combo[i,2]]
}

out_capacityXmon <- runSensitivity(parList = simPar_capacityXmon, nSim = 4000, nCores = 10)

time_capacityXmon <- out_capacityXmon[[2]]
out_capacityXmon <- out_capacityXmon[[1]]

out_capacityXmon2 <- delistSensitivity(out_capacityXmon)
# saveRDS(object = out_capacityXmon2, file="workspaces/capacityXmon_delisted.rds")
