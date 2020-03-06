###############################################################################
# How many MC simulations are needed to ensure stable performance results?
###############################################################################

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
simPar <- read.csv(here("data/baseSimPar.csv"), stringsAsFactors = F)
simPar <- simPar[simPar$scenario == "baseGreen",]

# Run three different "chains" of 1000 MCMC iterations to see how they converge
nChains <- 6

# Run 8000 iterations within each "chain"
nSim <- 10000

# Parallelize over 3 cores - 1 for each chain
registerDoParallel(cores=nChains)

ptime <- system.time({ 
	
	perfOut <- foreach (k = 1:nChains) %dopar% {
		
		perfL <- list(
			avgS = rep(NA, nSim),
			SR = matrix(NA, nrow = nSim, ncol = 3),
			HS = matrix(NA, nrow = nSim, ncol = 3))
		
		perf <- perfL # list(all = perfL, a = perfL, b = perfL, true.data = perfL)
		
		for(i in 1:nSim){
			out <- reconstrSim(simPar)
			
			# Only using full observation model
			perf[[1]][i] <- out$performance$currentMPE
			for(b in 2:3){ # SR benchmarks (b = 2), or HS benchmarks (b = 3)
				perf[[b]][i, 1:2] <- out$performance$benchMPE[b-1, ]
				perf[[b]][i, 3] <- out$performance$statusDiff[b-1]
			}
			
			# If using "all three" types of observed data
			# for(j in 1:4){ # for all three types of "observed data"
			# avgS (b = 1)
			# perf[[j]][[1]][i] <- out$performance[[j]]$currentMPE
			# for(b in 2:3){ # SR benchmarks (b = 2), or HS benchmarks (b = 3)
			# 	perf[[j]][[b]][i, 1:2] <- out$performance[[j]]$benchMPE[b-1, ]
			# 	perf[[j]][[b]][i, 3] <- out$performance[[j]]$statusDiff[b-1]
			# }
			# } #end j
		} # end nSim
		
		perf
	}
	
})[3]

ptime/60

# 10 mins for 10,000 MC iterations with decline in capacity etc. (March 29, 2019)

#-------------------------------------------------------------
# How does the bias evolve over the number of simulations

# Create function to return a cumulative summation even when there's NAs in the data
# Just treat NAs as zeros
# This is relevant for the SR benchmarks, since SGEN1 or SMSY may be NA under certain Ricker parameterisations
cumsumNA <- function(x){
	miss <- is.na(x)
	x[miss] <- 0
	cs <- cumsum(x)
	cs[miss] <- NA
	return(cs)
}

# What is the number of simulations that leads to perfomance metrics within 3%?
nSimEnough <- 6000
# Performance metrics are:
#		(1) MPE in avgS, 
#   (2) MPE in lowerSRbench, 
#   (3) MPE in upperSRbench, 
#   (4) MPE in lowerHSbench, 
#   (5) MPE in upperHSbench, 
#   (6) ppnWrong by SR,
#   (7) ppnWrong by HS.
cumMetric <- list(); length(cumMetric) <- 7

for(i in 1:7){ # for each metric
	cumMetric[[i]] <- matrix(NA, nrow = nSim, ncol = nChains)
}

for(k in 1:nChains) {
	# avgS
	cumMetric[[1]][, k] <- cumsum(perfOut[[k]][[1]])/c(1:nSim)
	# SRbench
	cumMetric[[2]][, k] <- cumsumNA(perfOut[[k]][[2]][, 1])/c(1:nSim)
	cumMetric[[3]][, k] <- cumsumNA(perfOut[[k]][[2]][, 2])/c(1:nSim)
	# HS bench
	cumMetric[[4]][, k] <- cumsum(perfOut[[k]][[3]][, 1])/c(1:nSim)
	cumMetric[[5]][, k] <- cumsum(perfOut[[k]][[3]][, 2])/c(1:nSim)
	# ppnWrong
	cumMetric[[6]][, k] <- cumsum(perfOut[[k]][[2]][, 3]>=4)/c(1:nSim)
	cumMetric[[7]][, k] <- cumsum(perfOut[[k]][[3]][, 3]>=4)/c(1:nSim)
}


# Calculate the maximum percent error among different "chains"
percError <- matrix(NA, nrow = nSim, ncol = 7)
for (i in 1:nSim) {
	for (j in 1:7) {
		percError[i, j] <- (max(cumMetric[[j]][i, ]) - min(cumMetric[[j]][i, ]))
	}}

# If we want to keep this percent error below 3%, how many simulations do we need?
nSimEnough <- numeric(7)
for(j in 1:7){
	# nSimEnough[j] <- which(cumsumNA(percError[, j] > 0.03) == max(cumsumNA(percError[, j] > 0.03), na.rm=TRUE))[1]
	nSimEnough[j] <- min(which(percError[, j] < 0.03))[1]
	
}

plot(1:nSim, percError[,1], "n", ylim=c(0, 0.2), ylab = "Max. % difference among 6 runs", las=1, bty="l", xlab="Number of simulations", yaxs="i")
for(j in 1:7) lines(seq(1, nSim, 10), percError[seq(1, nSim, 10),j], col=j)
abline(h = c(0.03), lty=3)

arrows(x0=nSimEnough, x1 = nSimEnough, y0 = -0.01, y1 = 0, length=0.08, col=c(1:7), xpd=NA)
max(nSimEnough)

legend("topright", title="Performance measure", col=1:7, lwd=1, c("MPE in AvgS", "MPE in lower SR", "MPE in upper SR", "MPE in lower HS", "MPE in upper HS", "Ppn. wrong by SR", "Ppn. wrong by HS"), bty="n")

# Plot MPE
par(mfcol=c(2,2), mar=c(3,5,2,2), oma=c(2,2,2,0), mgp=c(3,1,0))

i <- 4
# Plot variability ammong chains
plot(c(1:nSim), cumMetric[[i]][, 1], "n", xlab="", ylab="", bty="l",  las=1, ylim=range(cumMetric[[i]][10:nSim, ]), yaxt="n")
A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 1, format="f"), "%", sep=""), las=1)
for(k in 1:nChains) lines(c(1:nSim), cumMetric[[i]][, k], col=k, lwd=0.8)

