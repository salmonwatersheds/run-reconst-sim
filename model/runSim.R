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

source("model/plottingFns.R")

#Temporary inputs
# here <- here::here
simPar <- read.csv(here("data/baseSimPar.csv"), stringsAsFactors = F)

# set.seed(987) #a1
# # set.seed(98567) #a2
# # set.seed(6123478) #a3
# 
# a <- rnorm(simPar$nPop, simPar$a_mean, simPar$sigma_a)

# Plot three different "draws" of a parameter. Does it make a difference for performace results?
# par(mfrow=c(1,1), mar=c(4,4,2,1))
# aRange <- seq(1, 2, 0.02)
# plot(aRange, dnorm(aRange, simPar$a_mean, simPar$sigma_a), "l", las=1, ylab="Density", xlab="Productivity parameter (a)", bty="l")
# points(a1, rep(0, simPar$nPop), col=2, pch="I")
# points(a2, rep(0.15, simPar$nPop), col=3, pch="I")
# points(a3, rep(0.3, simPar$nPop), col=4, pch="I")
# legend("topright", pch="I", col=c(2:4), bty="n", c("First draw (a1)", "Second draw (a2)", "Third draw (a3)"))

# # Which is faster, rbinom or sample?
# n <- 10^6
# system.time(sample(c(0,1), size = n, replace = TRUE, prob = c(0.4, 0.6)))
# system.time(rbinom(n, size = 1, prob = 0.6))
# # Sample is clearly faster

###############################################################################
# Base case simulation
###############################################################################

# Run three different "chains" of 1000 MCMC iterations to see how they converge
nChains <- 4

# Run 8000 iterations within each "chain"
nSim <- 10000

# Parallelize over 3 cores - 1 for each chain
registerDoParallel(cores=nChains)

ptime <- system.time({ 
	
	perfOut <- foreach (k = 1:nChains) %dopar% {
		
		perfL <- list(
			SR = matrix(NA, nrow = nSim, ncol = 3),
			HS = matrix(NA, nrow = nSim, ncol = 3))
		
		perf <- list(all = perfL, a = perfL, b=perfL)
		
		for(i in 1:nSim){
			out <- reconstrSim(simPar)
			for(j in 1:3){ # for all three types of "observed data"
				for(b in 1:2){ # SR benchmarks (b = 1) or HS benchmarks (b = 2)
				perf[[j]][[b]][i, 1:2] <- out$performance[[1]]$MPE[b, ]
				perf[[j]][[b]][i, 3] <- out$performance[[1]]$statusDiff[b]
			}}
		}

		perf
	}

})[3]

ptime/60

# 8 mins for 8000 MC iterations

#-------------------------------------------------------------
# How does the bias evolve over the number of simulations
comparison <- 1 # Which set of 'observed' data to use. 1 = all, 2 = perfect obs, 3 = complete coverage
nSimEnough <- 6000

par(mfcol=c(2,2), mar=c(3,4,2,2), oma=c(2,2,2,0), mgp=c(3,1,0))

for(i in 1:2){ # for both metrics
	for(j in 1:2){ # for upper and lower bounds
		s <- 10
		Y <- c(
			cumsum(perfOut[[1]][[comparison]][[i]][s:nSim, j])/c(s:nSim), 
			cumsum(perfOut[[2]][[comparison]][[i]][s:nSim, j])/c(s:nSim), 
			cumsum(perfOut[[3]][[comparison]][[i]][s:nSim, j])/c(s:nSim))
		
		plot(c(1:nSim), cumsum(perfOut[[1]][[comparison]][[i]][, j])/c(1:nSim), "l", xlab="", ylab="", bty="l",  las=1, ylim=range(Y))
		for(k in 2:3) lines(c(1:nSim), cumsum(perfOut[[k]][[comparison]][[i]][, j])/c(1:nSim), lty=k)
		abline(v = nSimEnough, lty=2, col=4)
		
		y <- mean(mean(perfOut[[1]][[comparison]][[i]][1:nSimEnough, j]), mean(perfOut[[2]][[comparison]][[i]][1:nSimEnough, j]), mean(perfOut[[3]][[comparison]][[i]][1:nSimEnough, j]))
		abline(h =  y, col=2)
		text(nSimEnough + 1000, y, paste(formatC(round(y*100, 1), 1, format="f"), "%"), col=2, font=2, xpd=NA, pos=1)
		
		mtext(side=3, line=0.5, adj=0, matrix(c("a) Upper benchmark (Smsy)",  "c) Lower benchmark (Sgen1)", "b) Upper benchmark (S75th)",  "d) Lower benchmark (S25th)"), 2, 2)[j,i])
		
		if(j==1) mtext(side=3, line=2.5, c("Spawner-recruitment", "Historic spawners")[i], font=2)
	}}
mtext(side=2, outer=TRUE, "Mean percent error (MPE)")
mtext(side=1, "Number of simulations", outer=TRUE)

#-------------------------------------------------------------
# Proportion of simulations that are wrong

ppnWrong <- c(SR = mean(c(sum(perfOut[[1]][[comparison]]$SR[1:4000, 3]>=4)/4000, sum(perfOut[[2]][[comparison]]$SR[1:4000, 3]>=4)/4000, sum(perfOut[[3]][[comparison]]$SR[, 3]>=4)/nSim)), HS = mean(c(sum(perfOut[[1]][[comparison]]$HS[1:4000, 3]>=4)/4000, sum(perfOut[[2]][[comparison]]$HS[1:4000, 3]>=4)/4000, sum(perfOut[[3]][[comparison]]$HS[1:4000, 3]>=4)/4000)))

par(mfrow=c(1,2), mar=c(4,3,2,1), oma=c(0,3,0,0))
plot(c(1:nSim), cumsum(perfOut[[1]][[comparison]]$SR[, 3]>=4)/c(1:nSim), "l", ylim=c(0, 0.1), main="Stock-recruitment", xlab="Number of simulations", bty="l", las=1, ylab="")
mtext(side=2, "Proportion of simulations\nwith wrong status", line=3.5)
for(j in 2:3) lines(c(1:nSim), cumsum(perfOut[[j]][[comparison]]$SR[, 3]>=4)/c(1:nSim), lty=j)
abline(h = ppnWrong['SR'], col=2)
text(nSim-50, 0.08, pos=1, paste(formatC(round(ppnWrong['SR']*100, 1), 1, format="f"), "%"), col=2, font=2, xpd=NA)
abline(v = 4000, lty=2, col=4)

plot(c(1:nSim), cumsum(perfOut[[1]][[comparison]]$HS[, 3]>=4)/c(1:nSim), "l", ylim=c(0.1, 0.2), main="Historic spawners", xlab="Number of simulations", bty="l", las=1, ylab="")
for(j in 2:3) lines(c(1:nSim), cumsum(perfOut[[j]][[comparison]]$HS[, 3]>=4)/c(1:nSim), lty=j)
abline(h = ppnWrong['HS'], col=2)
text(nSim-50, ppnWrong['HS'], pos=1, paste(formatC(round(ppnWrong['HS']*100, 1), 1, format="f"), "%"), col=2, font=2, xpd=NA)
abline(v = 4000, lty=2, col=4)

#-------------------------------------------------------------
# The misclassification of status
par(mfrow=c(3,2))
for(j in 1:3){
	plotStatusDiff(perfOut[[j]][[comparison]]$SR[1:4000,3])
	mtext(side=3, "Stock-recruitment", font=2)
	
	plotStatusDiff(perfOut[[j]][[comparison]]$HS[1:4000,3])
	mtext(side=3, "Historic spawners", font=2)
}

# What is the impact of different "a" draws on conclusions?
# ppnWrong_a <- data.frame(
# 	a = rep(1:3, each = 2*3),
# 	metric = rep(rep(c("SR", "perc"), each = 3), 3),
# 	chain = rep(c(1:3), 2*3),
# 	ppnWrong = rep(NA, 2*3*3)
# )
# for(i in 1:2){ # for each metric
# 	for(j in 1:3){ # for each chain
# 		ppnWrong_a[ppnWrong_a$metric == c("SR", "perc")[i] & ppnWrong_a$chain == j, 4] <- c(sum(perfOut_a1[[j]][[i]][1:4000, 3]>=4)/4000, sum(perfOut_a2[[j]][[i]][1:4000, 3]>=4)/4000, sum(perfOut_a3[[j]][[i]][1:4000, 3]>=4)/4000)
# 	}}
# 
# fit <- lm(ppnWrong ~ a + chain + metric, data=ppnWrong_a)
# Very little impact: no significant effect of a or chain.

###############################################################################
# Over increasing bias in observation
###############################################################################

# Base case obs_bias = log(1/1.5) = - 0.4
obsBias <- log(1/seq(1, 5, 0.2)) 

# Run 1000 iterations within each "chain"
nSim <- 4000

# Parallelize over multiple cores
registerDoParallel(cores = 7)

ptime <- system.time({ 
	
	out_obsBias <- foreach (i = 1:length(obsBias)) %dopar% {
		
		simPar.i <- simPar
		simPar.i$obs_bias <- obsBias[i]
		
		perf <- list(
			SR = matrix(NA, nrow = nSim, ncol = 3),
			perc = matrix(NA, nrow = nSim, ncol = 3)
		)
		
		for(i in 1:nSim){
			out <- reconstrSim(simPar.i, a)
			for(b in 1:2){ # SR benchmarks (b = 1) or perc benchmarks (b = 2)
				perf[[b]][i, 1:2] <- out$performance$MPE[b, ]
				perf[[b]][i, 3] <- out$performance$statusDiff[b]
			}
		}
		
		perf
	}
	
})[3]

ptime/60
# 13.69 mins on 7 cores for 21 values

#----------------------------------------------------
# MPE in benchmarks

MPE_obsBias <- list(SR = matrix(NA, nrow = length(obsBias), ncol = 2), HS = matrix(NA, nrow = length(obsBias), ncol = 2, dimnames = list(NULL, c("SR", "HS"))))
for (i in 1:length(obsBias)) {
	MPE_obsBias[[1]][i, ] <- apply(out_obsBias[[i]]$SR[, 1:2], 2, mean)
	MPE_obsBias[[2]][i, ] <- apply(out_obsBias[[i]]$HS[, 1:2], 2, mean)
}

par(mfcol=c(1,2), mar=c(4,4,2,2), oma=c(0,0,0,0), mgp=c(3,1,0))
plot(obsBias, MPE_obsBias[[1]][, 1], "l", xlab="Bias in observation error", ylab="Mean percent error", main="Stock-recruitment", bty="l", las=1, col=2, ylim=c(-0.8, 0.6))
lines(obsBias, MPE_obsBias[[1]][, 2], lty=2, col=2)
legend(-1.6, 0.6, col=2, lty=c(1,2), c("Lower (Sgen)", "Upper (Smsy)"), bty="n")
abline(v = log(1/1.5), lty=2)
abline(v = 0, lty=3)
abline(h=0)

plot(obsBias, MPE_obsBias[[2]][, 1], "l", xlab="Bias in observation error", ylab="Mean percent error", main="Historic spawners", bty="l", las=1, col=4, ylim=c(-0.8, 0.6))
lines(obsBias, MPE_obsBias[[2]][, 2], lty=2, col=4)
legend(-1.6, 0.6, col=4, lty=c(1,2), c("Lower (S25th)", "Upper (S75th)"), bty="n")
abline(v = log(1/1.5), lty=2)
abline(v = 0, lty=3)
abline(h=0)


#----------------------------------------------------
# Proportion of simulations with wrong status

ppnWrong_obsBias <- matrix(NA, nrow = length(obsBias), ncol = 2, dimnames = list(NULL, c("SR", "perc")))
for (i in 1:length(obsBias)) {
	ppnWrong_obsBias[i,1] <- sum(out_obsBias[[i]]$SR[, 3]>=4)/nSim
	ppnWrong_obsBias[i,2] <- sum(out_obsBias[[i]]$perc[, 3]>=4)/nSim
}

par(mfrow = c(1,1), mar=c(4,4,4,1), oma=c(0,0,0,0))
plot(obsBias, ppnWrong_obsBias[, 1], "l", ylim = range(ppnWrong_obsBias), col = 2, lwd=1.5, las=1, bty="l", ylab = "Proportion of sim. with wrong status", xlab = "Bias in observation error")
lines(obsBias, ppnWrong_obsBias[, 2],  col=4, lwd=1.5)
legend('bottomleft', col=c(2,4), title="Metric", c("Stock-recruitment", "Historic spawners"), lwd=1.5, bty="n")
axis(side=3, at=obsBias[seq(1,21,5)], labels=seq(1, 5, 0.2)[seq(1,21,5)])
mtext(side = 3, "'true' Expansion Factor III", line=2.5)
abline(v = log(1/1.5), lty=2)
