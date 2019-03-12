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

#Temporary inputs
# here <- here::here
simPar <- read.csv(here("data/baseSimPar.csv"), stringsAsFactors = F)
simPar <- simPar[simPar$scenario == "base",]
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
# How many MC simulations are needed??
###############################################################################

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
		
		perf <- list(all = perfL, a = perfL, b = perfL, true.data = perfL)
		
		for(i in 1:nSim){
			out <- reconstrSim(simPar)
			for(j in 1:4){ # for all three types of "observed data"
				# avgS (b = 1)
				perf[[j]][[1]][i] <- out$performance[[j]]$currentMPE
				for(b in 2:3){ # SR benchmarks (b = 2), or HS benchmarks (b = 3)
					perf[[j]][[b]][i, 1:2] <- out$performance[[j]]$benchMPE[b-1, ]
					perf[[j]][[b]][i, 3] <- out$performance[[j]]$statusDiff[b-1]
			}}
		} # end nSim

		perf
	}

})[3]

ptime/60

# 8 mins for 8000 MC iterations

#-------------------------------------------------------------
# How does the bias evolve over the number of simulations
comparison <- 1 # Which set of 'observed' data to use. 1 = all, 2 = perfect obs, 3 = complete coverage

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
		cumMetric[[1]][, k] <- cumsum(perfOut[[k]][[comparison]][[1]])/c(1:nSim)
		# SRbench
		cumMetric[[2]][, k] <- cumsum(perfOut[[k]][[comparison]][[2]][, 1])/c(1:nSim)
		cumMetric[[3]][, k] <- cumsum(perfOut[[k]][[comparison]][[2]][, 2])/c(1:nSim)
		# HS bench
		cumMetric[[4]][, k] <- cumsum(perfOut[[k]][[comparison]][[3]][, 1])/c(1:nSim)
		cumMetric[[5]][, k] <- cumsum(perfOut[[k]][[comparison]][[3]][, 2])/c(1:nSim)
		# ppnWrong
		cumMetric[[6]][, k] <- cumsum(perfOut[[k]][[comparison]][[2]][, 3]>=4)/c(1:nSim)
		cumMetric[[7]][, k] <- cumsum(perfOut[[k]][[comparison]][[3]][, 3]>=4)/c(1:nSim)
	}

percError <- matrix(NA, nrow = nSim, ncol = 7)
for (i in 1:nSim) {
	for (j in 1:7) {
		percError[i, j] <- (max(cumMetric[[j]][i, ]) - min(cumMetric[[j]][i, ]))
	}}

nSimEnough <- numeric(7)
for(j in 1:7){
	nSimEnough[j] <- which(cumsum(percError[, j] > 0.031) == max(cumsum(percError[, j] > 0.031)))
}

plot(1:nSim, percError[,1], "l", ylim=c(0, 0.1), ylab = "Max. difference among 6 runs", las=1, bty="l", xlab="Number of simulations", yaxs="i")
for(j in 2:7) lines(1:nSim, percError[,j], col=j)
abline(h = 0)
abline(h = c(0.03), lty=3)

arrows(x0=nSimEnough, x1 = nSimEnough, y0 = 0.1, y1 = 0, length=0.08, col=c(1:7))
max(nSimEnough)


# Plot MPE
par(mfcol=c(2,2), mar=c(3,5,2,2), oma=c(2,2,2,0), mgp=c(3,1,0))

i <- 4
# Plot variability ammong chains
plot(c(1:nSim), cumMetric[[i]][, 1], "n", xlab="", ylab="", bty="l",  las=1, ylim=range(cumMetric[[i]][10:nSim, ]), yaxt="n")
A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 1, format="f"), "%", sep=""), las=1)
for(k in 1:nChains) lines(c(1:nSim), cumMetric[[i]][, k], col=k, lwd=0.8)
		

###############################################################################
# Base case simulations
###############################################################################

nSim <- 5000

a_mean <- c(0.4, 1.4, 2.4) # range for CCC = 1.3 - 1.5
targetHarvest <- c(0.10, 0.35, 0.80)
baseScen <- cbind(a_mean = rep(a_mean, each = length(targetHarvest)), targetHarvest = rep(targetHarvest, length(a_mean)))

simPar_base <- list(); length(simPar_base) <- (dim(baseScen)[1])
for (i in 1:(dim(baseScen)[1])) {
	simPar_base[[i]] <- simPar
	simPar_base[[i]]$a_mean <- baseScen[i,1]
	simPar_base[[i]]$targetHarvest <- baseScen[i,2]
	
}

out_base <- runSensitivity(basePar = simPar_base, sensitivity = NULL, nSim = nSim, nCores = 9, multiPar = TRUE)

time_base <- out_base[[2]]
out_base <- out_base[[1]]

#-------------------------------------------------------------
# The misclassification of status
i <- 6
comparison <- 1 #1 = all, 2 = A, 3 = B, 4 = true.data
# quartz(width = 4.5, height = 2.2, pointsize = 8)
par(mfrow=c(1,2))

plotStatusDiff(out_base[[i]][[comparison]]$HS[,3])
mtext(side=3, "Historic spawners", font=2, col=col2[1])

plotStatusDiff(out_base[[i]][[comparison]]$SR[, 3])
mtext(side=3, "Stock-recruitment", font=2, col=col2[2])


par(mfcol=c(3,3), mar=rep(0.5, 4), oma=c(5,4,2,3))
for(i in 1:9){
	j <- c(3:1, 6:4, 9:7)[i]
	
	R <- sum(
			length(which(out_base[[j]][[comparison]]$SR[, 3] == 3)), 
			length(which(out_base[[j]][[comparison]]$SR[, 3] == 6)), 
			length(which(out_base[[j]][[comparison]]$SR[, 3] == 8)))
	
	A <- sum(
			length(which(out_base[[j]][[comparison]]$SR[, 3] == 4)), 
			length(which(out_base[[j]][[comparison]]$SR[, 3] == 2)), 
			length(which(out_base[[j]][[comparison]]$SR[, 3] == 9)))
	
	G <- sum(
			length(which(out_base[[j]][[comparison]]$SR[, 3] == 5)), 
			length(which(out_base[[j]][[comparison]]$SR[, 3] == 7)), 
			length(which(out_base[[j]][[comparison]]$SR[, 3] == 1)))
	
barplot(c(G,A,R)/nSim, space = 0.1, col = c(g = "#8EB687", a = "#DFD98D", r = "#B66A64"), ylim=c(0,1), yaxt="n")
if(i <= 3) axis(side=2, at=c(0, 0.5, 1), las=1) else axis(side=2, at=c(0, 0.5, 1), labels=FALSE)

if(i == 3) mtext(side = 1, "Low\nproductivity", line=2.5)
if(i == 6) mtext(side = 1, "Moderate\nproductivity", line=2.5)
if(i == 9) mtext(side = 1, "High\nproductivity", line=2.5)

if(i == 7) mtext(side = 4, "High\nharvest", line=2)
if(i == 8) mtext(side = 4, "Moderate\nharvest", line=2)
if(i == 9) mtext(side = 4, "Low\nharvest", line=2)
}
mtext(side=2, outer=TRUE, "Proportion of simulations", line=2.5)

focalScen <- c(5,9,6)
ppnWrongSummary <- matrix(NA, nrow = 5, ncol = 6, dimnames = list(c("G", "A", "R", "Cautious", "Risky"), c("SR.base", "HS.base", "SR.amber", "HS.amber", "SR.red", "HS.red")))
for(i in 1:3){
	for(j in 1:2){
		ppnWrongSummary[1, (i-1)*2+j] <- length(which(out_base[[focalScen[i]]][[comparison]][[c(3,2)[j]]][, 3] == 1))/nSim
		ppnWrongSummary[2, (i-1)*2+j] <- length(which(out_base[[focalScen[i]]][[comparison]][[c(3,2)[j]]][, 3] == 2))/nSim
		ppnWrongSummary[3, (i-1)*2+j] <- length(which(out_base[[focalScen[i]]][[comparison]][[c(3,2)[j]]][, 3] == 3))/nSim
		
		ppnWrongSummary[4, (i-1)*2+j] <- (length(which(out_base[[focalScen[i]]][[comparison]][[c(3,2)[j]]][, 3] == 4)) + length(which(out_base[[focalScen[i]]][[comparison]][[c(3,2)[j]]][, 3] == 5)) + length(which(out_base[[focalScen[i]]][[comparison]][[c(3,2)[j]]][, 3] == 7)))/nSim
		
		ppnWrongSummary[5, (i-1)*2+j] <- (length(which(out_base[[focalScen[i]]][[comparison]][[c(3,2)[j]]][, 3] == 6)) + length(which(out_base[[focalScen[i]]][[comparison]][[c(3,2)[j]]][, 3] == 8)) + length(which(out_base[[focalScen[i]]][[comparison]][[c(3,2)[j]]][, 3] == 9)))/nSim
	}
}

par(mfrow=c(1,1), mar=c(5,4,2,15), oma=rep(0,4))
bp <- barplot(ppnWrongSummary, col=c(g = "#8EB687", a = "#DFD98D", r = "#B66A64", grey(0.8), 1), space = c(0.3, 0.05), names = rep(c("HS", "SR"), 3), las=1, ylab = "Proportion of simulations")
text(mean(bp[1:2]), -0.25, xpd=NA, "Moderate productivity\nModerate harvest", cex=0.8, font=2)
text(mean(bp[3:4]), -0.25, xpd=NA, "High productivity\nHigh harvest", cex=0.8)
text(mean(bp[5:6]), -0.25, xpd=NA, "Moderate productivity\nHigh harvest", cex=0.8)

legend(7.5, 0.8, fill=c(g = "#8EB687", a = "#DFD98D", r = "#B66A64", grey(0.8), 1), legend = c("Correctly assessed as green", "Correctly assessed as amber", "Correctly assessed as red", "Incorrect (conservative)", "Incorrect (risky)"), xpd=NA, bty="n")

# Why are there so many incorrect with HS but not SR for baseScen 5?
j <- 6

x <- data.frame(MPE = c(out_base[[j]][[comparison]]$avgS, out_base[[j]][[comparison]]$SR[,1], out_base[[j]][[comparison]]$SR[,2], out_base[[j]][[comparison]]$HS[,1], out_base[[j]][[comparison]]$HS[,2]), metric = c(rep(c("avgS", "SR.lower", "SR.upper", "HS.lower", "HS.upper"), each=nSim)))

tapply(x$MPE, x$metric, mean)

hist(x$MPE[x$metric=="HS.lower"])

col2 <- c("#203864", "#8FAADC")
par(mfrow=c(1,1), mar=c(5,6,2,1), oma=rep(0,4))
boxplot(MPE ~ metric, data=x, yaxt="n", outline=FALSE, bty="n",pars = list(boxwex = 0.5, staplewex = 0.5, outwex = 0.5), names=c("average\nspawners", "S25", "S75", "Smsy", "Sgen1"), col=c("white", rep(c(col2), each=2)))
A <- axis(side=2, labels = FALSE)
mtext(side=2, line=4, "Mean percent error")
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
abline(h = 0)

#-----------
quartz(width = 6.3, height = 3.3, pointsize=10)
par(mar=c(4,5,2,1))
xPlace <- c(1, 1.2, 1.3, 1.5, 1.6)
plot(xPlace, tapply(x$MPE, x$metric, mean), "n", xlab="", ylab="", bty="l", yaxt="n", ylim=range(tapply(x$MPE, x$metric, quantile, c(0.25, 0.75))), xaxt="n", xlim=c(0.9, 1.7))
axis(side=1, at=xPlace, labels=FALSE)
u <- par('usr')
text(xPlace, rep(u[3]-(u[4]-u[3])*0.1, 3), pos=1, c("avgS", "S25", "S75", "Smsy", "Sgen1"), xpd=NA)
abline(h=0)
segments(x0=xPlace, x1=xPlace, y0=tapply(x$MPE, x$metric, quantile, 0.25), y1=tapply(x$MPE, x$metric, quantile, 0.75), col=grey(0.8), lwd=10, lend=2)
points(xPlace, tapply(x$MPE, x$metric, mean), pch=19, cex=0.8)
points(xPlace, tapply(x$MPE, x$metric, median), pch=4)
A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Percent error", cex=par('cex'))

legend("bottomleft", bty="n", pch=c(19, 4, 15), col=c(1,1,grey(0.8)), c("mean", "median", "inter-quartile range"))
###############################################################################
# Over increasing bias in observation
###############################################################################
obsBias <- log(1/seq(1, 5, 0.2))
out_obsBias <- runSensitivity(basePar = simPar, 
															sensitivity = list("obs_bias", obsBias), 
															nSim = 5000,
															nCores = 10)

time_obsBias <- out_obsBias[[2]]
out_obsBias <- out_obsBias[[1]]

# 13.69 mins on 7 cores for 21 values
# 23 mins min on 7 cores for 21 values w/ 6000 simulations each
# 20 mins on 10 cores for 21 values w/ 5000 simulations each (Feb 26)

plotSensitivity(out_obsBias, sensitivityPar = obsBias, baseValue = simPar$obs_bias, sensitivityName = "Bias in observed spawners", comparison = 1, include.all.ppnBias = FALSE, newQuartz = FALSE)

#-----
# How do things change if the real exp factor is 3?
obsBias2 <- obsBias[seq(1, 20, 2)]
simPar2 <- simPar
simPar2$ExpFactor3 <- 3

out_obsBias2 <- runSensitivity(basePar = simPar2, 
															 sensitivity = list("obs_bias", obsBias2), 
															 nSim = 5000,
															 nCores = 10)

plotSensitivity(out_obsBias2[[1]], sensitivityPar = obsBias2, baseValue = log(1/simPar2$ExpFactor3), sensitivityName = "Bias in observed spawners", comparison = 1, include.all.ppnBias = FALSE, newQuartz = FALSE)

plot(sensitivityPar, MPE$avgS[,'mean'], "o", pch=19, ylim=c(-0.3,3), xlab="Bias in observed spawners", bty="l", ylab="", yaxt="n")
A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Percent error", cex=par('cex'))
abline(v =  log(1/simPar2$ExpFactor3), lty=2)
abline(h=0)
lines(sensitivityPar, MPE$Sgen1[,'mean'], "o", col=statusCols['r'], pch=15, lwd=2)
lines(sensitivityPar, MPE$Smsy[,'mean'], "o", col=statusCols['g'], pch=16, lwd=2)
legend("topleft", pch=c(19, 15, 16), col=c(1, statusCols['r'], statusCols['g']), lwd=2, c("avgS", "Sgen1 (lower)", "Smsy (upper)"), bty="n")

#-----
# How do things change if the real productivity/harvest gives red status?
simPar3 <- simPar
simPar3$targetHarvest <- 0.8

out_obsBias3 <- runSensitivity(basePar = simPar3, 
															 sensitivity = list("obs_bias", obsBias2), 
															 nSim = 5000,
															 nCores = 10)

plotSensitivity(out_obsBias3[[1]], sensitivityPar = obsBias2, baseValue = log(1/simPar3$ExpFactor3), sensitivityName = "Bias in observed spawners", comparison = 1, include.all.ppnBias = FALSE, newQuartz = TRUE)

plot(sensitivityPar, MPE$avgS[,'mean'], "o", pch=19, ylim=c(-1, 0.5), xlab="Bias in observed spawners", bty="l", ylab="", yaxt="n")
A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Percent error", cex=par('cex'))
abline(v =  log(1/simPar2$ExpFactor3), lty=2)
abline(h=0)
lines(sensitivityPar, MPE$Sgen1[,'mean'], "o", col=statusCols['r'], pch=15, lwd=2)
lines(sensitivityPar, MPE$Smsy[,'mean'], "o", col=statusCols['g'], pch=16, lwd=2)
legend("topleft", pch=c(19, 15, 16), col=c(1, statusCols['r'], statusCols['g']), lwd=2, c("avgS", "Sgen1 (lower)", "Smsy (upper)"), bty="n")


# #----------------------------------------------------
# # MPE in benchmarks
# comparison <- 1 # Which set of 'observed' data to use. 1 = all, 2 = perfect obs, 3 = complete coverage
# 
# MPE_obsBias <- list(
# 	avgS =  matrix(NA, nrow = length(obsBias), ncol = 1), 
# 	SR = matrix(NA, nrow = length(obsBias), ncol = 2), 
# 	HS = matrix(NA, nrow = length(obsBias), ncol = 2))
# 
# for (i in 1:length(obsBias)) {
# 	MPE_obsBias[[1]][i, ] <- mean(out_obsBias[[i]][[comparison]]$avgS)
# 	MPE_obsBias[[2]][i, ] <- apply(out_obsBias[[i]][[comparison]]$SR[, 1:2], 2, mean)
# 	MPE_obsBias[[3]][i, ] <- apply(out_obsBias[[i]][[comparison]]$HS[, 1:2], 2, mean)
# }
# 
# par(mfcol=c(1,2), mar=c(4,5,2,2), oma=c(0,0,0,0), mgp=c(3,1,0))
# 
# # # MPE in avgS
# # plot(obsBias, MPE_obsBias$avgS, "l", xlab="Bias in observation error", ylab="", main="Avg. spawners", bty="l", yaxt="n", ylim=c(-0.8, 0.8))
# # A <- axis(side=2, labels = FALSE)
# # mtext(side=2, line=4, "Mean percent error")
# # axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 1, format="f"), "%", sep=""), las=1)
# # abline(v = log(1/1.5), lty=2)
# # abline(v = 0, lty=3)
# # abline(h=0)
# 
# # MPE in SR benchmarks
# plot(obsBias, MPE_obsBias$SR[, 1], "n", xlab="Bias in observation error", ylab="", main="SR benchmarks", bty="l", yaxt="n", col=2, ylim=c(-0.8, 0.8))
# A <- axis(side=2, labels = FALSE)
# mtext(side=2, line=4, "Mean percent error")
# axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
# lines(obsBias, MPE_obsBias$avgS, col=grey(0.6))
# lines(obsBias, MPE_obsBias$SR[, 2], lty=2, col=2)
# lines(obsBias, MPE_obsBias$SR[, 1], col=2)
# legend(-1.6, 0.8, col=c(2,2,grey(0.6)), lty=c(1,2,1), c("Lower (Sgen)", "Upper (Smsy)", "AvgS"), bty="n")
# abline(v = log(1/1.5), lty=2)
# abline(v = 0, lty=3)
# abline(h=0)
# 
# # MPE in HS benchmarks
# plot(obsBias, MPE_obsBias$HS[, 1], "n", xlab="Bias in observation error", ylab="", main="HS benchmarks", bty="l", yaxt="n", col=4, ylim=c(-0.8, 1.5))
# A <- axis(side=2, labels = FALSE)
# axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
# mtext(side=2, line=4, "Mean percent error")
# lines(obsBias, MPE_obsBias$avgS, col=grey(0.6))
# lines(obsBias, MPE_obsBias$HS[, 2], lty=2, col=4)
# lines(obsBias, MPE_obsBias$HS[, 1], col=4)
# legend(-1.6, 1.7, col=c(4, 4, grey(0.6)), lty=c(1,2,1), c("Lower (S25th)", "Upper (S75th)", "AvgS"), bty="n")
# abline(v = log(1/1.5), lty=2)
# abline(v = 0, lty=3)
# abline(h=0)
# 
# 
# #----------------------------------------------------
# # Proportion of simulations with wrong status
# 
# dummy <- matrix(NA, nrow = length(obsBias), ncol = 2, dimnames = list(NULL, c("SR", "HS")))
# ppnWrong_obsBias <- list(all = dummy, a = dummy, b = dummy)
# for (i in 1:length(obsBias)) {
# 	for(j in 1:3){
# 	ppnWrong_obsBias[[j]][i,1] <- sum(out_obsBias[[i]][[j]]$SR[, 3]>=4)/nSim
# 	ppnWrong_obsBias[[j]][i,2] <- sum(out_obsBias[[i]][[j]]$HS[, 3]>=4)/nSim
# }}
# 
# quartz(width = 6, height = 4, pointsize = 10)
# par(mfrow = c(1,2), mar=c(4,4,5,1), oma=c(0,0,0,0))
# for(i in 1:2){
# 	plot(obsBias, ppnWrong_obsBias[[1]][, i], "l", ylim = range(ppnWrong_obsBias), col = c(2,4)[i], lwd=1.5, las=1, bty="l", ylab = "Proportion of sim. with wrong status", xlab = "Bias in observation error")
# 	lines(obsBias, ppnWrong_obsBias[[2]][, i],  col=c(2,4)[i], lty=2)
# 	lines(obsBias, ppnWrong_obsBias[[3]][, i],  col=c(2,4)[i], lty=3)
# 	
# 	axis(side=3, at=obsBias[seq(1,21,5)], labels=seq(1, 5, 0.2)[seq(1,21,5)])
# 	mtext(side = 3, "'true' Expansion Factor III", line=2.5)
# 	mtext(side = 3, c("Stock-recruitment", "Historic spawners")[i], line=4, font=2)
# 	abline(v = log(1/1.5), lty=2)
# 	
# 	if(i==1) legend("topright", lty=c(3,2,1), col=2, c("Obs. bias", "Incomplete monitoring", "Both"))
# }


###############################################################################
# Over increasing bias in catch
###############################################################################
catchBias <- log(seq(0.5, 1.5, 0.1))
out_catchBias <- runSensitivity(basePar = simPar, 
															sensitivity = list("catch_bias", catchBias), 
															nSim = 5000,
															nCores = 10)

time_catchBias <- out_catchBias[[2]]
out_catchBias <- out_catchBias[[1]]

plotSensitivity(out_catchBias, sensitivityPar = catchBias, baseValue = simPar$catch_bias, sensitivityName = "Bias in observed catch", comparison = 1, include.all.ppnBias = FALSE, newQuartz = FALSE)

# 13 mins on 10 cores for 11 values w/ 5000 simulations each (Feb 26)
#----------------------------------------------------
# # MPE in benchmarks
# comparison <- 1 # Which set of 'observed' data to use. 1 = all, 2 = perfect obs, 3 = complete coverage
# 
# MPE_catchBias <- list(
# 	avgS =  matrix(NA, nrow = length(catchBias), ncol = 1), 
# 	SR = matrix(NA, nrow = length(catchBias), ncol = 2), 
# 	HS = matrix(NA, nrow = length(catchBias), ncol = 2))
# 
# for (i in 1:length(catchBias)) {
# 	MPE_catchBias[[1]][i, ] <- mean(out_catchBias[[i]][[comparison]]$avgS)
# 	MPE_catchBias[[2]][i, ] <- apply(out_catchBias[[i]][[comparison]]$SR[, 1:2], 2, mean)
# 	MPE_catchBias[[3]][i, ] <- apply(out_catchBias[[i]][[comparison]]$HS[, 1:2], 2, mean)
# }
# 
# par(mfcol=c(1,2), mar=c(4,5,2,2), oma=c(0,0,0,0), mgp=c(3,1,0))
# 
# # # MPE in avgS
# # plot(obsBias, MPE_obsBias$avgS, "l", xlab="Bias in observation error", ylab="", main="Avg. spawners", bty="l", yaxt="n", ylim=c(-0.8, 0.8))
# # A <- axis(side=2, labels = FALSE)
# # mtext(side=2, line=4, "Mean percent error")
# # axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 1, format="f"), "%", sep=""), las=1)
# # abline(v = log(1/1.5), lty=2)
# # abline(v = 0, lty=3)
# # abline(h=0)
# 
# # MPE in SR benchmarks
# plot(catchBias, MPE_catchBias$SR[, 1], "n", xlab="Bias in observed catch", ylab="", main="SR benchmarks", bty="l", yaxt="n", col=2, ylim=c(0, 0.3))
# A <- axis(side=2, labels = FALSE)
# mtext(side=2, line=4, "Mean percent error")
# axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
# lines(catchBias, MPE_catchBias$avgS, "o", col=grey(0.6), cex=0.5, pch=21, bg="white")
# lines(catchBias, MPE_catchBias$SR[, 2], "o", lty=2, col=2, cex=0.5, pch=21, bg="white")
# lines(catchBias, MPE_catchBias$SR[, 1], "o", col=2, cex=0.5, pch=19)
# legend(-1.6, 0.8, col=c(2,2,grey(0.6)), lty=c(1,2,1), c("Lower (Sgen)", "Upper (Smsy)", "AvgS"), bty="n")
# abline(v = 0, lty=3)
# abline(h=0)
# 
# # MPE in HS benchmarks
# plot(catchBias, MPE_catchBias$HS[, 1], "n", xlab="Bias in observed catch", ylab="", main="HS benchmarks", bty="l", yaxt="n", col=4, ylim=c(0.05, 1.7))
# A <- axis(side=2, labels = FALSE)
# axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
# mtext(side=2, line=4, "Mean percent error")
# lines(catchBias, MPE_catchBias$avgS, "o", col=grey(0.6), cex=0.5, pch=21, bg="white")
# lines(catchBias, MPE_catchBias$HS[, 2], "o", lty=2, col=4, cex=0.5, pch=21, bg="white")
# lines(catchBias, MPE_catchBias$HS[, 1], "o", col=4, cex=0.5, pch=19)
# legend(-1.6, 1.7, col=c(4, 4, grey(0.6)), lty=c(1,2,1), c("Lower (S25th)", "Upper (S75th)", "AvgS"), bty="n")
# abline(v = 0, lty=3)
# abline(h=0)
# 

# #----------------------------------------------------
# # Proportion of simulations with wrong status
# 
# dummy <- matrix(NA, nrow = length(catchBias), ncol = 2, dimnames = list(NULL, c("SR", "HS")))
# ppnWrong_catchBias <- list(all = dummy, a = dummy, b = dummy)
# for (i in 1:length(catchBias)) {
# 	for(j in 1:3){
# 		ppnWrong_catchBias[[j]][i,1] <- sum(out_catchBias[[i]][[j]]$SR[, 3]>=4)/nSim
# 		ppnWrong_catchBias[[j]][i,2] <- sum(out_catchBias[[i]][[j]]$HS[, 3]>=4)/nSim
# 	}}
# 
# quartz(width = 6, height = 4, pointsize = 10)
# par(mfrow = c(1,2), mar=c(4,4,5,1), oma=c(0,0,0,0))
# for(i in 1:2){
# 	plot(catchBias, ppnWrong_catchBias[[1]][, i], "l", ylim = range(ppnWrong_obsBias), col = c(2,4)[i], lwd=1.5, las=1, bty="l", ylab = "Proportion of sim. with wrong status", xlab = "Bias in observed catch")
# 	lines(catchBias, ppnWrong_catchBias[[2]][, i],  col=c(2,4)[i], lty=2)
# 	lines(catchBias, ppnWrong_catchBias[[3]][, i],  col=c(2,4)[i], lty=3)
# 	
# 	mtext(side = 3, c("Stock-recruitment", "Historic spawners")[i], line=4, font=2)
# 	abline(v = 0, lty=3)
# 	
# 	if(i==1) legend("topright", lty=c(3,2,1), col=2, c("Obs. bias", "Incomplete monitoring", "Both"))
# }


###############################################################################
# Over increasing correlation among subpopulations
###############################################################################
# correlPop <- c(-0.02, -0.01, 0, 0.01, 0.02, 0.04, seq(0.06, 0.98, 0.04), 0.99)# seq(-0.02, 0.99, 0.02)
correlPop <- c(-0.02, 0, 0.02, seq(0.05, 0.95, 0.05), 0.99)# seq(-0.02, 0.99, 0.02)

out_correlPop <- runSensitivity(basePar = simPar, 
															sensitivity = list("correlPop", correlPop), 
															nSim = 5000,
															nCores = 10)

time_correlPop <- out_correlPop[[2]]
out_correlPop <- out_correlPop[[1]]

# 28 mins for 6000 simulations at 23 levels of correlPop
# 
plotSensitivity(out_correlPop, sensitivityPar = correlPop, baseValue = simPar$correlPop, sensitivityName = "Correlation among subpopulations", comparison = 1, include.all.ppnBias = FALSE, newQuartz = TRUE)

###############################################################################
# Over increasing interannual variability in age-at-return
###############################################################################
ageErr <- seq(0.1, 0.9, 0.1)

out_ageErr <- runSensitivity(basePar = simPar, 
																sensitivity = list("ageErr", ageErr), 
																nSim = 5000,
																nCores = 9)

time_ageErr <- out_ageErr[[2]]
out_ageErr <- out_ageErr[[1]]

# 8 mins for 5000 simulations on 9 cores with at 9 levels of ageErr
# 
plotSensitivity(out_ageErr, sensitivityPar = ageErr, baseValue = simPar$ageErr, sensitivityName = "Interannual variability in age-at-maturity", comparison = 1, include.all.ppnBias = FALSE, newQuartz = TRUE)

###############################################################################
# Different numbers of subpopulations and proportion indicator
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

out_nPop <- runSensitivity(basePar = simPar_nPop, sensitivity = NULL, nSim = 5000, nCores = 5, multiPar = TRUE)

time_nPop <- out_nPop[[2]]
out_nPop <- out_nPop[[1]]

#------------------------

MPE_nPop <- list(
	avgS =  matrix(NA, nrow = length(simPar_nPop), ncol = 4, dimnames=list(NULL, c("mean", "median", "25th", "75th"))), 
	S25 = matrix(NA, nrow = length(simPar_nPop), ncol = 4, dimnames=list(NULL, c("mean", "median", "25th", "75th"))),
	S75 = matrix(NA, nrow = length(simPar_nPop), ncol = 4, dimnames=list(NULL, c("mean", "median", "25th", "75th"))),
	Sgen1 = matrix(NA, nrow = length(simPar_nPop), ncol = 4, dimnames=list(NULL, c("mean", "median", "25th", "75th"))), 
	Smsy = matrix(NA, nrow = length(simPar_nPop), ncol = 4, dimnames=list(NULL, c("mean", "median", "25th", "75th")))
)

for (i in 1:length(simPar_nPop)) {
	
	MPE_nPop$avgS[i, 'mean'] <- mean(out_nPop[[i]][[comparison]]$avgS)
	MPE_nPop$avgS[i, 'median'] <- quantile(out_nPop[[i]][[comparison]]$avgS, 0.5)
	MPE_nPop$avgS[i, '25th'] <- quantile(out_nPop[[i]][[comparison]]$avgS, 0.25)
	MPE_nPop$avgS[i, '75th'] <- quantile(out_nPop[[i]][[comparison]]$avgS, 0.75)
	
	MPE_nPop$S25[i, 'mean'] <- mean(out_nPop[[i]][[comparison]]$HS[,1])
	MPE_nPop$S25[i, 'median'] <- quantile(out_nPop[[i]][[comparison]]$HS[,1], 0.5)
	MPE_nPop$S25[i, '25th'] <- quantile(out_nPop[[i]][[comparison]]$HS[,1], 0.25)
	MPE_nPop$S25[i, '75th'] <- quantile(out_nPop[[i]][[comparison]]$HS[,1], 0.75)
	
	MPE_nPop$S75[i, 'mean'] <- mean(out_nPop[[i]][[comparison]]$HS[,2])
	MPE_nPop$S75[i, 'median'] <- quantile(out_nPop[[i]][[comparison]]$HS[,2], 0.5)
	MPE_nPop$S75[i, '25th'] <- quantile(out_nPop[[i]][[comparison]]$HS[,2], 0.25)
	MPE_nPop$S75[i, '75th'] <- quantile(out_nPop[[i]][[comparison]]$HS[,2], 0.75)
	
	MPE_nPop$Sgen1[i, 'mean'] <- mean(out_nPop[[i]][[comparison]]$SR[,1])
	MPE_nPop$Sgen1[i, 'median'] <- quantile(out_nPop[[i]][[comparison]]$SR[,1], 0.5)
	MPE_nPop$Sgen1[i, '25th'] <- quantile(out_nPop[[i]][[comparison]]$SR[,1], 0.25)
	MPE_nPop$Sgen1[i, '75th'] <- quantile(out_nPop[[i]][[comparison]]$SR[,1], 0.75)
	
	MPE_nPop$Smsy[i, 'mean'] <- mean(out_nPop[[i]][[comparison]]$SR[,2])
	MPE_nPop$Smsy[i, 'median'] <- quantile(out_nPop[[i]][[comparison]]$SR[,2], 0.5)
	MPE_nPop$Smsy[i, '25th'] <- quantile(out_nPop[[i]][[comparison]]$SR[,2], 0.25)
	MPE_nPop$Smsy[i, '75th'] <- quantile(out_nPop[[i]][[comparison]]$SR[,2], 0.75)
}

#----------------------------------------------------
# Proportion of simulations with wrong status

dummy <- matrix(NA, nrow = length(simPar_nPop), ncol = 2, dimnames = list(NULL, c("HS", "SR")))
ppnWrong_nPop <- list(all = dummy, a = dummy, b = dummy)
ppnRisky_nPop <- ppnWrong_nPop
for (i in 1:length(simPar_nPop)) {
	for(j in 1:3){
		ppnWrong_nPop[[j]][i,2] <- sum(out_nPop[[i]][[j]]$SR[, 3]>=4)/nSim
		ppnWrong_nPop[[j]][i,1] <- sum(out_nPop[[i]][[j]]$HS[, 3]>=4)/nSim
		
		ppnRisky_nPop[[j]][i,2] <- sum(out_nPop[[i]][[j]]$SR[, 3]>5 & out_nPop[[i]][[j]]$SR[, 3] != 7)/nSim
		ppnRisky_nPop[[j]][i,1] <- sum(out_nPop[[i]][[j]]$SR[, 3]>5 & out_nPop[[i]][[j]]$SR[, 3] != 7)/nSim
		
	}}

#----------------------------------------------------
# Plot
cols <- c(ind = "#475A83", nonInd = "#C2B642")
x <- c(1,2,2.5,3.5,4)
quartz(width = 7, height = 5, pointsize=10)

layout(matrix(c(6,2,4,1,2,4,1,3,5,10,3,5,9,7,8, 9,7,8), 6, 3, byrow=TRUE))
par(mar=c(4,5,2,2), oma=c(0,0,1,0), mgp=c(3,1,0))

plot(x, MPE_nPop$avgS[, 'mean'], "n", xlab="", ylab="", bty="l", yaxt="n", ylim=range(MPE_nPop$avgS), xaxt="n", xlim=c(0.5,4.5))
axis(side=1, at=x, labels=FALSE)
u <- par('usr')
text(c(1,2.25,3.75), rep(u[3]-(u[4]-u[3])*0.1, 3), pos=1, c("n=35", "n=10", "n=140"), xpd=NA)
abline(h=0)
segments(x0=x, x1=x, y0=MPE_nPop$avgS[, '25th'], y1=MPE_nPop$avgS[, '75th'], col=c(grey(0.8), rep(c("#C2B642", "#475A83"), 2)), lwd=3, lend=2)
points(x, MPE_nPop$avgS[, 'mean'], pch=19, cex=0.8)
points(x, MPE_nPop$avgS[, 'median'], pch=4)
mtext(side=3, "Average spawners", font=2, line=1)
A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Percent error", cex=par('cex'))

for(j in c(3,2,5,4)){
	plot(x, MPE_nPop[[j]][, 'mean'], "n", xlab="", ylab="", bty="l", yaxt="n", ylim=range(MPE_nPop[[j]]), xaxt="n", xlim=c(0.5,4.5))
	axis(side=1, at=x, labels=FALSE)
	u <- par('usr')
	text(c(1,2.25,3.75), rep(u[3]-(u[4]-u[3])*0.1, 3), pos=1, c("n=35", "n=10", "n=140"), xpd=NA)
	abline(h=0)
	segments(x0=x, x1=x, y0=MPE_nPop[[j]][, '25th'], y1=MPE_nPop[[j]][, '75th'], col=c(grey(0.8), rep(c("#C2B642", "#475A83"), 2)), lwd=3, lend=2)
	points(x, MPE_nPop[[j]][, 'mean'], pch=19, cex=0.8)
	points(x, MPE_nPop[[j]][, 'median'], pch=4, cex=0.8)
	
	A <- axis(side=2, labels = FALSE)
	axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
	mtext(side=2, line=4, "Percent error", cex=par('cex'))
	if(j == 3) mtext(side=3, "Historic spawners", font=2, line=1)
	if(j == 5){
		mtext(side=3, "Stock-recruitment", font=2, line=1)
		mtext(side=4, "upper", line=1)
	}
	if(j==4) mtext(side=4, "lower", line=1)
}

plot(1,1,"n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
legend("center", pch=c(19, 4, 15), col=c(1,1,grey(0.8)), c("mean", "median", "inter-quartile range"), xpd=NA, bty="n")

for(i in 1:2){
	ylims <- range(ppnWrong_nPop[[1]][, i])
	plot(x, ppnWrong_nPop[[1]][, i], "n", ylim = c(0, max(ylims)), las=1, bty="l", ylab = "", xlab = "", yaxs="i", xaxt="n", xlim=c(0.5,4.5))
	axis(side=1, at=x, labels=FALSE)
	u <- par('usr')
	text(c(1,2.25,3.75), rep(u[3]-(u[4]-u[3])*0.1, 3), pos=1, c("n=35", "n=10", "n=140"), xpd=NA)
	mtext(side=2, line=4, "Proportion incorrect", cex=par('cex'))
	segments(x0=x, x1=x, y0=0, y1=ppnWrong_nPop[[1]][, i], col=c(grey(0.8), rep(c("#C2B642", "#475A83"), 2)), lwd=10, lend=2)
	segments(x0=x, x1=x, y0=ppnWrong_nPop[[1]][, i], y1=ppnWrong_nPop[[1]][, i]-ppnRisky_nPop[[1]][, i], lwd=10, lend=2)
}


###############################################################################
# With no change in monitoring coverage
###############################################################################
# correlPop <- c(-0.02, -0.01, 0, 0.01, 0.02, 0.04, seq(0.06, 0.98, 0.04), 0.99)# seq(-0.02, 0.99, 0.02)
simPar_mon <- list(simPar, simPar)

# In second simulation, remove decline in coverage
simPar_mon[[2]]['ppnChange_ind'] <- 0
simPar_mon[[2]]['ppnChange_nonInd'] <- 0

out_mon <- runSensitivity(basePar = simPar_mon, sensitivity = NULL, nSim = 5000, nCores = 2, multiPar = TRUE)

time_mon <- out_mon[[2]]
out_mon <- out_mon[[1]]

# par(mfrow=c(1,2))
# plotStatusDiff(out_mon[[1]][[comparison]]$SR[1:nSimEnough,3])
# plotStatusDiff(out_mon[[1]][[comparison]]$HS[1:nSimEnough,3])
# 
# plotStatusDiff(out_mon[[2]][[comparison]]$SR[1:nSimEnough,3])
# plotStatusDiff(out_mon[[2]][[comparison]]$HS[1:nSimEnough,3])

comparison <- 1 # Which set of 'observed' data to use. 1 = all, 2 = perfect obs, 3 = complete coverage

MPE_mon <- list(
	avgS = list(mean = matrix(NA, nrow = 3, ncol = 2), 
						lower = matrix(NA, nrow = 3, ncol = 2),
						upper = matrix(NA, nrow = 3, ncol = 2)),
	SR = list(mean = matrix(NA, nrow = 6, ncol = 2), 
						lower = matrix(NA, nrow = 6, ncol = 2),
						upper = matrix(NA, nrow = 6, ncol = 2)),
	HS = list(mean = matrix(NA, nrow = 6, ncol = 2), 
						lower = matrix(NA, nrow = 6, ncol = 2),
						upper = matrix(NA, nrow = 6, ncol = 2)))

for (i in 1:2) {# for each monitoring scenario (decline & no decline)
	for(j in 1:3){ # for each comparison
		MPE_mon$avgS$mean[j,i] <- mean(out_mon[[i]][[j]]$avgS) 
		MPE_mon$SR$mean[j,i] <- mean(out_mon[[i]][[j]]$SR[, 1]) #lower
		MPE_mon$SR$mean[j+3,i] <- mean(out_mon[[i]][[j]]$SR[, 2]) #upper
		MPE_mon$HS$mean[j,i] <- mean(out_mon[[i]][[j]]$HS[, 1]) #lower
		MPE_mon$HS$mean[j+3,i] <- mean(out_mon[[i]][[j]]$HS[, 2]) #upper
		
		MPE_mon$avgS$lower[j,i] <- quantile(out_mon[[i]][[j]]$avgS, 0.25) 
		MPE_mon$SR$lower[j,i] <- quantile(out_mon[[i]][[j]]$SR[, 1], 0.25) #lower
		MPE_mon$SR$lower[j+3,i] <- quantile(out_mon[[i]][[j]]$SR[, 2], 0.25) #upper
		MPE_mon$HS$lower[j,i] <- quantile(out_mon[[i]][[j]]$HS[, 1], 0.25) #lower
		MPE_mon$HS$lower[j+3,i] <- quantile(out_mon[[i]][[j]]$HS[, 2], 0.25) #upper
		
		MPE_mon$avgS$upper[j,i] <- quantile(out_mon[[i]][[j]]$avgS, 0.75) 
		MPE_mon$SR$upper[j,i] <- quantile(out_mon[[i]][[j]]$SR[, 1], 0.75) #lower
		MPE_mon$SR$upper[j+3,i] <- quantile(out_mon[[i]][[j]]$SR[, 2], 0.75) #upper
		MPE_mon$HS$upper[j,i] <- quantile(out_mon[[i]][[j]]$HS[, 1], 0.75) #lower
		MPE_mon$HS$upper[j+3,i] <- quantile(out_mon[[i]][[j]]$HS[, 2], 0.75) #upper
	}}


MPE_mon.avgS <- data.frame(decl = rep(c("decl", "noDecl"), each=nSim*3), comp = rep(rep(c("both", "a", "b"), each=nSim), 2), MPE = c(out_mon[[1]][[1]]$avgS, out_mon[[1]][[2]]$avgS, out_mon[[1]][[3]]$avgS, out_mon[[2]][[1]]$avgS, out_mon[[2]][[2]]$avgS, out_mon[[2]][[3]]$avgS))
MPE_mon.avgS$comp <- relevel(MPE_mon.avgS$comp, ref="both")

MPE_mon.SR <- data.frame(decl = rep(rep(c("decl", "noDecl"), each=nSim*3), 2), comp = rep(rep(rep(c("both", "a", "b"), each=nSim), 2), 2), benchmark = rep(c("lower", "upper"), each = nSim*3*2), MPE = c(out_mon[[1]][[1]]$SR[,1], out_mon[[1]][[2]]$SR[,1], out_mon[[1]][[3]]$SR[,1], out_mon[[2]][[1]]$SR[,1], out_mon[[2]][[2]]$SR[,1], out_mon[[2]][[3]]$SR[,1], out_mon[[1]][[1]]$SR[,2], out_mon[[1]][[2]]$SR[,2], out_mon[[1]][[3]]$SR[,2], out_mon[[2]][[1]]$SR[,2], out_mon[[2]][[2]]$SR[,2], out_mon[[2]][[3]]$SR[,2]))
MPE_mon.SR$comp <- relevel(MPE_mon.SR$comp, ref="both")

												 	
MPE_mon.HS <- data.frame(decl = rep(rep(c("decl", "noDecl"), each=nSim*3), 2), comp = rep(rep(rep(c("both", "a", "b"), each=nSim), 2), 2), benchmark = rep(c("lower", "upper"), each = nSim*3*2), MPE = c(out_mon[[1]][[1]]$HS[,1], out_mon[[1]][[2]]$HS[,1], out_mon[[1]][[3]]$HS[,1], out_mon[[2]][[1]]$HS[,1], out_mon[[2]][[2]]$HS[,1], out_mon[[2]][[3]]$HS[,1], out_mon[[1]][[1]]$HS[,2], out_mon[[1]][[2]]$HS[,2], out_mon[[1]][[3]]$HS[,2], out_mon[[2]][[1]]$HS[,2], out_mon[[2]][[2]]$HS[,2], out_mon[[2]][[3]]$HS[,2]))
MPE_mon.HS$comp <- relevel(MPE_mon.HS$comp, ref="both")

par(mfrow=c(2,1), mar=c(1, 5, 2, 1), oma=c(8, 0, 0, 0))
bp <- boxplot(MPE ~ decl*comp*benchmark, data = MPE_mon.SR, las=2, col=c(2, "#FF000050"), names = NULL, xaxt="n", yaxt="n")
A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 1, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Mean percent error")
abline(h = 0)
abline(v = seq(2.5, 12, 2), lwd=0.5)
abline(v = 6.5, lwd=2)
mtext(side=3, line = 1, adj=0, "a) Stock-recruitment", font=2)
legend("topright", fill=c(2, "#FF000050"), c("Decline in coverage", "No decline in coverage"))

bp <- boxplot(MPE ~ decl*comp*benchmark, data = MPE_mon.HS, las=2, col=c(4, "#0000FF50"), xaxt="n", yaxt="n")
A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 1, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Mean percent error")
abline(h = 0)
abline(v = seq(2.5, 12, 2), lwd=0.5)
abline(v = 6.5, lwd=2)
text(seq(1.5, 12, 2), rep(-0.35, 6), rep(c("Both", "No bias", "Complete\ncoverage"), 2), xpd=NA)
mtext(side=3, line = 1, adj=0, "b) Historic spawners", font=2)

text(c(3.5, 9.5), rep(-0.55, 2), c("Lower benchmark", "Upper benchmark"), xpd=NA, font=2)
segments(x0=6.5, x1=6.5, y0=-0.2, y1=-0.6, lwd=2, xpd=NA)


dummy <- matrix(NA, nrow = 2, ncol = 2, dimnames = list(NULL, c("SR", "HS")))
ppnWrong_mon <- list(all = dummy, a = dummy, b = dummy)
for (i in 1:2) {
	for(j in 1:3){
		ppnWrong_mon[[j]][i,1] <- sum(out_mon[[i]][[j]]$SR[, 3]>=4)/nSim
		ppnWrong_mon[[j]][i,2] <- sum(out_mon[[i]][[j]]$HS[, 3]>=4)/nSim
	}}


#-------------
meanX <- c(
	tapply(MPE_mon.avgS$MPE[MPE_mon.avgS$comp == "both"], MPE_mon.avgS$decl[MPE_mon.avgS$comp == "both"], mean),
	tapply(MPE_mon.HS$MPE[MPE_mon.HS$comp == "both" & MPE_mon.HS$benchmark=="lower"], MPE_mon.HS$decl[MPE_mon.HS$comp == "both"& MPE_mon.HS$benchmark=="lower"], mean),
	tapply(MPE_mon.HS$MPE[MPE_mon.HS$comp == "both" & MPE_mon.HS$benchmark=="upper"], MPE_mon.HS$decl[MPE_mon.HS$comp == "both"& MPE_mon.HS$benchmark=="upper"], mean),
	tapply(MPE_mon.SR$MPE[MPE_mon.SR$comp == "both" & MPE_mon.SR$benchmark=="lower"], MPE_mon.SR$decl[MPE_mon.SR$comp == "both"& MPE_mon.SR$benchmark=="lower"], mean),
	tapply(MPE_mon.SR$MPE[MPE_mon.SR$comp == "both" & MPE_mon.SR$benchmark=="upper"], MPE_mon.SR$decl[MPE_mon.SR$comp == "both"& MPE_mon.SR$benchmark=="upper"], mean))
medianX <- c(
	tapply(MPE_mon.avgS$MPE[MPE_mon.avgS$comp == "both"], MPE_mon.avgS$decl[MPE_mon.avgS$comp == "both"], median),
	tapply(MPE_mon.HS$MPE[MPE_mon.HS$comp == "both" & MPE_mon.HS$benchmark=="lower"], MPE_mon.HS$decl[MPE_mon.HS$comp == "both"& MPE_mon.HS$benchmark=="lower"], median),
	tapply(MPE_mon.HS$MPE[MPE_mon.HS$comp == "both" & MPE_mon.HS$benchmark=="upper"], MPE_mon.HS$decl[MPE_mon.HS$comp == "both"& MPE_mon.HS$benchmark=="upper"], median),
	tapply(MPE_mon.SR$MPE[MPE_mon.SR$comp == "both" & MPE_mon.SR$benchmark=="lower"], MPE_mon.SR$decl[MPE_mon.SR$comp == "both"& MPE_mon.SR$benchmark=="lower"], median),
	tapply(MPE_mon.SR$MPE[MPE_mon.SR$comp == "both" & MPE_mon.SR$benchmark=="upper"], MPE_mon.SR$decl[MPE_mon.SR$comp == "both"& MPE_mon.SR$benchmark=="upper"], median))

lowerX <- c(
	tapply(MPE_mon.avgS$MPE[MPE_mon.avgS$comp == "both"], MPE_mon.avgS$decl[MPE_mon.avgS$comp == "both"], quantile, 0.25),
	tapply(MPE_mon.HS$MPE[MPE_mon.HS$comp == "both" & MPE_mon.HS$benchmark=="lower"], MPE_mon.HS$decl[MPE_mon.HS$comp == "both"& MPE_mon.HS$benchmark=="lower"], quantile, 0.25),
	tapply(MPE_mon.HS$MPE[MPE_mon.HS$comp == "both" & MPE_mon.HS$benchmark=="upper"], MPE_mon.HS$decl[MPE_mon.HS$comp == "both"& MPE_mon.HS$benchmark=="upper"], quantile, 0.25),
	tapply(MPE_mon.SR$MPE[MPE_mon.SR$comp == "both" & MPE_mon.SR$benchmark=="lower"], MPE_mon.SR$decl[MPE_mon.SR$comp == "both"& MPE_mon.SR$benchmark=="lower"], quantile, 0.25),
	tapply(MPE_mon.SR$MPE[MPE_mon.SR$comp == "both" & MPE_mon.SR$benchmark=="upper"], MPE_mon.SR$decl[MPE_mon.SR$comp == "both"& MPE_mon.SR$benchmark=="upper"], quantile, 0.25))


upperX <- c(
	tapply(MPE_mon.avgS$MPE[MPE_mon.avgS$comp == "both"], MPE_mon.avgS$decl[MPE_mon.avgS$comp == "both"], quantile, 0.75),
	tapply(MPE_mon.HS$MPE[MPE_mon.HS$comp == "both" & MPE_mon.HS$benchmark=="lower"], MPE_mon.HS$decl[MPE_mon.HS$comp == "both"& MPE_mon.HS$benchmark=="lower"], quantile, 0.75),
	tapply(MPE_mon.HS$MPE[MPE_mon.HS$comp == "both" & MPE_mon.HS$benchmark=="upper"], MPE_mon.HS$decl[MPE_mon.HS$comp == "both"& MPE_mon.HS$benchmark=="upper"], quantile, 0.75),
	tapply(MPE_mon.SR$MPE[MPE_mon.SR$comp == "both" & MPE_mon.SR$benchmark=="lower"], MPE_mon.SR$decl[MPE_mon.SR$comp == "both"& MPE_mon.SR$benchmark=="lower"], quantile, 0.75),
	tapply(MPE_mon.SR$MPE[MPE_mon.SR$comp == "both" & MPE_mon.SR$benchmark=="upper"], MPE_mon.SR$decl[MPE_mon.SR$comp == "both"& MPE_mon.SR$benchmark=="upper"], quantile, 0.75))

quartz(width = 6.3, height = 3.3, pointsize=10)
par(mar=c(4,5,2,1))
xPlace <- c(1, 1.02, 1.2, 1.22, 1.3, 1.32, 1.5,1.52, 1.6,1.62)

plot(xPlace, meanX, "n", xlab="", ylab="", bty="l", yaxt="n", ylim=range(lowerX, upperX), xaxt="n", xlim=c(0.9, 1.72))
axis(side=1, at=c(1.01, 1.21, 1.31, 1.51, 1.61), labels=FALSE)
u <- par('usr')
text(c(1.01, 1.21, 1.31, 1.51, 1.61), rep(u[3]-(u[4]-u[3])*0.1, 3), pos=1, c("avgS", "S25", "S75", "Smsy", "Sgen1"), xpd=NA)
abline(h=0)
segments(x0=xPlace, x1=xPlace, y0=lowerX, y1=upperX, col=rep(grey(c(0.8, 0.4)), 5), lwd=10, lend=2)
points(xPlace, meanX, pch=19, cex=0.8)
points(xPlace, medianX, pch=4)
A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Percent error", cex=par('cex'))

legend("topright", bty="n", pch=c(19, 4, 15), col=c(1,1,grey(0.8), grey(0.3)), c("mean", "median", "inter-quartile range"))
