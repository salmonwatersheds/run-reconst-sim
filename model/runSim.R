library(mvtnorm)
library(here)
library(data.table) # for shift function
library(gsl) # for lambert_W0 function

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
# cuCustomCorrMat <- read.csv(here("data/baseCorrMatrix.csv"), stringsAsFactors=F)

set.seed(987)
a <- rnorm(simPar$nPop, simPar$a_mean, simPar$sigma_a)



plot(1:nYears, spawners[, 1], "n", ylim=range(spawners))
for(i in 1:nPop) lines(1:nYears, spawners[, i], col="#00000020")


i<-3
plot(1:nYears, spawners[, i], "o", pch=19)
points(1:nYears, obsSpawners[, i], col=2)
abline(h=a[i]/simPar$b, lty=3)


plot(1:nYears, trueCatch, "o", pch=19)
points(1:nYears, obsCatch, col=2)

# # Which is faster, rbinom or sample?
# n <- 10^6
# system.time(sample(c(0,1), size = n, replace = TRUE, prob = c(0.4, 0.6)))
# system.time(rbinom(n, size = 1, prob = 0.6))
# # Sample is clearly faster

# What is the long run avg ppn sampled for indicator and non-indicator?
ppnSampled <- cbind(apply(sampled[,1:simPar$nIndicator], 1, sum)/simPar$nIndicator, 
			apply(sampled[,(simPar$nIndicator+1):nPop], 1, sum)/simPar$nNonIndicator)

apply(ppnSampled, 2, mean)
simPar$propSampled

# How do diffferent expansion factors compare to true spawners?

plot(1:simYears, apply(spawners[(simPar$gen + 3):nYears, ], 1, sum), "o", ylab = "Spawners", xlab = "Years", las = 1)
points(1:simYears, spawnersExp1, col = 2, pch = 19, cex = 0.8)

lines(1:simYears, apply(spawners[(simPar$gen + 3):nYears, 1:simPar$nIndicator], 1, sum), col=2)
points(1:simYears, spawnersExp2, pch = 19, cex = 0.8)

lines(1:simYears, spawnersExp3, col=3)


plot(1:simYears, spawnersExp3, "o")
abline(h=quantile(spawnersExp3, c(0.25, 0.75)), col=c(2,3))
points(simYears, spawnersExp3[simYears], col=2, pch=19, cex=0.5)
abline(h=c(Sgen1, Smsy), lty=2, col=c(2,3))

# Use whole function
nSim <- 1000
SRBias <- matrix(NA, nrow = nSim, ncol = 2)
percBias <- matrix(NA, nrow = nSim, ncol = 2)
SRStatus <- numeric(nSim)
percStatus <- numeric(nSim)

for(i in 1:nSim){
	out <- reconstrSim(simPar, a)
	SRBias[i, ] <- out$performance$benchBias['SR', ]
	percBias[i, ] <- out$performance$benchBias['perc', ]
	SRStatus[i] <- out$performance$statusDiff[1]
	percStatus[i] <- out$performance$statusDiff[2]
}

par(mfrow=c(1,2))
plotStatusDiff(SRStatus); mtext(side=3, "Stock-recruit metric")
plotStatusDiff(percStatus); mtext(side=3, "Percentile metric")

length(which(SRStatus>=4))/nSim
length(which(percStatus>=4))/nSim

apply(SRBias, 2, mean)
apply(percBias, 2, mean)
apply(abs(SRBias), 2, mean)
apply(abs(percBias), 2, mean)


correlPop <- -0.1
sigMat <- matrix(as.numeric(simPar$sigma_u), nrow = 1, ncol = nPop) 
covMat <- t(sigMat) %*% sigMat # Calculate shared variance
corMat <- covMat * correlPop # Correct based on correlation
diag(corMat) <- as.numeric(simPar$sigma_u^2) # Add variance
rmvnorm(1, rep(0, nPop), sigma = corMat)
