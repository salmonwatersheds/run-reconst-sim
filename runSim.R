library(mvtnorm)
library(here)

source("populationSubmodFns.R")
source("obsSubmodFns.R")
source("expansionFactors.R")


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
