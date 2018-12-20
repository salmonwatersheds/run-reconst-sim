library(mvtnorm)
library(here)

source("populationSubmodFns.R")


plot(1:nYears, spawners[, 1], "n", ylim=range(spawners))
for(i in 1:nPop) lines(1:nYears, spawners[, i], col="#00000020")
