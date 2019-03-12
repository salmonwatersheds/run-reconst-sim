library(here)
library(lme4)

library(matrixcalc)
library(corpcor)

###############################################################################
# Load data
###############################################################################

CU_names <- read.csv("data/CentralCoastChum_CUs.csv")

datSR <- read.csv("data/NCC_chum_streams_SR_data.csv")
datSR <- subset(datSR, is.element(datSR$CU, CU_names[1:9,'CU_findex']))

# What are the central coast chum CU names?
CU_names <- as.character(CU_names[1:9,'CU_name'])

###############################################################################
# Prune dataset
###############################################################################

datSR$y <- log(datSR$Recruits/datSR$Spawners)

datSR <- subset(datSR, is.na(datSR$y)==FALSE) # Goes from 23664 to 6104 observations

# Remove levels of River no longer in dataset
datSR$River <- as.factor(as.character(datSR$River)) 

# How many streams have at least 10 observations?
length(unique(datSR$River)) # total number of rivers = 260
length(which(tapply(datSR$y, datSR$River, length) >=10)) # total number of rivers with >=10 obs = 181
# Probably worth restricting to those with > 10 obs...
keepRivers <- names(which(tapply(datSR$y, datSR$River, length) >=10))
datSR <- subset(datSR, is.element(datSR$River, keepRivers))

###############################################################################
# Fit SR model
###############################################################################


datSR$Spawners_scaled <- datSR$Spawners * 10^-5

fit <- lmer(y ~ River + Spawners_scaled:River - 1 + (1|CU), data = datSR)

summary(fit)

###############################################################################
# Productivity
###############################################################################
keepRivers
CU.rivers <- character(length(keepRivers))

datLoc <- data.frame(name = keepRivers, CU = character(length(keepRivers)), xLONG = rep(NA, length(keepRivers)), yLAT = rep(NA, length(keepRivers)), indicator = rep(NA, length(keepRivers)))
datLoc$CU <- as.character(datLoc$CU)
datLoc$indicator <- as.character(datLoc$indicator)
for(i in 1:length(keepRivers)){
	datLoc$CU[i] <- as.character(datSR$CU[which(datSR$River==keepRivers[i])[1]])
	datLoc$xLONG[i] <- datSR$xLONG[which(datSR$River==keepRivers[i])[1]]
	datLoc$yLAT[i] <- datSR$yLAT[which(datSR$River==keepRivers[i])[1]]
	datLoc$indicator[i] <- as.character(datSR$Indicator[which(datSR$River==keepRivers[i])[1]])
	
}
	
datLoc$prod <- as.numeric(summary(fit)$coefficients[1:length(keepRivers),'Estimate'])
datLoc$densDep <- -as.numeric(summary(fit)$coefficients[(length(keepRivers)+1):(2*length(keepRivers)),'Estimate']) * 10^-5

# Histogram of productivity and density dependence
hist(datLoc$prod, breaks = seq(0, 3, 0.1), main="River-level productivity", xlab="Productivity (a)", las=1)
# for(i in 1:length(unique(CU.rivers))){
# 	hist(datLoc$prod[which(CU.rivers==unique(CU.rivers)[i])], col = "#FF000030", breaks=seq(0, 3, 0.1), main=unique(CU.rivers)[i])
# }
abline(v = mean(datLoc$prod), col=2, lty=2, lwd=2)

mean(datLoc$prod)
sd(datLoc$prod) # Standard deviation is about half of among CU

quantile(datLoc$prod, c(0.025, 0.975))

plot(datLoc$xLONG, datLoc$yLAT, pch=c(1,19)[as.numeric(datLoc$indicator=="N") + 1], col=c(2, 4)[(prod < mean(prod))+1], cex=prod/mean(prod))

plot(datLoc$xLONG, datLoc$yLAT, pch=c(1,19)[as.numeric(datLoc$indicator=="N") + 1], col=c(2, 4)[(densDep < mean(densDep))+1], cex=densDep/mean(densDep))


a_RE <- ranef(fit)$CU[,1]
a_mean <- summary(fit)$coefficients[1,1]
sigma_a <- 0.11865#summary(fit)$varcor
hist(a_mean + a_RE, xlim = a_mean + c(-3, 4.5)*sigma_a, col=grey(0.8), border="white", main="", las=1, freq=FALSE)
lines(seq(1, 2, 0.01), dnorm(seq(1, 2, 0.01), a_mean, sigma_a))
abline(v = a_mean, lty=2, lwd=1.5)

# Distribution assumed by Holt et al. (2018)
# 0.5 to 2.0, mean of 1.0

# Values from Dorner et al. (2008) for central coast chum
a_Dorner08 <- c(1.6, 1.32, 1.94)
abline(v = a_Dorner08, col=4, lty=2)

# Density dependence 

hist(log(1/datLoc$densDep), main="", xlab="log Smax (log(1/b))", las=1, , breaks=seq(3, 13, 0.5))
hist(log(1/datLoc$densDep[datLoc$indicator=="Y"]), col="#FF000030", border="#FF000080", add=TRUE, breaks=seq(3, 13, 0.5))
hist(log(1/datLoc$densDep[datLoc$indicator=="N"]), col="#0000FF30", border="#0000FF80", add=TRUE, breaks=seq(3, 13, 0.5))
abline(v = c(mean(log(1/datLoc$densDep[datLoc$indicator=="Y"])), mean(log(1/datLoc$densDep[datLoc$indicator=="N"]))), lwd=2, col=c(2,4), lty=2)

# Max spawner abundance for indicator and non-indicator
avgSmax <- c(mean(log(1/datLoc$densDep[datLoc$indicator=="Y"])), mean(log(1/datLoc$densDep[datLoc$indicator=="N"])))
sdSmax <- c(sd(log(1/datLoc$densDep[datLoc$indicator=="Y"])), sd(log(1/datLoc$densDep[datLoc$indicator=="N"])))

x <- seq(0, 14, 0.1)
lines(x, dnorm(x, avgSmax[1], sdSmax[1])*length(which(datLoc$indicator=="Y")), col=2, lty=2)
lines(x, dnorm(x, avgSmax[2], sdSmax[2])*length(which(datLoc$indicator=="N")), col=4, lty=2)
legend("topleft", fill=c(2,4), c("indicator", "non-indicator"), bty="n")

text(9.5, 25, expression(mu == 7.95), col=2, adj=0)
text(9.5, 22, expression(sigma == 1.44), col=2, adj=0)
text(4.5, 25, expression(mu == 6.95), col=4, adj=0)
text(4.5, 22, expression(sigma == 1.18), col=4, adj=0)


c(sd(log(1/datLoc$densDep[datLoc$indicator=="Y"])), sd(log(1/datLoc$densDep[datLoc$indicator=="N"])))

###############################################################################
# Residuals
###############################################################################
eps <- resid(fit)
hist(eps)

plot(1:10, eps[1:10], "n", xlim = range(datSR$BroodYear, na.rm=TRUE), ylim = range(eps), xlab = "Brood year", ylab = "Residual", las=1, bty="l")
abline(h = 0)
for(i in 1:length(keepRivers)){
	lines(datSR$BroodYear[datSR$River == keepRivers[i]], eps[datSR$River == keepRivers[i]], "o", col="#00000030")
}

eps_scaled <- eps / max(abs(eps))

# dataframe with residuals by year by River
datRes <- matrix(NA, nrow= length(1954:2010), ncol = length(keepRivers) +1)
datRes[,1] <- c(1954:2010)
datResRaw <- datRes
for(i in 1:length(keepRivers)){
	I <- which(datSR$River == keepRivers[i])
	datRes[is.element(datRes[, 1], datSR$BroodYear[I]), i+1] <- eps_scaled[I] 
	datResRaw[is.element(datRes[, 1], datSR$BroodYear[I]), i+1] <- eps[I] 
}

# Autocorrelation in residuals across subpopulations
eps1 <- c(as.matrix(datRes[1:(dim(datRes)[1] - 1), 2:(length(keepRivers) + 1)]))
eps2 <- c(as.matrix(datRes[2:(dim(datRes)[1]), 2:(length(keepRivers) + 1)]))
rho <- cor(eps1[is.na(eps1*eps2) == FALSE], eps2[is.na(eps1*eps2) == FALSE])


par(mar=c(4,10,1,1))
plot(1,1,"n", xlim = range(datSR$BroodYear, na.rm=TRUE), yaxt="n", ylim=c(0,2*length(keepRivers)), xlab = "Brood year", ylab = "", bty="n", main= "Scaled residuals")
abline(v = seq(1954, 2010, 5), col=grey(0.8), lty=3)
for(i in 1:length(keepRivers)){
	points(c(1954:2010), rep((i*2 - 1), length(c(1954:2010))), cex = 3*abs(datRes[,i+1]), col=c(2,4)[c(datRes[,i+1]>=0)+1], bg = c("#FF000040", "#0000FF40")[c(datRes[,i+1]>=0)+1], pch=21)
}

###############################################################################
# Calculate correlation among stocks
###############################################################################

resMat <- as.matrix(datRes[,2:(length(keepRivers)+1)])

# Order the corMat by latitude
o <- order(datLoc$CU, datLoc$yLAT)
corrMat <- matrix(NA, nrow =length(keepRivers), ncol = length(keepRivers))
for(i in 1:length(keepRivers)){
	for(j in 1:length(keepRivers)){
		corrMat[i,j] <- cor(resMat[is.na(resMat[,o[i]]+resMat[,o[j]]) == FALSE, o[i]], resMat[is.na(resMat[,o[i]]+resMat[,o[j]]) == FALSE, o[j]])
	}
}

par(mar=c(2,2,2,1))
plot(0:(length(keepRivers)+1), 0:(length(keepRivers)+1), "n", xlim=c(0.5, ((length(keepRivers)+0.5))), ylim=c(0.5, (length(keepRivers)+0.5)), xaxt="n", yaxt="n", xlab="", ylab="", bty="n", main="Correlation in residuals between CUs")
for(i in 1:length(keepRivers)){
	for(j in 1:length(keepRivers)){
		if(i > j){
			points(i, (length(keepRivers)+1)-j, cex=2*abs(corrMat[i,j]), col=c(2,4)[as.numeric(corrMat[i,j]>=0)+1], bg = c("#FF000040", "#0000FF40")[c(corrMat[i,j]>=0)+1], pch=21)
	}}
}
for(i in 2:length(unique(datLoc$CU))){
	I <- which(datLoc$CU[o] == sort(unique(datLoc$CU))[i])[1]
	Im <- which(datLoc$CU[o] == sort(unique(datLoc$CU))[i-1])[1]
	segments(
		x0 = I - 0.5, 
		x1 = length(keepRivers), 
		y0 = (length(keepRivers) + 1) - (I - 0.5), 
		y1 = (length(keepRivers) + 1) - (I - 0.5), 
		lwd=1.5)
	
	ty <- (length(keepRivers) + 1) -  ((I - 0.5) - (I-Im)/2)
	text(((I - 0.5) - (I-Im)/2), ty, sort(unique(datLoc$CU))[i-1], xpd=NA, pos=2)
}
text(162, 18, sort(unique(datLoc$CU))[length(unique(datLoc$CU))], xpd=NA, pos=2)

legend("bottomleft", pt.cex = 2*abs(corrMat[i,j]), bg="white")

# Mean correlation
correlPop <- mean(corrMat[lower.tri(corrMat)], na.rm=TRUE)

correlPop.withinCU <- numeric(length(unique(datLoc$CU)))
for(i in 1:length(unique(datLoc$CU))){
	if(length(which(datLoc$CU[o]==unique(datLoc$CU)[i])) > 1){
		ind <- which(datLoc$CU[o]==unique(datLoc$CU)[i])
		correlPop.withinCU[i] <- mean(lower.tri(corrMat[ind, ind]), na.rm=TRUE)
}}
names(correlPop.withinCU) <- unique(datLoc$CU)

mean(correlPop.withinCU[1:7])

# Variance
sigma_u <- summary(fit)$sigma
sd(eps)

# Is within pop different?
sdPop <- apply(datResRaw[, 2:(length(keepRivers)+1)], 2, sd, na.rm=TRUE)
mean(sdPop) # Slightly lower 

#-------------
correlPop <- 2*correlPop
nPop <- 35
sigMat <- matrix(as.numeric(sigma_u), nrow = 1, ncol = nPop) 
varMat <- t(sigMat) %*% sigMat # Calculate shared variance
covMat <- correlPop * varMat # Correct based on correlation
diag(covMat) <- as.numeric(sigma_u^2) # Add variance

is.positive.definite(covMat)

covMat2 <- make.positive.definite(covMat)
is.positive.definite(covMat2)
isSymmetric(covMat2)


rmvnorm(1, rep(0, nPop), sigma = covMat2)
