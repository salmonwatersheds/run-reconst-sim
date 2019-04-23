library(here)
library(lme4)

library(matrixcalc)
library(corpcor)

###############################################################################
# Setup and model fitting
###############################################################################

CU_names <- read.csv("data/CentralCoastChum_CUs.csv")

datSR <- read.csv("data/NCC_chum_streams_SR_data.csv")
datSR <- subset(datSR, is.element(datSR$CU, CU_names[1:9,'CU_findex']))

# What are the central coast chum CU names?
CU_names <- as.character(CU_names[1:9,'CU_name'])

#------------------------------------------------------------------------------
# Prune dataset
#------------------------------------------------------------------------------

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
nR <- length(keepRivers) # Number of rivers in final dataset

# Create dataframe containing river-level information
datLoc <- data.frame(name = keepRivers, CU = character(nR), xLONG = rep(NA, nR), yLAT = rep(NA, nR), indicator = rep(NA, nR))
datLoc$CU <- as.character(datLoc$CU)
datLoc$indicator <- as.character(datLoc$indicator)
for(i in 1:nR){
	datLoc$CU[i] <- as.character(datSR$CU[which(datSR$River==keepRivers[i])[1]])
	datLoc$xLONG[i] <- datSR$xLONG[which(datSR$River==keepRivers[i])[1]]
	datLoc$yLAT[i] <- datSR$yLAT[which(datSR$River==keepRivers[i])[1]]
	datLoc$indicator[i] <- as.character(datSR$Indicator[which(datSR$River==keepRivers[i])[1]])
}

#------------------------------------------------------------------------------
# Fit SR model
#------------------------------------------------------------------------------

# Scale spawners to avoid error
datSR$Spawners_scaled <- datSR$Spawners * 10^-5

# Fit model with random effect on productivity by CU and river-level
# estimates of productivity and density dependence
fit <- lmer(y ~ River + Spawners_scaled:River - 1 + (1|CU), data = datSR)

summary(fit)


# Include productivity and density dependence estimates in river-level data frame
datLoc$prod <- as.numeric(summary(fit)$coefficients[1:nR,'Estimate'])
datLoc$densDep <- -as.numeric(summary(fit)$coefficients[(nR+1):(2*nR),'Estimate']) * 10^-5

###############################################################################
# Productivity
###############################################################################
quartz(width = 6.3, height = 2.8, pointsize=10)
par(mfrow=c(1,2), mar=c(4,4,2,1))
# Histogram of productivity over all rivers
hist(datLoc$prod, breaks = seq(0, 3, 0.1), main="", xlab=expression(paste("Productivity (", italic(a), ")", sep="")), las=1, col=grey(0.8), border="white")
abline(v = mean(datLoc$prod), col='dodgerblue', lwd=2)
mtext(side=3, line=0.5, adj=0, "a) River-level productivity parameters")
# # Add indicator and non-indicator
# hist(datLoc$prod[datLoc$indicator=="Y"], col="#FF000030", border="#FF000080", add=TRUE, breaks=seq(0, 3, 0.1))
# hist(datLoc$prod[datLoc$indicator=="N"], col="#0000FF30", border="#0000FF80", add=TRUE, breaks=seq(0, 3, 0.1))
# Is there a significant difference in productivity between indicator and non-indicator?
t.test(datLoc$prod ~ datLoc$indicator)
# No (p = 0.7359)

mean(datLoc$prod)
sd(datLoc$prod) # Standard deviation is about half of among CU

# What range of mean productivities should be explored in base case?
quantile(datLoc$prod, c(0.1, 0.9))

# Distribution assumed by Holt et al. (2018)
# 0.5 to 2.0, mean of 1.0

# Values from Dorner et al. (2008) for central coast chum
a_Dorner08 <- c(1.6, 1.32, 1.94)
abline(v = a_Dorner08, col=4, lty=2)

###############################################################################
#  Density dependence 
###############################################################################

# Histogram of b is hard to interpret because it's so skewed:
hist(datLoc$densDep, breaks = seq(0,  0.025, 0.0001), xlim=c(0, 0.008), main="", xlab="Density dependence (b)")

# Histogram of Smax = 1/b
hist(1/datLoc$densDep, breaks = seq(0,  2e5, 1000), xlim=c(0, 0.5e5), main="", xlab="Smax (1/b)")
hist(1/datLoc$densDep[datLoc$indicator=="Y"], col="#FF000030", border="#FF000080", add=TRUE, breaks = seq(0,  2e5, 1000))
hist(1/datLoc$densDep[datLoc$indicator=="N"], col="#0000FF30", border="#0000FF80", add=TRUE, breaks = seq(0,  2e5, 1000))

# Histogram of log Smax = log(1/b)
hist(log(1/datLoc$densDep), main="", xlab="log Smax (log(1/b))", las=1, , breaks=seq(3, 13, 0.5))
hist(log(1/datLoc$densDep[datLoc$indicator=="Y"]), col="#FF000030", border="#FF000080", add=TRUE, breaks=seq(3, 13, 0.5))
hist(log(1/datLoc$densDep[datLoc$indicator=="N"]), col="#0000FF30", border="#0000FF80", add=TRUE, breaks=seq(3, 13, 0.5))
abline(v = c(mean(log(1/datLoc$densDep[datLoc$indicator=="Y"])), mean(log(1/datLoc$densDep[datLoc$indicator=="N"]))), lwd=2, col=c(2,4), lty=2)

# log Smax for indicator and non-indicator
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

###############################################################################
# Temporal autocorrelation in residuals
###############################################################################

# Extract all residuals from model fit
eps <- resid(fit)
hist(eps)

# Plot residuals over time for each river
plot(1:10, eps[1:10], "n", xlim = range(datSR$BroodYear, na.rm=TRUE), ylim = range(eps), xlab = "Brood year", ylab = "Residual", las=1, bty="l")
abline(h = 0)
for(i in 1:nR){
	lines(datSR$BroodYear[datSR$River == keepRivers[i]], eps[datSR$River == keepRivers[i]], "o", col="#00000030")
}
# Looks like increasing synchrony, but may just be a decline in monitoring!

# Scale residuals for plotting purposes
eps_scaled <- eps / max(abs(eps))

# Create dataframe with residuals by year by River
datRes <- matrix(NA, nrow= length(1954:2010), ncol = nR +1)
datRes[,1] <- c(1954:2010)
datResRaw <- datRes
for(i in 1:nR){
	I <- which(datSR$River == keepRivers[i])
	datRes[is.element(datRes[, 1], datSR$BroodYear[I]), i+1] <- eps_scaled[I] 
	datResRaw[is.element(datRes[, 1], datSR$BroodYear[I]), i+1] <- eps[I] 
}

# Autocorrelation in residuals across all rivers
eps1 <- c(as.matrix(datRes[1:(dim(datRes)[1] - 1), 2:(nR + 1)]))
eps2 <- c(as.matrix(datRes[2:(dim(datRes)[1]), 2:(nR + 1)]))
rho <- cor(eps1[is.na(eps1*eps2) == FALSE], eps2[is.na(eps1*eps2) == FALSE])

# Plot autocorrelation over time for each river (hard to interpret with so many rivers!)
par(mar=c(4,10,1,1))
plot(1,1,"n", xlim = range(datSR$BroodYear, na.rm=TRUE), yaxt="n", ylim=c(0,2*nR), xlab = "Brood year", ylab = "", bty="n", main= "Scaled residuals")
abline(v = seq(1954, 2010, 5), col=grey(0.8), lty=3)
for(i in 1:nR){
	points(c(1954:2010), rep((i*2 - 1), length(c(1954:2010))), cex = 3*abs(datRes[,i+1]), col=c(2,4)[c(datRes[,i+1]>=0)+1], bg = c("#FF000040", "#0000FF40")[c(datRes[,i+1]>=0)+1], pch=21)
}

###############################################################################
# Correlation among subpopulations in residuals
###############################################################################

# Create matrix of river-by-river residuals
resMat <- as.matrix(datRes[,2:(nR+1)])

# Order the corMat by CU
o <- order(datLoc$CU, datLoc$yLAT)
corrMat <- matrix(NA, nrow =nR, ncol = nR)
for(i in 1:nR){
	for(j in 1:nR){
		corrMat[i,j] <- cor(resMat[is.na(resMat[,o[i]]+resMat[,o[j]]) == FALSE, o[i]], resMat[is.na(resMat[,o[i]]+resMat[,o[j]]) == FALSE, o[j]])
	}
}

# Plot corMat (large plot)
par(mar=c(2,2,2,1))
plot(0:(nR+1), 0:(nR+1), "n", xlim=c(0.5, ((nR+0.5))), ylim=c(0.5, (nR+0.5)), xaxt="n", yaxt="n", xlab="", ylab="", bty="n", main="Correlation in residuals between CUs")
for(i in 1:nR){
	for(j in 1:nR){
		if(i > j){
			points(i, (nR+1)-j, cex=2*abs(corrMat[i,j]), col=c(2,4)[as.numeric(corrMat[i,j]>=0)+1], bg = c("#FF000040", "#0000FF40")[c(corrMat[i,j]>=0)+1], pch=21)
	}}
}

for(i in 2:length(unique(datLoc$CU))){
	I <- which(datLoc$CU[o] == sort(unique(datLoc$CU))[i])[1]
	Im <- which(datLoc$CU[o] == sort(unique(datLoc$CU))[i-1])[1]
	segments(
		x0 = I - 0.5, 
		x1 = nR, 
		y0 = (nR + 1) - (I - 0.5), 
		y1 = (nR + 1) - (I - 0.5), 
		lwd=1.5)
	
	ty <- (nR + 1) -  ((I - 0.5) - (I-Im)/2)
	text(((I - 0.5) - (I-Im)/2), ty, sort(unique(datLoc$CU))[i-1], xpd=NA, pos=2)
}
text(162, 18, sort(unique(datLoc$CU))[length(unique(datLoc$CU))], xpd=NA, pos=2)

legend("bottomleft", pt.cex = 2*abs(corrMat[i,j]), bg="white")

# Calculate the mean correlation

# a) Across all CUs
correlPop <- mean(corrMat[lower.tri(corrMat)], na.rm=TRUE)

# b) Within each CU
# *Only one subpopulations in CM-14 so ignored
correlPop.withinCU <- numeric(length(unique(datLoc$CU)))
for(i in 1:length(unique(datLoc$CU))){
	if(length(which(datLoc$CU[o]==unique(datLoc$CU)[i])) > 1){
		ind <- which(datLoc$CU[o]==unique(datLoc$CU)[i])
		correlPop.withinCU[i] <- mean(lower.tri(corrMat[ind, ind]), na.rm=TRUE)
}}
names(correlPop.withinCU) <- unique(datLoc$CU)

mean(correlPop.withinCU[1:7])

###############################################################################
# Variance in residuals within river
###############################################################################

# Variance overall
summary(fit)$sigma
sd(eps, na.rm=TRUE) # Different because of random effect?

# Within subpopulation
sdPop <- apply(datResRaw[, 2:(nR+1)], 2, sd, na.rm=TRUE)
sigma_u <- mean(sdPop) # Slightly lower 

# ###############################################################################
# # Test range in autocorrelation where variance-covariance is positive definite
# ###############################################################################
# nPop <- 35
# correlPop.all <- seq(-0.1, 0.99, 0.01)
# isPosDef <- numeric(correlPop.all)
# 
# for(i in 1:length(correlPop.all)){
# 	correlPop <- 2*correlPop
# 
# sigMat <- matrix(as.numeric(sigma_u), nrow = 1, ncol = nPop) 
# varMat <- t(sigMat) %*% sigMat # Calculate shared variance
# covMat <- correlPop * varMat # Correct based on correlation
# diag(covMat) <- as.numeric(sigma_u^2) # Add variance
# 
# is.positive.definite(covMat)
# 
# 
# rmvnorm(1, rep(0, nPop), sigma = covMat2)
