library(here)
library(lme4)

library(matrixcalc)
library(corpcor)

###############################################################################
# Load data
###############################################################################

CU_names <- read.csv("data/CentralCoastChum_CUs.csv")
CU_names <- CU_names[,1]

datSR <- read.csv("data/PSE_chumRecruitsPerSpawner.csv")
datSR <- subset(datSR, is.element(datSR$location, CU_names))

# What are the central coast chum CU names?
CUs <- unique(datSR$location)

###############################################################################
# Create dataframe with x = spawners and y = log(recruits/spawner)
###############################################################################
datFit <- data.frame(CU = rep(NA, 480), year = rep(NA, 480), x = rep(NA, 480), y = rep(NA, 480))
I <- 1

for(i in 1:length(CUs)){
	
	datSR.i <- subset(datSR, datSR$location == CUs[i])
	yr <- unique(datSR.i$year)

	datFit$CU[I:(I + length(yr) - 1)] <- rep(CUs[i], length(yr))
	datFit$year[I:(I + length(yr) - 1)] <- yr

	for(j in 1:length(yr)){
		if(length(datSR.i$datavalue[which(datSR.i$parameter == "Spawners" & datSR.i$year == yr[j])]) > 0){
			datFit$x[c(I:(I + length(yr) - 1))[j]] <- datSR.i$datavalue[which(datSR.i$parameter == "Spawners" & datSR.i$year == yr[j])]
			if(length(datSR.i$datavalue[which(datSR.i$parameter == "Recruits" & datSR.i$year == yr[j])] > 0)){
				datFit$y[c(I:(I + length(yr) - 1))[j]] <- log(datSR.i$datavalue[which(datSR.i$parameter == "Recruits" & datSR.i$year == yr[j])]/datFit$x[c(I:(I + length(yr) - 1))[j]])
			}
		}
	}

	I <- (I + length(yr))	
}
datFit <- datFit[is.na(datFit$x)==FALSE,]
hist(datFit$x * 10^-5)

datFit$x_scaled <- datFit$x * 10^-5

###############################################################################
# Fit SR models
###############################################################################

fit <- lmer(y ~ x_scaled + (1 + x_scaled| CU), data = datFit)

summary(fit)

###############################################################################
# Productivity
###############################################################################

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
abline(v = a_Dorner08, col=2, lty=2)

###############################################################################
# Residuals
###############################################################################
eps <- resid(fit)
hist(eps)

plot(1:10, eps[1:10], "n", xlim = range(datFit$year, na.rm=TRUE), ylim = range(eps), xlab = "Year", ylab = "Residual", las=1, bty="l")
abline(h = 0)
for(i in 1:8){
	lines(datFit$year[datFit$CU == as.numeric(CUs)[i]], eps[datFit$CU == as.numeric(CUs)[i]], "o", col=i)
}

eps_scaled <- eps / max(abs(eps))

# dataframe with residuals by year by CU
datRes <- data.frame(year = 1954:2012, CU1 = rep(NA, 59), CU2= rep(NA, 59), CU3 = rep(NA, 59), CU4 = rep(NA, 59), CU5= rep(NA, 59), CU6 = rep(NA, 59), CU7 = rep(NA, 59), CU8 = rep(NA, 59))
datResRaw <- datRes
for(i in 1:8){
	I <- which(datFit$CU == as.numeric(CUs)[i])
	datRes[is.element(datRes$year, datFit$year[I]),i+1] <- eps_scaled[I] 
	datResRaw[is.element(datResRaw$year, datFit$year[I]),i+1] <- eps[I] 
}

# Autocorrelation in residuals across subpopulations
eps1 <- c(as.matrix(datRes[1:58,2:9]))
eps2 <- c(as.matrix(datRes[2:59,2:9]))
rho <- cor(eps1[is.na(eps1*eps2) == FALSE], eps2[is.na(eps1*eps2) == FALSE])

###############################################################################
# Plotting
###############################################################################

# par(mar=c(4,8,1,1))
# plot(1,1,"n", xlim = range(datFit$year, na.rm=TRUE), yaxt="n", ylim=c(0,16), xlab = "Year", ylab = "", bty="n")
# abline(v = seq(1955, 2012, 5), col=grey(0.8), lty=3)
# abline(h = c(1:8)*2 -1, col=grey(0.8), lty=3)
# for(i in 1:8){
# 	polygon(
# 		x = c(datFit$year[datFit$CU == as.numeric(CUs)[i]], max(datFit$year[datFit$CU == as.numeric(CUs)[i]]), min(datFit$year[datFit$CU == as.numeric(CUs)[i]])), 
# 		y = (i*2 - 1) + c(eps[datFit$CU == as.numeric(CUs)[i]], 0, 0),
# 		border = NA,
# 		col = "#FF000030")
# }
# 
# text(rep(1952, 8), c(1:8)*2 -1, CUs, adj=1, xpd=NA)

par(mar=c(4,10,1,1))
plot(1,1,"n", xlim = range(datFit$year, na.rm=TRUE), yaxt="n", ylim=c(0,16), xlab = "Year", ylab = "", bty="n", main= "Scaled residuals")
abline(v = seq(1955, 2012, 5), col=grey(0.8), lty=3)
abline(h = c(1:8)*2 -1, col=grey(0.8), lty=3)
for(i in 1:8){
	points(c(1954:2012), rep((i*2 - 1), 59), cex = 3*abs(datRes[,i+1]), col=c(2,4)[c(datRes[,i+1]>=0)+1], bg = c("#FF000040", "#0000FF40")[c(datRes[,i+1]>=0)+1], pch=21)
}
text(rep(1952, 8), c(1:8)*2 -1, CUs, adj=1, xpd=NA)

###############################################################################
# Calculate correlation among stocks
###############################################################################

resMat <- as.matrix(datRes[,2:9])

corrMat <- matrix(NA, nrow =8, ncol = 8)
for(i in 1:8){
	for(j in 1:8){
		corrMat[i,j] <- cor(resMat[is.na(resMat[,i]+resMat[,j]) == FALSE,i], resMat[is.na(resMat[,i]+resMat[,j]) == FALSE,j])
	}
}

par(mar=c(2,12,2,1))
plot(0:9, 0:9, "n", xlim=c(0.5, 8.5), ylim=c(0.5, 8.5), xaxt="n", yaxt="n", xlab="", ylab="", bty="n", main="Correlation in residuals between CUs")
for(i in 1:8){
	for(j in 1:8){
		points(i, 9-j, cex=5*abs(corrMat[i,j]), col=c(2,4)[as.numeric(corrMat[i,j]>=0)+1], bg = c("#FF000040", "#0000FF40")[c(corrMat[i,j]>=0)+1], pch=21)
	}
}
text(rep(0, 8), c(8:1), paste(CUs, " (", c(1:8), ")", sep=""), adj=1, xpd=NA)
text(1:8, rep(0, 8), paste("(", c(1:8), ")", sep=""), xpd=NA)

# Mean correlation
correlPop <- mean(corrMat[lower.tri(corrMat)])

# Variance
sigma_u <- summary(fit)$sigma
sd(eps)

# Is within pop different?
sdPop <- apply(as.matrix(datResRaw[,2:9]), 2, sd, na.rm=TRUE)
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

###############################################################################
# Density-dependence parameter
###############################################################################

# b = 1/Smax - Smax is unfished spawner abundance though...
# What is the range of spawner abundances?
# THe CU-level data really tell us nothing about this, since spawner abundances
# are aggregated at the level of the CU. Need to look at the nuSEDS data.

# THe following data are from the NCC Salmon Database from LGL
# Downloaded from the GitHub folder nccdbv2/r-package/code-dev/output/ncc-streams/
# production -> test_out.xlsx on Febraury 8, 2019

nccEsc <- read.csv("data/nccdbv2_NCCStreamEscapement.csv")
nccEsc <- subset(nccEsc, nccEsc$SPP == "CM" & is.element(nccEsc$CU_fname, c("CM::Douglas-Gardner", "CM::Bella Coola River-Late", "CM::Bella Coola-Dean Rivers", "CM::Hecate Lowlands", "CM::Mussel-Kynoch","CM::Rivers Inlet", "CM::Smith Inlet", "CM::Spiller-Fitz Hugh-Burke", "CM::Wannock")))

# Look at escapement since 1954, when we started to get catch data
nccEsc2 <- t(as.matrix(nccEsc[,90:152]))

hist(log(nccEsc2), xlab="log(escapement)", yaxt="n", main="")
axis(side=2, at=c(0, 1000, 2000), las=1)
hist(log(nccEsc2[,which(nccEsc$IsIndicator == "Y")]), col="#FF000040", border=NA, add=TRUE)
hist(log(nccEsc2[,which(nccEsc$IsIndicator == "N")]), col="#0000FF40", border=NA, add=TRUE)
legend("topright", fill=c("#FF000040", "#0000FF40"), c("indicator", "non-indicator"))
abline(v = mean(nccEsc2[,which(nccEsc$IsIndicator == "Y")]))

#------------------------------------------------------------------------------
# Max escapement per stream
maxEsc <- apply(nccEsc2, 2, max, na.rm=TRUE)
hist(log(maxEsc), main="", xlab=expression(log(max(S[obs]))), las=1)
hist(log(maxEsc[which(nccEsc$Indicator == "Y")]), col="#FF000040", border=NA, add=TRUE)
hist(log(maxEsc[which(nccEsc$Indicator == "N")]), col="#0000FF40", border=NA, add=TRUE)
abline(v = mean(log(maxEsc[which(nccEsc$Indicator == "Y")])), col=2, lwd=2)
abline(v = mean(log(maxEsc[which(nccEsc$Indicator == "N")]), na.rm=TRUE), col=4, lwd=2)

# Max spawner abundance for indicator and non-indicator
avgSmax <- c(mean(log(maxEsc[which(nccEsc$Indicator == "Y")])), mean(log(maxEsc[which(nccEsc$IsIndicator == "N")]), na.rm=TRUE))
sdSmax <- c(sd(log(maxEsc[which(nccEsc$Indicator == "Y")])), sd(log(maxEsc[which(nccEsc$IsIndicator == "N")]), na.rm=TRUE))

x <- seq(0, 14, 0.1)
lines(x, dnorm(x, avgSmax[1], sdSmax[1])*length(maxEsc[which(nccEsc$Indicator == "Y")]), col=2, lty=2)
lines(x, dnorm(x, avgSmax[2], sdSmax[2])*sum(is.na(maxEsc[which(nccEsc$IsIndicator == "N")])==FALSE), col=4, lty=2)
legend("topleft", fill=c(2,4), c("indicator", "non-indicator"), bty="n")

text(9.5, 80, expression(mu == 9.1), col=2, adj=0)
text(9.5, 70, expression(sigma == 1.38), col=2, adj=0)
text(5.5, 80, expression(mu == 7.2), col=4, adj=0)
text(5.5, 70, expression(sigma == 2.07), col=4, adj=0)

sum(is.na(maxEsc[which(nccEsc$Indicator == "Y")])==FALSE)
sum(is.na(maxEsc[which(nccEsc$Indicator == "N")])==FALSE)

# 
Smax <- rlnorm(200, avgSmax[1], sdSmax[1])
mean(log(Smax))
sd(log(Smax))

hist(log(Smax), col="#00FF0060", border=NA, add=TRUE)

#------------------------------------------------------------------------------
# Mean escapement per stream
meanEsc <- apply(nccEsc2, 2, mean, na.rm=TRUE)
hist(log(meanEsc), main="", xlab="log(Smax)")
hist(log(meanEsc[which(nccEsc$Indicator == "Y")]), col="#FF000040", border=NA, add=TRUE)
hist(log(meanEsc[which(nccEsc$Indicator == "N")]), col="#0000FF40", border=NA, add=TRUE)
abline(v = mean(log(meanEsc[which(nccEsc$Indicator == "Y")])), col=2, lwd=2)
abline(v = mean(log(meanEsc[which(nccEsc$Indicator == "N")]), na.rm=TRUE), col=4, lwd=2)

# Max spawner abundance for indicator and non-indicator
avgSmean <- c(mean(log(meanEsc[which(nccEsc$Indicator == "Y")])), mean(log(maxEsc[which(nccEsc$IsIndicator == "N")]), na.rm=TRUE))
sdSmean <- c(sd(log(meanEsc[which(nccEsc$Indicator == "Y")])), sd(log(maxEsc[which(nccEsc$IsIndicator == "N")]), na.rm=TRUE))

x <- seq(0, 14, 0.1)
lines(x, dnorm(x, avgSmean[1], sdSmean[1])*100, col=2, lty=2)
lines(x, dnorm(x, avgSmean[2], sdSmean[2])*100, col=4, lty=2)
legend("topleft", fill=c(2,4), c("indicator", "non-indicator"), bty="n")

sum(is.na(maxEsc[which(nccEsc$Indicator == "Y")])==FALSE)
sum(is.na(maxEsc[which(nccEsc$Indicator == "N")])==FALSE)

# 
Smax <- rlnorm(200, avgSmax[1], sdSmax[1])
mean(log(Smax))
sd(log(Smax))

hist(log(Smax), col="#00FF0060", border=NA, add=TRUE)
