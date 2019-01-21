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
nPop <- 10
sigMat <- matrix(as.numeric(sigma_u), nrow = 1, ncol = nPop) 
varMat <- t(sigMat) %*% sigMat # Calculate shared variance
covMat <- correlPop * varMat # Correct based on correlation
diag(covMat) <- as.numeric(sigma_u^2) # Add variance

is.positive.definite(covMat)

covMat2 <- make.positive.definite(covMat)
is.positive.definite(covMat2)
isSymmetric(covMat2)


rmvnorm(1, rep(0, nPop), sigma = covMat2)
