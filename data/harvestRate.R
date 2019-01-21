library(here)
library(RColorBrewer)
###############################################################################
# Load data
###############################################################################

CU_names <- read.csv("data/CentralCoastChum_CUs.csv")
CU_names <- as.character(CU_names[1:8,1]) #Ignore Wannock (no data)

datER <- read.csv("data/PSE_chumExploitationRate.csv")
datER <- subset(datER, is.element(datER$location, CU_names))

# Are all CUs in ER data?
CU_names[is.element(CU_names, datER$location) == FALSE]

###############################################################################
# Summarize data
###############################################################################
meanER <- mean(datER$datavalue)
sdER <- sd(datER$datavalue)

meanER_period1 <- mean(datER$datavalue[datER$year <= 1990])
meanER_period2 <- mean(datER$datavalue[datER$year > 1990])

# Perhaps makes more sense to look at sd among CUs within each year?
mean(tapply(datER$datavalue, datER$year, sd))
# This is slightly lower than looking across all time.

###############################################################################
# Plot data
###############################################################################

par(mfrow=c(4,2), mar=c(2,2,1,1), oma=c(3,3,0,0))
for(i in 1:8){
	datER.i <- subset(datER, datER$location == CU_names[i])
	plot(datER.i$year, datER.i$datavalue, "o", xlim=c(1954, 2014), ylim=range(datER$datavalue, na.rm=TRUE), bty="n", xlab="", ylab="", las=1, main=CU_names[i])
}

# All together
colCU <- rainbow(8)#brewer.pal(8, "Accent")
par(mfrow=c(1,1), mar=c(4,4,2,1), oma=rep(0, 4))
plot(datER.i$year, datER.i$datavalue, "n", xlim=c(1954, 2014), ylim=range(datER$datavalue, na.rm=TRUE), bty="l", xlab="", ylab="", las=1)
for(i in 1:8){
	datER.i <- subset(datER, datER$location == CU_names[i])
	lines(datER.i$year, datER.i$datavalue, "o", col=colCU[i])
}
segments(x0=1954, x1=1990, y0=meanER_period1, y1=meanER_period1, lty=2, l)
segments(x0=1991, x1=2014, y0=meanER_period2, y1=meanER_period2, lty=2)
legend("topright", col=colCU, lwd=1, CU_names, ncol=2)
