library(here)
library(RColorBrewer)

# indCol <- c(ind = "#FF0000", nonInd = "#0000FF")
indCol <- c(ind = "#000000", nonInd = "#666666")

###############################################################################
# Load data
###############################################################################

datMon <- read.csv("data/English2016_monitoringCoverage.csv")
datMon$num_nonindicator <- datMon$num_total - datMon$num_indicator

totIndicator <- rep(c(679, 642), 32)
tot <- 961 + totIndicator

datMon$ppn_nonindicator <- (datMon$num_nonindicator)/(tot-totIndicator)

datMon$ppn_total <- datMon$num_total/tot

###############################################################################
# 1) What was the decline in monitoring since the mid-1980s?
# This is the base case scenario (scenario B in paper)
###############################################################################

# Pre-decline and during decline
period1 <- which(datMon$year > (2013-50) & datMon$year <= 1986)
period2 <- which(datMon$year > 1986)

preDecl <- c(mean(datMon$ppn_indicator[period1]), mean(datMon$ppn_nonindicator[period1]), mean(datMon$ppn_total[period1]))

fitInd <- lm((datMon$ppn_indicator[period2] - preDecl[1]) ~ 0 + c(1:length(period2)))
fitNonInd <- lm((datMon$ppn_nonindicator[period2] - preDecl[2]) ~ 0 + c(1:length(period2)))
fit <- lm((datMon$ppn_total[period2] - preDecl[3]) ~ 0 + c(1:length(period2)))

#------------------------------------------------------------------------------
# Plotting
#------------------------------------------------------------------------------

# par(mfrow=c(1,2), mar=c(3,4,2,1))
# # Reproduce plots from English (2016)
# plot(datMon$year, datMon$num_total, "l", las=1, bty="l", xlab="", ylab="Number", ylim=c(0, 1800))
# # lines(datMon$year, tot, lty=3)
# # lines(datMon$year, totIndicator, lty=3, col=indCol['ind'])
# # abline(h = seq(200, 1800, 200), col="#00000030")
# abline(v = 1986, lty=2)
# lines(datMon$year, datMon$num_indicator, col=indCol['ind'])
# legend("topleft", col=c(1, indCol['ind']), c("Total", "Indicator"), bty="n", lwd=1)
# mtext(side=3, adj = 0, "a) Number of stream-species monitored", line=1)


# plot(datMon$year, datMon$ppn_indicator, "l", col=paste(indCol['ind'], 50, sep=""), las=1, bty="l", xlab="", ylab="Proportion", ylim=c(0, 1))
# # abline(h = seq(0, 1, 0.1), col="#00000030")
# lines(datMon$year, datMon$ppn_nonindicator, col=paste(indCol['nonInd'], 50, sep=""))
# abline(v = 1986, lty=2)
# legend("bottomleft", col=indCol, lwd=c(1, 1), c("Indicator", "Non-indicator"), bty="n")
# mtext(side=3, adj = 0, "b) Proportion of stream-species monitored", line=1)
# 
# # Decline lines
# lines(datMon$year[period1], rep(preDecl[1], length(period1)), col=indCol['ind'], lty=2)
# lines(datMon$year[period2], preDecl[1] + predict(fitInd), col=indCol['ind'], lty=2)
# 
# lines(datMon$year[period1], rep(preDecl[2], length(period1)), col=indCol['nonInd'], lty=2)
# lines(datMon$year[period2], preDecl[2] + predict(fitNonInd), col=indCol['nonInd'], lty=2)
# 
# # What is the total ppn change in indicator and non-indicator?
# tail(predict(fitInd), 1)
# tail(predict(fitNonInd), 1)
# tail(predict(fit), 1)

###############################################################################
# 2) What was the decline in monitoring in the last 5 years of data?
# Scenario C in paper
###############################################################################

# Pre-decline and during decline
period1.b <- which(datMon$year > (2013 - 50) & datMon$year <= 2013 - 5)
period2.b <- which(datMon$year > 2013 - 5)

preDecl.b <- c(ind = mean(datMon$ppn_indicator[period1.b]), ninInd = mean(datMon$ppn_nonindicator[period1.b]), tot = mean(datMon$ppn_total[period1.b]))

fitInd.b <- lm((datMon$ppn_indicator[period2.b] - preDecl.b[1]) ~ 0 + c(1:length(period2.b)))
fitNonInd.b <- lm((datMon$ppn_nonindicator[period2.b] - preDecl.b[2]) ~ 0 + c(1:length(period2.b)))
fit.b <- lm((datMon$ppn_total[period2.b] - preDecl.b[3]) ~ 0 + c(1:length(period2.b)))

# # What is the total ppn change in indicator and non-indicator?
# tail(predict(fitInd.b), 1)
# tail(predict(fitNonInd.b), 1)
# tail(predict(fit.b), 1)

# # Decline lines
# lines(datMon$year[c(period1.b, period2.b)], c(rep(preDecl.b[1], length(period1.b)), preDecl.b[1] + predict(fitInd.b)), col=indCol['ind'])
# lines(datMon$year[c(period1.b, period2.b)], c(rep(preDecl.b[2], length(period1.b)), preDecl.b[2] + predict(fitNonInd.b)), col=indCol['nonInd'])

###############################################################################
# 3) Look at decline in monitoring for chum salmon only using river-level SR data
# Scenario D in paper
###############################################################################

datSR <- read.csv("data/NCC_chum_streams_SR_data.csv")

# Create dummy variable that is 1 if there is a spawners estimate
datSR$monitored <- rep(0, dim(datSR)[1])
datSR$monitored[is.na(datSR$Spawners) == FALSE] <- 1


numMonitored <- cbind(
	all = tapply(datSR$monitored, datSR$BroodYear, sum),
	ind = tapply(datSR$monitored[datSR$Indicator=="Y"], datSR$BroodYear[datSR$Indicator=="Y"], sum),
	nonInd = tapply(datSR$monitored[datSR$Indicator=="N"], datSR$BroodYear[datSR$Indicator=="N"], sum))

ppnMonitored <- numMonitored / cbind(
	all = tapply(datSR$monitored, datSR$BroodYear, length),
	ind = tapply(datSR$monitored[datSR$Indicator=="Y"], datSR$BroodYear[datSR$Indicator=="Y"], length),
	nonInd = tapply(datSR$monitored[datSR$Indicator=="N"], datSR$BroodYear[datSR$Indicator=="N"], length))

year <- as.numeric(rownames(numMonitored))

# plot(year, numMonitored[,1], "o", ylim=c(0, 420))
# abline(v = 1986.5, col=2)
# points(year, numMonitored[,2], "o", col=2)
# points(year, numMonitored[,3], "o", col=4)

# # As a proportion of total
# plot(year, ppnMonitored[,1], "n", ylim=c(0, 1))
# abline(v = 1986.5, col=2)
# points(year, ppnMonitored[,2],  col=2, cex=0.8)
# points(year, ppnMonitored[,3],  col=4, cex=0.8)

tot <- cbind(all = tapply(datSR$monitored, datSR$BroodYear, length), ind = tapply(datSR$monitored[datSR$Indicator=="Y"], datSR$BroodYear[datSR$Indicator=="Y"], length), nonInd = tapply(datSR$monitored[datSR$Indicator=="N"], datSR$BroodYear[datSR$Indicator=="N"], length))[1,]

#------------------------------------------------------------------------------
# What was the decline in monitoring FOR CHUM?
#------------------------------------------------------------------------------

# Pre-decline and during decline
period1.c <- which(year > (2013-50) & year <= 1986)
period2.c <- which(year > 1986)

preDecl.c <- c(tot = mean(ppnMonitored[period1.c, 1]), ind = mean(ppnMonitored[period1.c, 2]), nonInd = mean(ppnMonitored[period1.c, 3]))

fitInd.c <- lm((ppnMonitored[period2.c, 2] - preDecl.c[2]) ~ 0 + c(1:length(period2.c)))
fitNonInd.c <- lm((ppnMonitored[period2.c, 3] - preDecl.c[3]) ~ 0 + c(1:length(period2.c)))
fit.c <- lm((ppnMonitored[period2.c, 1] - preDecl.c[1]) ~ 0 + c(1:length(period2.c)))

# Decline lines
lines(year[c(period1.c, period2.c)], c(rep(preDecl.c[2], length(period1.c)), preDecl.c[2] + predict(fitInd.c)), col=indCol['ind'], lwd=2, lty=2)
lines(year[c(period1.c, period2.c)], c(rep(preDecl.c[3], length(period1.c)), preDecl.c[3] + predict(fitNonInd.c)), col=indCol['nonInd'], lwd=2, lty=2)

# # What is the total ppn change in indicator and non-indicator?
# tail(predict(fitInd), 1)
# tail(predict(fitNonInd), 1)
# tail(predict(fit), 1)

# Declines in chum salmon monitoring have actually been more extreme than overall...
