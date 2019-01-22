library(here)
library(RColorBrewer)

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
# What was the decline in monitoring?
###############################################################################

# Pre-decline and during decline
period1 <- which(datMon$year > (2013-50) & datMon$year <= 1986)
period2 <- which(datMon$year > 1986)

preDecl <- c(mean(datMon$ppn_indicator[period1]), mean(datMon$ppn_nonindicator[period1]), mean(datMon$ppn_total[period1]))

fitInd <- lm((datMon$ppn_indicator[period2] - preDecl[1]) ~ 0 + c(1:length(period2)))
fitNonInd <- lm((datMon$ppn_nonindicator[period2] - preDecl[2]) ~ 0 + c(1:length(period2)))
fit <- lm((datMon$ppn_total[period2] - preDecl[3]) ~ 0 + c(1:length(period2)))

###############################################################################
# Plotting
###############################################################################

par(mfrow=c(1,2), mar=c(3,4,2,1))
# Reproduce plots from English (2016)
plot(datMon$year, datMon$num_total, "l", col=2, lwd=1.5, las=1, bty="l", xlab="", ylab="Number", ylim=c(0, 1800))
lines(datMon$year, tot, lty=3, col=2)
lines(datMon$year, totIndicator, lty=3, col=4)
abline(h = seq(200, 1800, 200), col="#00000030")
abline(v = 1986, lty=3)
lines(datMon$year, datMon$num_indicator, col=4, lwd=1.5)
legend("topleft", col=c(2,4), lwd=1.5, c("Total", "Indicator"))
mtext(side=3, adj = 0, "a) Number of stream-species monitored", line=1)


plot(datMon$year, datMon$ppn_indicator, "l", col=4, lwd=1.5, las=1, bty="l", xlab="", ylab="Proportion", ylim=c(0, 1))
abline(h = seq(0, 1, 0.1), col="#00000030")
lines(datMon$year, datMon$ppn_nonindicator, col=3)
abline(v = 1986, lty=3)
legend("bottomleft", col=c(3,4), lwd=c(1, 1.5), lty=c(2,1), c("Non-Indicator", "Indicator"))
mtext(side=3, adj = 0, "b) Proportion of stream-species monitored", line=1)

# Decline lines
lines(datMon$year[period1], rep(preDecl[1], length(period1)), col=4, lwd=2, lty=2)
lines(datMon$year[period2], preDecl[1] + predict(fitInd), col=4, lwd=2, lty=2)

lines(datMon$year[period1], rep(preDecl[2], length(period1)), col=3, lwd=2, lty=2)
lines(datMon$year[period2], preDecl[2] + predict(fitNonInd), col=3, lwd=2, lty=2)

# What is the total ppn change in indicator and non-indicator?
tail(predict(fitInd), 1)
tail(predict(fitNonInd), 1)
tail(predict(fit), 1)
