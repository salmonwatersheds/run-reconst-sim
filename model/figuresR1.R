library(gplots)
statusCols <- c(g = "#8EB687", a = "#DFD98D", r = "#B66A64")

###############################################################################
# Base case simulations
###############################################################################

# out_base <- readRDS("workspaces/base100Mon_delisted_baseGreen.rds")

all3 <- TRUE

if(all3 == FALSE){
	quartz(width = 6.3, height= 3, pointsize = 10)
	I <- 1
	layout(matrix(c(2,1,1,2,1,1,7,3,3,5,4,4,5,4,4,8,6,6), nrow=3))
} else {
	quartz(width = 6, height= 7, pointsize = 9)
	par(oma=c(0,3.5,2,0))
	I <- 1:3
	layout(matrix(c(rep(c(2,1,1,8,7,7,14,13,13), 2), 19,3,3,21,9,9,23,15,15, rep(c(5,4,4,11,10,10,17,16,16),2), 20,6,6,22,12,12,24,18,18), nrow=9))
}

for(i in I){
for(k in 1:2){
	# Base case results
	statusDiff <- out_base[[i]][[k+1]][, 3]
	# statusDiff <- out_base[[i]]$HS[, 3]
	statusDiffMat <- matrix(NA, nrow = 3, ncol = 3, dimnames=list(c("R", "A", "G"), c("R", "A", "G")))
	statusDiffMat[1, 1] <- length(which(statusDiff == 3)) #/ nSim
	statusDiffMat[2, 2] <- length(which(statusDiff == 2)) #/ nSim
	statusDiffMat[3, 3] <- length(which(statusDiff == 1)) #/ nSim
	statusDiffMat[1, 2] <- length(which(statusDiff == 4)) #/ nSim
	statusDiffMat[1, 3] <- length(which(statusDiff == 5)) #/ nSim
	statusDiffMat[2, 1] <- length(which(statusDiff == 6))
	statusDiffMat[2, 3] <- length(which(statusDiff == 7))
	statusDiffMat[3, 1] <- length(which(statusDiff == 8))
	statusDiffMat[3, 2] <- length(which(statusDiff == 9))
	
	# Plot percentages in each category
	par(mar=c(5,5,0,0))
	plot(1:4, 1:4, "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
	for(s in 1:2){
		axis(side = s, at = c(1:4), labels=FALSE)
	}
	axis(side = 1, at = seq(1.5, 3.5, 1), c("Red", "Amber", "Green"), tck = 0, las=1)
	axis(side = 2, at = rev(seq(1.5, 3.5, 1)), c("Red", "Amber", "Green"), tck = 0, las=1)
	
	mtext(side = 2, line = 4, "Estimated status", cex=0.8)
	mtext(side = 1, line = 3.5, "True status", cex=0.8)
	for(h in 1:3){
		for(j in 1:3){
			if(j < h){ # Optimistic misclassifications
				polygon(x = c(j, j+1, j+1, j), y = 4 - c(h , h , h-1, h-1), border = NA, col = 1)
				fnt <- 2 
				fntcol <- "white"
			}else if (j > h){ # Pessimistic misclassifications
				polygon(x = c(j, j+1, j+1, j), y = 4 - c(h , h , h-1, h-1), border = NA, col = grey(0.8))
				fnt <- 1
				fntcol <- 1
			} else { # Correct status
				polygon(x = c(j, j+1, j+1, j), y = 4 - c(h , h , h-1, h-1), border = NA, col = statusCols[4-h])
				fnt <- 1
				fntcol <- 1
			}
			text(j + 0.5, 4 - h + 0.5, paste(formatC(round(statusDiffMat[h,j]/length(statusDiff)*100, 1), 1, format="f"), "%"), font=fnt, col=fntcol)
			# text(j + 0.5, 4 - h + 0.5, statusDiffMat[h,j])
		}
	}
	
	if(k==1 & all3 == TRUE){
		mtext(side=2, line = 6, c("High productivity\nHCR", "Low productivity\nModerate harvest", "Low productivity\nHigh harvest")[i], cex=0.8, font=2)
	}
	
	par(mar=c(0,5,4,0))
	TR <- apply(statusDiffMat, 2, sum)/length(statusDiff)
	bp <- barplot(TR, col=statusCols[c('r', 'a', 'g')], yaxt="n", names.arg = FALSE)
	text(bp, TR, pos=3, paste(formatC(round(TR*100, 1), 1, format="f"), "%"), xpd=NA)
	
	if(all3 == FALSE){
		mtext(side = 3, adj = 0, line=2, c("a) Spawner-recruitment", "b) Historical spawners")[k], cex=0.8)
	} else {
		mtext(side = 3, adj = 0, line=2, paste(letters[(i-1)*2 +k], ")", sep=""), cex=0.8)
		if(i==1) mtext(side=3, line=3.5, c("Spawner-recruitment", "Historical spawners")[k], cex=0.8, font = 2)
	}
	
	par(mar=c(5,0,0,4))
	OBS <- apply(statusDiffMat, 1, sum)/length(statusDiff)
	bp <- barplot(rev(OBS), col=rev(statusCols[c('r', 'a', 'g')]), yaxt="n", horiz=TRUE, xaxt="n", names.arg=FALSE)
	text(rev(OBS), bp, pos=4, paste(formatC(round(rev(OBS)*100, 1), 1, format="f"), "%"), xpd=NA)
	
	} # end k
	
} # end i

#------------------------------------------------------------------------------
# Why misclassifications? RB in benchmarks
#------------------------------------------------------------------------------

delistedOut <- delistSensitivity(out_base)

quartz(width = 6, height= 3, pointsize = 9)
# layout(matrix(c(1,2,3,3), nrow=2, byrow=2))
par(mfrow=c(1,2), mar=c(3, 4, 1, 1), oma=c(2,3.5,2,2))

# RBs

x <- c(1:3)
jit <- 0.1
ptCex <- 1

plotCI(x, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", xlim=c(0.5, 3.5))

axis(side=1, at=x, labels=FALSE)
u <- par('usr')
text(x, rep(u[3] - 0.1*(u[4]-u[3]), 3), c("high", "low", "low"), xpd=NA)
text(0, u[3] - 0.1*(u[4]-u[3]), "productivity:", xpd=NA, adj=1, font=2)
text(x, rep(u[3] - 0.2*(u[4]-u[3]), 3), c("HCR", "moderate", "high"), xpd=NA)
text(0, u[3] - 0.2*(u[4]-u[3]), "target harvest:", xpd=NA, adj=1, font=2)
abline(h=0)

plotCI(x - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)

plotCI(x + jit, delistedOut$RB$Smsy[, 'median'], li = delistedOut$RB$Smsy[, '25th'], ui=delistedOut$RB$Smsy[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])

plotCI(x, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Relative bias")

legend("topright", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN1])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), bty="n", pt.cex = ptCex, pt.bg = c(statusCols['r'], 1, statusCols['g']))

mtext(side=3, line=0.5, adj=0, "a) Spawner-recruitment")

# Historical spawners

yAll <- c(delistedOut$RB$S25[, '25th'], delistedOut$RB$S25[, '75th'], delistedOut$RB$S55[, '25th'], delistedOut$RB$S50[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th'])

plotCI(x, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", ylim=range(yAll), xlim=c(0.5, 3.5))

axis(side=1, at=x, labels=FALSE)
u <- par('usr')
text(x, rep(u[3] - 0.1*(u[4]-u[3]), 3), c("high", "low", "low"), xpd=NA)
text(x, rep(u[3] - 0.2*(u[4]-u[3]), 3), c("HCR", "moderate", "high"), xpd=NA)
abline(h=0)

plotCI(x - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)

plotCI(x + jit, delistedOut$RB$S50[, 'median'], li = delistedOut$RB$S50[, '25th'], ui=delistedOut$RB$S50[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])

plotCI(x, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)

mtext(side=3, line=0.5, adj=0, "b) Historical spawners")

legend("topright", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[25])), expression(italic(S[AVG])), expression(paste(italic(S[50])))), bty="n", pt.cex = ptCex, pt.bg = c(statusCols['r'], 1, statusCols['g']))

#------------------------------------------------------------------------------
# Comparing to 100% monitoring coverage of indicator, non-indicator, and both
#------------------------------------------------------------------------------

sensitivityPar <- c(1:4)
delistedOut <- readRDS("workspaces/base100Mon_delisted_baseGreen.rds")

quartz(width = 6, height= 3, pointsize = 9)
# layout(matrix(c(1,2,3,3), nrow=2, byrow=2))
par(mfrow=c(1,2), mar=c(3, 4, 1, 1), oma=c(2,2,2,2))

x <- c(1:4)
jit <- 0.1
ptCex <- 1
# Stock-recruitment

yAll <- c(delistedOut$RB$Sgen1[, '25th'], delistedOut$RB$Sgen1[, '75th'], delistedOut$RB$Smsy[, '25th'], delistedOut$RB$Smsy[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th'])


plotCI(x - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", xlim=c(0.5, 4.5), ylim=range(yAll))

u <- par('usr')
axis(side=1, at=x, labels=FALSE)
text(x, u[3] - 0.15*(u[4]-u[3]), c("Base", "100%\nind.", "100%\nnon-ind", "100%\nboth"), xpd=NA, cex=0.8)

abline(h=0)

plotCI(x - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)

plotCI(x + jit, delistedOut$RB$Smsy[, 'median'], li = delistedOut$RB$Smsy[, '25th'], ui=delistedOut$RB$Smsy[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])

plotCI(x, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Relative bias")

#legend(3, u[4] + 0.22*(u[4]-u[3]), pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN1])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), pt.cex = ptCex, bg="white", xpd=NA)

mtext(side=3, line=0.5, adj=0, "a) Spawner-recruitment")

# Historical spawners
yAll <- c(delistedOut$RB$S25[, '25th'], delistedOut$RB$S25[, '75th'], delistedOut$RB$S50[, '25th'], delistedOut$RB$S50[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th'])

plotCI(x - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", ylim=range(yAll), xlim=c(0.5, 4.5))

u <- par('usr')
axis(side=1, at=x, labels=FALSE)
text(x, u[3] - 0.15*(u[4]-u[3]), c("Base", "100%\nind.", "100%\nnon-ind", "100%\nboth"), xpd=NA, cex=0.8)


abline(h=0)

plotCI(x - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)

plotCI(x + jit, delistedOut$RB$S50[, 'median'], li = delistedOut$RB$S50[, '25th'], ui=delistedOut$RB$S50[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])

plotCI(x, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)

mtext(side=3, line=0.5, adj=0, "b) Historical spawners")

# legend(4.1, u[4] + 0.22*(u[4]-u[3]), pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[S25])), expression(italic(S[AVG])), expression(paste(italic(S[50])))), pt.cex = 0.8, bg="white", xpd=NA)

###############################################################################
# Monitoring coverage bias
###############################################################################

#------------------------------------------------------------------------------
# Example of each scenario
#------------------------------------------------------------------------------
# indCol <- c(ind = "#000000", nonInd = "#9C9C9C")
indCol <- c(ind = "#475A83", nonInd = "#C2B642")
indCol2 <- c(ind = "#475A8380", nonInd = "#C2B64280")

quartz(width = 6.5, height= 2.8, pointsize = 10)
# par(mfrow=c(1,2), mar=c(3,4,2,1))
layout(matrix(c(1,1,2,2,3), nrow=1))
par(mar=c(3,4,2,1), oma=c(0,0,1,0))

# Proportion monitored from English (2016)
plot(datMon$year, datMon$ppn_indicator, "n", las=1, bty="l", xlab="", ylab="Proportion monitored", ylim=c(0, 1))

points(datMon$year, datMon$ppn_nonindicator, col=indCol2['nonInd'], cex=0.6, pch=19)
points(datMon$year, datMon$ppn_indicator, col=indCol2['ind'], cex=0.6, pch=4)

abline(v = 1986+0.5)

# legend("bottomleft", col=indCol, pch = c(4, 1), c("Indicator", "Non-indicator"), bty="n", pt.cex=0.5)
mtext(side=3, adj = 0, "a) All systems", line=1)

# Constant
lines(datMon$year[c(period1, period2)], rep(preDecl[1], length(c(period1, period2))), col=indCol['ind'], lwd=2)
lines(datMon$year[c(period1, period2)], rep(preDecl[2], length(c(period1, period2))), col=indCol['nonInd'], lwd=2)

# Decline lines - base case
lines(datMon$year[c(period1, period2)], c(rep(preDecl[1], length(period1)), preDecl[1] + predict(fitInd)), col=indCol['ind'])
lines(datMon$year[c(period1, period2)], c(rep(preDecl[2], length(period1)), preDecl[2] + predict(fitNonInd)), col=indCol['nonInd'])

# Declines recent
lines(datMon$year[c(period1.b, period2.b)], c(rep(preDecl.b[1], length(period1.b)), preDecl.b[1] + predict(fitInd.b)), col=indCol['ind'], lty=2)
lines(datMon$year[c(period1.b, period2.b)], c(rep(preDecl.b[2], length(period1.b)), preDecl.b[2] + predict(fitNonInd.b)), col=indCol['nonInd'], lty=2)

# Proportion monitored from chumSR data
plot(year, ppnMonitored[,1], "n", las=1, bty="l", xlab="", ylab="Proportion monitored", ylim=c(0, 1))
abline(v = 1986.5)
points(year, ppnMonitored[,2],  col=indCol2['ind'], cex=0.6, pch=4)
points(year, ppnMonitored[,3],  col=indCol2['nonInd'], cex=0.6, pch=19)

# Decline lines
lines(year[c(period1.c, period2.c)], c(rep(preDecl.c[2], length(period1.c)), preDecl.c[2] + predict(fitInd.c)), col=indCol['ind'], lty=3, lwd=1.5)
lines(year[c(period1.c, period2.c)], c(rep(preDecl.c[3], length(period1.c)), preDecl.c[3] + predict(fitNonInd.c)), col=indCol['nonInd'], lty=3, lwd=1.5)


mtext(side=3, adj = 0, "b) Central coast chum rivers", line=1)

# legend
par(mar=c(0,0,0,0))
plot(1,1,"n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

legend("top", pch = c(4, 19, NA, NA, NA, NA), col=c(indCol2, rep(grey(0.6), 4)), lty=c(NA, NA, 1, 1, 2, 3), lwd=c(NA, NA, 2, 1, 1, 1.5), c("Indicator", "Non-indicator", "Scenario A", "Scenario B", "Scenario C" ,"Scenario D"), xpd=NA, bty="n")

#------------------------------------------------------------------------------
# Sensitivity results on performance
#------------------------------------------------------------------------------

sensitivityPar <- c(1:4)
delistedOut <- out_mon2
delistedOut2 <- out_monDecl
baseValue <- 2

quartz(width = 6, height= 4.5, pointsize = 10)
# layout(matrix(c(1,2,3,3), nrow=2, byrow=2))
par(mfrow=c(2,2), mar=c(3, 4, 1, 1), oma=c(2,2,2,0))

# Ppn Wrong

for(j in 1:2){
	y <- rbind(delistedOut$ppn[[j]][1,], delistedOut2$ppn[[j]][1,], delistedOut$ppn[[j]][2,], delistedOut2$ppn[[j]][2,], delistedOut$ppn[[j]][3,], delistedOut2$ppn[[j]][3,], delistedOut$ppn[[j]][4,], delistedOut2$ppn[[j]][4,])
						 
	bp <- barplot(t(y), col=c(statusCols, grey(0.8), 1), las=1, ylab = "", space=c(0.5, 0.1), border=NA, xaxs="i", yaxs="i")
	abline(h = 0)
	text(bp, 0.02, srt=90, rep(c("Base", "Decline"), 4), adj=0, col="white")
	text(apply(matrix(bp, nrow=2), 2, mean), rep(-0.1, 4), LETTERS[1:4], xpd=NA)
	points(apply(matrix(bp, nrow=2), 2, mean)[2], 1.05, pch=8, xpd=NA)
	if(j == 1) mtext(side=2, line = 4, "Proportion of simulations")
	mtext(side=3, line=2, c("Stock-recruitment", "Historical spawners")[j])
	mtext(side=3, line=0.5, adj=0, c("a)", "b)")[j])
}

# RBs

inter <- function(out1, out2){
	y <- rbind(out1[1,], out2[1,], out1[2,], out2[2,], out1[3,], out2[3,], out1[4,], out2[4,])
	return(y)
}
# x <- c(seq(0.85, 3.85, 1), seq(1.15, 4.15, 1))
x <- sort(c(seq(0.95, 3.95, 1), seq(1.05, 4.05, 1)))
jit <- 0.25

plotCI(x - jit, inter(delistedOut$RB$Sgen1, delistedOut2$RB$Sgen1)[, 'median'], "n", li = inter(delistedOut$RB$Sgen1, delistedOut2$RB$Sgen1)[, '25th'], ui= inter(delistedOut$RB$Sgen1, delistedOut2$RB$Sgen1)[, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim=c(-0.2, 1.1), xlim=c(0.5, 4.5), xaxt="n")
axis(side=1, at=1:4, labels = LETTERS[1:4])

polygon(x=c(1.5, 2.5, 2.5, 1.5), y=c(-0.25, -0.25, 1.15, 1.15), col="#00000020", border=NA)
# abline(v = baseValue, lwd=10, col="#00000020")
abline(h=0)

plotCI(x - jit, inter(delistedOut$RB$Sgen1, delistedOut2$RB$Sgen1)[, 'median'], li = inter(delistedOut$RB$Sgen1, delistedOut2$RB$Sgen1)[, '25th'], ui= inter(delistedOut$RB$Sgen1, delistedOut2$RB$Sgen1)[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=0.8, pch=25, pt.bg=c("white", statusCols['r']), add=TRUE)

plotCI(x + jit, inter(delistedOut$RB$Smsy, delistedOut2$RB$Smsy)[, 'median'], li = inter(delistedOut$RB$Smsy, delistedOut2$RB$Smsy)[, '25th'], ui= inter(delistedOut$RB$Smsy, delistedOut2$RB$Smsy)[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=0.8, pch=24, pt.bg=c("white", statusCols['g']))

plotCI(x, inter(delistedOut$RB$avgS, delistedOut2$RB$avgS)[, 'median'], li = inter(delistedOut$RB$avgS, delistedOut2$RB$avgS)[, '25th'], ui= inter(delistedOut$RB$avgS, delistedOut2$RB$avgS)[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=0.8, pch=21, pt.bg=c("white", 1))

points(rep(baseValue, 2), par('usr')[c(3,4)], pch=8, xpd=NA)

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Relative bias")

# axis(side=3, at=log(1/c(1:5)), labels=1:5)
# mtext(side=3, line=3, "Corresponding Expansion Factor III", cex=0.8)
mtext(side=3, line=0.5, adj=0, "c)")
# legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN1])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), bty="n", pt.cex = 0.8)

# Historical spawners
plotCI(x - jit, inter(delistedOut$RB$S25, delistedOut2$RB$S25)[, 'median'], "n", li = inter(delistedOut$RB$S25, delistedOut2$RB$S25)[, '25th'], ui= inter(delistedOut$RB$S25, delistedOut2$RB$S25)[, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim=c(-0.2, 3.2), xlim=c(0.5, 4.5), xaxt="n")
axis(side=1, at=1:4, labels = LETTERS[1:4])

u <- par('usr')
polygon(x=c(1.5, 2.5, 2.5, 1.5), y=rep(u[3:4], each=2), col="#00000020", border=NA)
# abline(v = baseValue, lwd=10, col="#00000020")
abline(h=0)

plotCI(x - jit, inter(delistedOut$RB$S25, delistedOut2$RB$S25)[, 'median'], li = inter(delistedOut$RB$S25, delistedOut2$RB$S25)[, '25th'], ui= inter(delistedOut$RB$S25, delistedOut2$RB$S25)[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=0.8, pch=25, pt.bg=c("white", statusCols['r']), add=TRUE)

plotCI(x + jit, inter(delistedOut$RB$S50, delistedOut2$RB$S50)[, 'median'], li = inter(delistedOut$RB$S50, delistedOut2$RB$S50)[, '25th'], ui= inter(delistedOut$RB$S50, delistedOut2$RB$S50)[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=0.8, pch=24, pt.bg=c("white", statusCols['g']))

plotCI(x, inter(delistedOut$RB$avgS, delistedOut2$RB$avgS)[, 'median'], li = inter(delistedOut$RB$avgS, delistedOut2$RB$avgS)[, '25th'], ui= inter(delistedOut$RB$avgS, delistedOut2$RB$avgS)[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=0.8, pch=21, pt.bg=c("white", 1))

points(rep(baseValue, 2), par('usr')[c(3,4)], pch=8, xpd=NA)

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
# mtext(side=2, line=4, "Relative bias")

# axis(side=3, at=log(1/c(1:5)), labels=1:5)
# mtext(side=3, line=3, "Corresponding Expansion Factor III", cex=0.8)
mtext(side=3, line=0.5, adj=0, "d)")
# legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[25])), expression(italic(S[AVG])), expression(italic(S[50]))), bty="n", pt.cex = 0.8)


mtext(side=1, outer=TRUE, "Monitoring scenario", line=0)



###############################################################################
# Number of indicator & non-indicator streams
###############################################################################
# Cases:
# 1) Base: 15 ind + 20 non-ind = 35
# 2) SmallLow: 3 in + 7 non-ind = 10
# 3) SmallHigh: 8 ind + 2 non-ind = 10
# 4) LargeLow: 42 ind + 58 non-ind = 140
# 5) LargeHigh: 119 ind + 21 non-ind = 140

sensitivityPar <- c(1:5)
delistedOut <- out_nPop2
baseValue <- 1

quartz(width = 6, height= 4.5, pointsize = 10)
# layout(matrix(c(1,2,3,3), nrow=2, byrow=2))
par(mfrow=c(2,2), mar=c(3, 4, 1, 1), oma=c(3,2,2,0))


# Ppn Wrong


for(j in 1:2){
	bp <- barplot(t(delistedOut$ppn[[j]]), 
								col=c(statusCols, grey(0.8), 1), las=1, ylab = "", 
								space=c(0.5, 0.5, 0.1, 0.5, 0.1),
								border=NA, yaxs="i")
	abline(h = 0)
	
	# text(bp, 0.02, srt=90, rep(c("Base", "Decline"), 4), adj=0, col="white")
	
	text(bp, rep(-0.1, 5), LETTERS[1:5], xpd=NA)
	
	points(bp[1], 1.05, pch=8, xpd=NA) # Base case
	
	if(j == 1) mtext(side=2, line = 4, "Proportion of simulations")
	
	mtext(side=3, line=2, c("Stock-recruitment", "Historical spawners")[j])
	mtext(side=3, line=0.5, adj=0, c("a)", "b)")[j])
}

uBP <- par('usr')

# RBs

x <- bp
jit <- 0.1

plotCI(x - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", xlim=uBP[1:2])

axis(side=1, at=x, LETTERS[1:5])
u <- par('usr')
polygon(x=x[1] + c(-2, 2, 2, -2)*jit, y = rep(u[3:4], each=2), col="#00000020", border=NA)
abline(h=0)

plotCI(x - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=0.8, pch=25, pt.bg=c("white", rep(c("white", statusCols['r']), 2)), add=TRUE)

plotCI(x + jit, delistedOut$RB$Smsy[, 'median'], li = delistedOut$RB$Smsy[, '25th'], ui=delistedOut$RB$Smsy[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=0.8, pch=24, pt.bg=c("white", rep(c("white", statusCols['g']), 2)))

plotCI(x, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=0.8, pch=21, pt.bg=c("white", rep(c("white", 1), 2)))

points(rep(x[1], 2), par('usr')[c(3,4)], pch=8, xpd=NA)

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Relative bias")

legend("topright", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN1])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), bty="n", pt.cex = 0.8)

mtext(side=3, line=0.5, adj=0, "c)")

u <- par('usr')
y <- u[3] - 0.25*(u[4] - u[3])
segments(x0=x[1]-jit, x1=x[1]+jit, y0=y, y1=y, xpd=NA); text(x[1], y, pos=1, "Base", xpd=NA)
segments(x0=x[2] - jit, x1=x[3] + jit, y0=y, y1=y, xpd=NA); text(mean(x[2:3]), y, pos=1, "Small CUs", xpd=NA)
segments(x0=x[4] - jit, x1=x[5] + jit, y0=y, y1=y, xpd=NA); text(mean(x[4:5]), y, pos=1, "Large CUs", xpd=NA)

# Historical spawners

yAll <- c(delistedOut$RB$S25[, '25th'], delistedOut$RB$S25[, '75th'], delistedOut$RB$S55[, '25th'], delistedOut$RB$S50[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th'])

plotCI(x - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", ylim=range(yAll), xlim=uBP[1:2])

axis(side=1, at=x, LETTERS[1:5])
u <- par('usr')
polygon(x=x[1] + c(-2, 2, 2, -2)*jit, y = rep(u[3:4], each=2), col="#00000020", border=NA)
abline(h=0)

plotCI(x - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=0.8, pch=25, pt.bg=c("white", rep(c("white", statusCols['r']), 2)), add=TRUE)

plotCI(x + jit, delistedOut$RB$S50[, 'median'], li = delistedOut$RB$S50[, '25th'], ui=delistedOut$RB$S50[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=0.8, pch=24, pt.bg=c("white", rep(c("white", statusCols['g']), 2)))

plotCI(x, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=0.8, pch=21, pt.bg=c("white", rep(c("white", 1), 2)))

points(rep(x[1], 2), par('usr')[c(3,4)], pch=8, xpd=NA)

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)

mtext(side=3, line=0.5, adj=0, "d)")

u <- par('usr')
y <- u[3] - 0.25*(u[4] - u[3])
segments(x0=x[1]-jit, x1=x[1]+jit, y0=y, y1=y, xpd=NA); text(x[1], y, pos=1, "Base", xpd=NA)
segments(x0=x[2] - jit, x1=x[3] + jit, y0=y, y1=y, xpd=NA); text(mean(x[2:3]), y, pos=1, "Small CUs", xpd=NA)
segments(x0=x[4] - jit, x1=x[5] + jit, y0=y, y1=y, xpd=NA); text(mean(x[4:5]), y, pos=1, "Large CUs", xpd=NA)

mtext(side=1, outer=TRUE, "Scenario for number of indicator & non-indicator streams", line=1)

###############################################################################
# Change in capacity
###############################################################################

declCap <- amberHab + redHab
jit <- diff(declCap)[1]/10

delistedOut <- out_capacity2
baseValue <- simPar$greenHab

quartz(width = 6, height= 4.5, pointsize = 10)
# layout(matrix(c(1,2,3,3), nrow=2, byrow=2))
par(mfrow=c(2,2), mar=c(3, 4, 1, 1), oma=c(2,2,2,0))


# Ppn Wrong


for(j in 1:2){
	
	bp <- barplot(t(delistedOut$ppn[[j]]), col=c(statusCols, grey(0.8), 1), las=1, ylab = "", space=0, border=NA, xaxs="i", yaxs="i")
	axis(side = 1, at = bp[seq(1, 11, 2)], labels=declCap[seq(1, 11, 2)])
	points(bp[which(greenHab == baseValue)], 1.05, pch=8, xpd=NA)
	if(j == 1) mtext(side=2, line = 4, "Proportion of simulations")
	mtext(side=3, line=2, c("Stock-recruitment", "Historical spawners")[j])
	mtext(side=3, line=0.5, adj=0, c("a)", "b)")[j])
}

# RBs

ylims <- c(-0.2, 1.3)

plotCI(declCap - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim = ylims)

abline(v = baseValue, lwd=10, col="#00000020")
abline(h=0)

plotCI(declCap - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=0.8, pch=25, pt.bg="white", add=TRUE)
plotCI(declCap + jit, delistedOut$RB$Smsy[, 'median'], li = delistedOut$RB$Smsy[, '25th'], ui=delistedOut$RB$Smsy[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=0.8, pch=24, pt.bg="white")
plotCI(declCap, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=0.8, pch=21, pt.bg="white")

points(rep(baseValue, 2), par('usr')[c(3,4)], pch=8, xpd=NA)
A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Relative bias")

mtext(side=3, line=0.5, adj=0, "c)")
legend("topright", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN1])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), bty="n", pt.cex = 0.8)

# Historical spawners
# ylims <- range(c(delistedOut$RB$S25[, '25th'], delistedOut$RB$S25[, '75th'], delistedOut$RB$S50[, '25th'], delistedOut$RB$S50[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th']))
ylims <- c(-0.2, 4)
plotCI(declCap - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim =ylims)

abline(v = baseValue, lwd=10, col="#00000020")
abline(h=0)

plotCI(declCap - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=0.8, pch=25, pt.bg="white", add=TRUE)
plotCI(declCap + jit, delistedOut$RB$S50[, 'median'], li = delistedOut$RB$S50[, '25th'], ui=delistedOut$RB$S50[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=0.8, pch=24, pt.bg="white")
plotCI(declCap, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=0.8, pch=21, pt.bg="white")

points(rep(baseValue, 2), par('usr')[c(3,4)], pch=8, xpd=NA)

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)

mtext(side=3, line=0.5, adj=0, "d)")

legend("topright", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[25])), expression(italic(S[AVG])), expression(italic(S[50]))), bty="n", pt.cex = 0.8)


mtext(side=1, outer=TRUE, "Proportion of subpopulations with declining capacity", line=0)

###############################################################################
# Observation bias
###############################################################################

sensitivityPar <- obsBias
delistedOut <- out_obsBias2
baseValue <- simPar$obs_bias

quartz(width = 6, height= 4.5, pointsize = 10)
# layout(matrix(c(1,2,3,3), nrow=2, byrow=2))
par(mfrow=c(2,2), mar=c(3, 4, 1, 1), oma=c(2,2,2,0))


# Ppn Wrong


for(j in 1:2){
	
	bp <- barplot(t(delistedOut$ppn[[j]]), col=c(statusCols, grey(0.8), 1), las=1, ylab = "", space=0.1, border=NA, xaxs="i", yaxs="i")
	
	# axis(side = 1 , at = predict(lm(bp ~ sensitivityPar), newdata = data.frame(sensitivityPar = seq(-1.5, 0, 0.5))), labels=formatC(seq(-1.5, 0, 0.5), 1, format="f"))
	axis(side = 1 , at = predict(lm(bp ~ sensitivityPar), newdata = data.frame(sensitivityPar = seq(-1.6, 0, 0.4))), labels=formatC(seq(-1.6, 0, 0.4), 1, format="f"))
	
	points(predict(lm(bp ~ sensitivityPar), newdata = data.frame(sensitivityPar = baseValue)), 1.05, pch=8, xpd=NA)
	if(j == 1) mtext(side=2, line = 4, "Proportion of simulations")
	mtext(side=3, line=2, c("Stock-recruitment", "Historical spawners")[j])
	mtext(side=3, line=0.5, adj=0, c("a)", "b)")[j])
	}


# RBs

jit <- 0.02

plotCI(sensitivityPar - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xlim=c(-1.7, 0.1), xaxt="n")
axis(side=1, at=seq(-1.6, 0, 0.4))

abline(v = baseValue, lwd=16, col=grey(0.8))

abline(h=0)

plotCI(sensitivityPar - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=0.8, pch=25, pt.bg="white", add=TRUE)

plotCI(sensitivityPar + jit, delistedOut$RB$Smsy[, 'median'], li = delistedOut$RB$Smsy[, '25th'], ui=delistedOut$RB$Smsy[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=0.8, pch=24, pt.bg="white")

plotCI(sensitivityPar, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=0.8, pch=21, pt.bg="white")

points(rep(baseValue, 2), par('usr')[c(3,4)], pch=8, xpd=NA)

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Relative bias")

# axis(side=3, at=log(1/c(1:5)), labels=1:5)
# mtext(side=3, line=3, "Corresponding Expansion Factor III", cex=0.8)
mtext(side=3, line=0.5, adj=0, "c)")
legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN1])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), bty="n", pt.cex = 0.8)

# Historical spawners
plotCI(sensitivityPar - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim=range(c(delistedOut$RB$S25[, '25th'], delistedOut$RB$S25[, '75th'], delistedOut$RB$S50[, '25th'], delistedOut$RB$S50[, '75th'])), xlim=c(-1.7, 0.1), xaxt="n")
axis(side=1, at=seq(-1.6, 0, 0.4))

abline(v = baseValue, lwd=16, col=grey(0.8))
abline(h=0)

plotCI(sensitivityPar - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=0.8, pch=25, pt.bg="white", add=TRUE)
plotCI(sensitivityPar + jit, delistedOut$RB$S50[, 'median'], li = delistedOut$RB$S50[, '25th'], ui=delistedOut$RB$S50[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=0.8, pch=24, pt.bg="white")
plotCI(sensitivityPar, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=0.8, pch=21, pt.bg="white")
points(rep(baseValue, 2), par('usr')[c(3,4)], pch=8, xpd=NA)
A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
# mtext(side=2, line=4, "Relative Bias")

# axis(side=3, at=log(1/c(1:5)), labels=1:5)
# mtext(side=3, line=3, "Corresponding Expansion Factor III", cex=0.8)
mtext(side=3, line=0.5, adj=0, "d)")
legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[25])), expression(italic(S[AVG])), expression(italic(S[50]))), bty="n", pt.cex = 0.8)


mtext(side=1, outer=TRUE, expression(paste("Observation bias of spawners (", bar(delta), ")")), line=1)


###############################################################################
# Catch bias
###############################################################################

# Don't need to include historical spawners here because catch bias does not 
# affect spawner estimates, only SR relationship.

sensitivityPar <- catchBias
delistedOut <- out_catchBias2
baseValue <- simPar$catch_bias

quartz(width = 3.2, height= 4.5, pointsize = 10)
# layout(matrix(c(1,2,3,3), nrow=2, byrow=2))
par(mfrow=c(2,1), mar=c(3, 4, 1, 1), oma=c(1,1,2,0))


# Ppn Wrong


j <- 1
# for(j in 1:2){
	
	bp <- barplot(t(delistedOut$ppn[[j]]), col=c(statusCols, grey(0.8), 1), las=1, ylab = "", space=0.1, border=NA, xaxs="i", yaxs="i")
	
	axis(side = 1, at = bp[seq(1, length(bp), 2)], labels=formatC(seq(-1, 1, 0.4), 1, format="f"))
	
	points(bp[findInterval(baseValue, sensitivityPar)], 1.05, pch=8, xpd=NA)
	if(j == 1) mtext(side=2, line = 4, "Proportion of simulations")
	mtext(side=3, line=2, c("Stock-recruitment", "Historical spawners")[j])
	mtext(side=3, line=0.5, adj=0, c("a)", "b)")[j])
# }
	

# RBs

	jit <- 0.015
	
plotCI(sensitivityPar - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", ylim=c(-0.3, 1.4))
axis(side=1, at = seq(-1, 1, 0.4))
abline(v = baseValue, lwd=10, col="#00000020")
abline(h=0)
plotCI(sensitivityPar - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=0.8, pch=25, pt.bg="white", add=TRUE)
plotCI(sensitivityPar + jit, delistedOut$RB$Smsy[, 'median'], li = delistedOut$RB$Smsy[, '25th'], ui=delistedOut$RB$Smsy[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=0.8, pch=24, pt.bg="white")
plotCI(sensitivityPar, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=0.8, pch=21, pt.bg="white")
points(rep(baseValue, 2), par('usr')[c(3,4)], pch=8, xpd=NA)
A <- axis(side=2, labels = FALSE, at=c(0, 0.4, 0.8, 1.2))
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Relative bias")

# axis(side=3, at=log(1/c(1:5)), labels=1:5)
# mtext(side=3, line=3, "Corresponding Expansion Factor III", cex=0.8)
if(j==1){
	mtext(side=3, line=0.5, adj=0, "b)") 
	mtext(side=1, expression(paste("Observation bias of catch (", bar(chi), ")")), line=3)
	
} else {
	mtext(side=3, line=0.5, adj=0, "c)")
}
# legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN1])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), bty="n", pt.cex = 0.8)

# # Historical spawners
# ylims <- range(c(delistedOut$RB$S25[, '25th'], delistedOut$RB$S25[, '75th'], delistedOut$RB$S50[, '25th'], delistedOut$RB$S50[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th']))
# 
# plotCI(sensitivityPar - 0.01, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim=ylims)
# abline(v = baseValue, lwd=10, col="#00000020")
# abline(h=0)
# plotCI(sensitivityPar - 0.01, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=0.8, pch=25, pt.bg="white", add=TRUE)
# plotCI(sensitivityPar + 0.01, delistedOut$RB$S50[, 'median'], li = delistedOut$RB$S50[, '25th'], ui=delistedOut$RB$S50[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=0.8, pch=24, pt.bg="white")
# plotCI(sensitivityPar, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=0.8, pch=21, pt.bg="white")
# points(rep(baseValue, 2), par('usr')[c(3,4)], pch=8, xpd=NA)
# A <- axis(side=2, labels = FALSE)
# axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
# # mtext(side=2, line=4, "Relative Bias")
# 
# # axis(side=3, at=log(1/c(1:5)), labels=1:5)
# # mtext(side=3, line=3, "Corresponding Expansion Factor III", cex=0.8)
# mtext(side=3, line=0.5, adj=0, "d)")
# # legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[25])), expression(italic(S[AVG])), expression(italic(S[50]))), bty="n", pt.cex = 0.8)
# 

# mtext(side=1, outer=TRUE, expression(paste("Catch bias (", bar(chi), ")")), line=1)

###############################################################################
# Capacity X Monitoring
###############################################################################
#------------------------------------------------------------------------------
# Scenarios
#------------------------------------------------------------------------------
par(mfrow=c(1,1), mar=c(4,4,2,5))
plot(1:50, 1:50, "n", ylim=c(0.3, 0.9), las=1, xlab="Year of simulation", ylab="Proportion monitored", bty="l")
segments(1, 0.8, 50, 0.8)
segments(1, 0.4, 50, 0.4, lty=2)
lines(c(1:50), c(rep(0.8, 40), seq(0.8, 0.4, length.out=10)), lty=3, col=statusCols['r'])
lines(c(1:50), seq(0.8, 0.4, length.out=50), lty=2, lwd=2, col=statusCols['r'])
text(50, 0.8, "A", xpd=NA, pos=3)
text(46, 0.65, "D", col=statusCols['r'])
text(10, 0.4, pos=3, "B")
text(25, 0.6, "C", col=statusCols['r'], pos=3)

#------------------------------------------------------------------------------
# Sensitivity of performance
#------------------------------------------------------------------------------

declCap <- amberHab + redHab
jit <- 0.2

delistedOut <- out_capacityXmon2

quartz(width = 8, height= 4.5, pointsize = 10)
# layout(matrix(c(1,2,3,3), nrow=2, byrow=2))
par(mfrow=c(2,2), mar=c(3, 4, 1, 1), oma=c(2,2,2,0))


# Ppn Wrong


for(j in 1:2){
	
	bp <- barplot(t(delistedOut$ppn[[j]]), col=c(statusCols, grey(0.8), 1), las=1, ylab = "", border=NA, xaxs="i", yaxs="i", space=rep(c(0.6, rep(0.1, 5)), 4))
	abline(h=0)
	
	axis(side=1, at=bp, labels=FALSE, tck=-0.02, line=-0.1)
	uBP <- par('usr')
	text(c(mean(bp[1:6]), mean(bp[7:12]), mean(bp[13:18]), mean(bp[19:24])), rep(uBP[3], 4), pos=1, LETTERS[1:4], xpd=NA)
	
	if(j == 1) mtext(side=2, line = 4, "Proportion of simulations")
	mtext(side=3, line=2, c("Stock-recruitment", "Historical spawners")[j])
	mtext(side=3, line=0.5, adj=0, c("a)", "b)")[j])
}

# RBs

ylims <- range(c(delistedOut$RB$Sgen1[, '25th'], delistedOut$RB$Sgen1[, '75th'], delistedOut$RB$Smsy[, '25th'], delistedOut$RB$Smsy[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th']))

plotCI(bp - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim = ylims, xaxt="n")

mtext(side=3, "Monitoring scenario", line=1.5)
text(c(mean(bp[1:6]), mean(bp[7:12]), mean(bp[13:18]), mean(bp[19:24])), rep(ylims[2], 4), pos=3, LETTERS[1:4], xpd=NA)

abline(h=0)

plotCI(bp - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=0.8, pch=25, pt.bg="white", add=TRUE)
plotCI(bp + jit, delistedOut$RB$Smsy[, 'median'], li = delistedOut$RB$Smsy[, '25th'], ui=delistedOut$RB$Smsy[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=0.8, pch=24, pt.bg="white")
plotCI(bp, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=0.8, pch=21, pt.bg="white")

abline(v=c(mean(bp[6:7]), mean(bp[12:13]), mean(bp[18:19])))

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Relative bias")

axis(side=1, at=bp, labels=FALSE, tck=-0.02)
axis(side=1, at=bp[seq(1,24,2)], labels=rep(c(0,40,80), 4))

mtext(side=3, line=0.5, adj=0, "c)")


# Historical spawners
ylims <- range(c(delistedOut$RB$S25[, '25th'], delistedOut$RB$S25[, '75th'], delistedOut$RB$S50[, '25th'], delistedOut$RB$S50[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th']))

plotCI(bp - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim =ylims, xaxt="n")

abline(h=0)
mtext(side=3, "Monitoring scenario", line=1.5)
text(c(mean(bp[1:6]), mean(bp[7:12]), mean(bp[13:18]), mean(bp[19:24])), rep(ylims[2], 4), pos=3, LETTERS[1:4], xpd=NA)
abline(v=c(mean(bp[6:7]), mean(bp[12:13]), mean(bp[18:19])))


plotCI(bp - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=0.8, pch=25, pt.bg="white", add=TRUE)
plotCI(bp + jit, delistedOut$RB$S50[, 'median'], li = delistedOut$RB$S50[, '25th'], ui=delistedOut$RB$S50[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=0.8, pch=24, pt.bg="white")
plotCI(bp, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=0.8, pch=21, pt.bg="white")

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)

axis(side=1, at=bp, labels=FALSE, tck=-0.02)
axis(side=1, at=bp[seq(1,24,2)], labels=rep(c(0,40,80), 4))

mtext(side=3, line=0.5, adj=0, "d)")


mtext(side=1, outer=TRUE, "Proportion of subpopulations with declining capacity", line=0)

###############################################################################
# Capacity X Real monitoring scenarios
###############################################################################

redHab <- c(seq(0, 50, 25), 100)
declCap <- redHab
jit <- 0.2

delistedOut <- readRDS("workspaces/capacityXmon_delisted_baseGreen.rds")

quartz(width = 8, height= 5, pointsize = 10)
# layout(matrix(c(1,2,3,3), nrow=2, byrow=2))
par(mfrow=c(2,2), mar=c(4, 4, 1, 1), oma=c(1,2,2,0))

#------------------------------------------------------------------------------
# Ppn Wrong
#------------------------------------------------------------------------------

for(j in 1:2){
	
	bp <- barplot(t(delistedOut$ppn[[j]]), col=c(statusCols, grey(0.8), 1), las=1, ylab = "", border=NA, xaxs="i", yaxs="i", space=rep(c(0.6, rep(0.1, 3)), 4))
	
	abline(h=0)
	
	uBP <- par('usr')
	monLab <- c(mean(bp[2:3]), mean(bp[6:7]), mean(bp[10:11]), mean(bp[14:15]))
	# axis(side = 3, at = c(mean(bp[4:5]), mean(bp[8:9]), mean(bp[12:13])), labels=FALSE)
	text(monLab, rep(uBP[4], 4), pos=3, xpd=NA, LETTERS[1:4])
	
	# axis(side=1, at=bp, labels=FALSE, tck=-0.02)
	axis(side=1, at=bp, labels=rep(redHab, 4), tck=-0.02, las=2, mgp=c(3, 0.5, -0.04))
	
	if(j == 1) mtext(side=2, line = 4, "Proportion of simulations")
	mtext(side=3, line=2, c("Stock-recruitment", "Historical spawners")[j])
	mtext(side=3, line=1, adj=0, c("a)", "b)")[j])
}
#------------------------------------------------------------------------------
# RBs
#------------------------------------------------------------------------------
ptCex <- 1

ylims <- range(c(delistedOut$RB$Sgen1[, '25th'], delistedOut$RB$Sgen1[, '75th'], delistedOut$RB$Smsy[, '25th'], delistedOut$RB$Smsy[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th']))
# ylims <- c(-0.2, 1.5)
plotCI(bp - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim = ylims, xaxt="n")

u <- par('usr')
# mtext(side=3, "Monitoring scenario", line=1.5)
axis(side = 3, at = c(mean(bp[4:5]), mean(bp[8:9]), mean(bp[12:13])), labels=FALSE)
text(monLab, rep(u[4], 4), pos=3, xpd=NA, LETTERS[1:4])
		 
abline(h=0)

plotCI(bp - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)

plotCI(bp + jit, delistedOut$RB$Smsy[, 'median'], li = delistedOut$RB$Smsy[, '25th'], ui=delistedOut$RB$Smsy[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])

plotCI(bp, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)

abline(v=c(mean(bp[4:5]), mean(bp[8:9]), mean(bp[12:13])))

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Relative bias")

axis(side=1, at=bp, labels=rep(redHab, 4), tck=-0.02, mgp=c(3, 0.5, 0), las=2)

mtext(side=3, line=0.5, adj=0, "c)")

# legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN1])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), pt.cex = ptCex, xpd=NA, ncol=3, bg="white", pt.bg = c(statusCols['r'], 1, statusCols['g']))

#--------------------------
# Historical spawners
ylims <- range(c(delistedOut$RB$S25[, '25th'], delistedOut$RB$S25[, '75th'], delistedOut$RB$S50[, '25th'], delistedOut$RB$S50[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th']))
# ylims <- c(-0.2, 4.5)
plotCI(bp - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim =ylims, xaxt="n")

abline(h=0)

u <- par('usr')
# mtext(side=3, "Monitoring scenario", line=1.5)
axis(side = 3, at = c(mean(bp[4:5]), mean(bp[8:9]), mean(bp[12:13])), labels=FALSE)
text(monLab, rep(u[4], 4), pos=3, xpd=NA, LETTERS[1:4])


plotCI(bp - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], pch=25, pt.bg=statusCols['r'], add=TRUE, cex = ptCex)
plotCI(bp + jit, delistedOut$RB$S50[, 'median'], li = delistedOut$RB$S50[, '25th'], ui=delistedOut$RB$S50[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])
plotCI(bp, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)

abline(v=c(mean(bp[4:5]), mean(bp[8:9]), mean(bp[12:13])))

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)

axis(side=1, at=bp, labels=rep(redHab, 4), tck=-0.02, las=2, mgp=c(3, 0.5, 0))

mtext(side=3, line=0.5, adj=0, "d)")

# legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[25])), expression(italic(S[AVG])), expression(paste(italic(S[50])))), pt.cex = ptCex, xpd=NA, ncol = 3, bg="white", pt.bg = c(statusCols['r'], 1, statusCols['g']))

mtext(side=1, outer=TRUE, "Percentage of subpopulations with severe declines in capacity", line=-1)

#-----------------------------------------------------------------------------
# Summary of misclassifications
#-----------------------------------------------------------------------------
i <- 4
layout(matrix(c(2,1,1,2,1,1,7,3,3,5,4,4,5,4,4,8,6,6), nrow=3))
for(k in 1:2){
	# Base case results
	# statusDiff <- out_capacityXmon[[i]][[k+1]][, 3] # Under monitoring scenarios A-D
	statusDiff <- out_capacityXmon_100cov[[i]][[k+1]][, 3] # with 100% monitoring coverage
	# statusDiff <- out_base[[i]]$HS[, 3]
	statusDiffMat <- matrix(NA, nrow = 3, ncol = 3, dimnames=list(c("R", "A", "G"), c("R", "A", "G")))
	statusDiffMat[1, 1] <- length(which(statusDiff == 3)) #/ nSim
	statusDiffMat[2, 2] <- length(which(statusDiff == 2)) #/ nSim
	statusDiffMat[3, 3] <- length(which(statusDiff == 1)) #/ nSim
	statusDiffMat[1, 2] <- length(which(statusDiff == 4)) #/ nSim
	statusDiffMat[1, 3] <- length(which(statusDiff == 5)) #/ nSim
	statusDiffMat[2, 1] <- length(which(statusDiff == 6))
	statusDiffMat[2, 3] <- length(which(statusDiff == 7))
	statusDiffMat[3, 1] <- length(which(statusDiff == 8))
	statusDiffMat[3, 2] <- length(which(statusDiff == 9))
	
	# Plot percentages in each category
	par(mar=c(5,5,0,0))
	plot(1:4, 1:4, "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
	for(s in 1:2){
		axis(side = s, at = c(1:4), labels=FALSE)
	}
	axis(side = 1, at = seq(1.5, 3.5, 1), c("Red", "Amber", "Green"), tck = 0, las=1)
	axis(side = 2, at = rev(seq(1.5, 3.5, 1)), c("Red", "Amber", "Green"), tck = 0, las=1)
	
	mtext(side = 2, line = 4, "Estimated status", cex=0.8)
	mtext(side = 1, line = 3.5, "True status", cex=0.8)
	for(h in 1:3){
		for(j in 1:3){
			if(j < h){ # Optimistic misclassifications
				polygon(x = c(j, j+1, j+1, j), y = 4 - c(h , h , h-1, h-1), border = NA, col = 1)
				fnt <- 2 
				fntcol <- "white"
			}else if (j > h){ # Pessimistic misclassifications
				polygon(x = c(j, j+1, j+1, j), y = 4 - c(h , h , h-1, h-1), border = NA, col = grey(0.8))
				fnt <- 1
				fntcol <- 1
			} else { # Correct status
				polygon(x = c(j, j+1, j+1, j), y = 4 - c(h , h , h-1, h-1), border = NA, col = statusCols[4-h])
				fnt <- 1
				fntcol <- 1
			}
			text(j + 0.5, 4 - h + 0.5, paste(formatC(round(statusDiffMat[h,j]/length(statusDiff)*100, 1), 1, format="f"), "%"), font=fnt, col=fntcol)
			# text(j + 0.5, 4 - h + 0.5, statusDiffMat[h,j])
		}
	}
	
	par(mar=c(0,5,4,0))
	TR <- apply(statusDiffMat, 2, sum)/length(statusDiff)
	bp <- barplot(TR, col=statusCols[c('r', 'a', 'g')], yaxt="n", names.arg = FALSE)
	text(bp, TR, pos=3, paste(formatC(round(TR*100, 1), 1, format="f"), "%"), xpd=NA)
	mtext(side = 3, adj = 0, line=2, c("a) Stock-recruitment", "b) Historical spawners")[k], cex=0.8)
	
	par(mar=c(5,0,0,4))
	OBS <- apply(statusDiffMat, 1, sum)/length(statusDiff)
	bp <- barplot(rev(OBS), col=rev(statusCols[c('r', 'a', 'g')]), yaxt="n", horiz=TRUE, xaxt="n", names.arg=FALSE)
	text(rev(OBS), bp, pos=4, paste(formatC(round(rev(OBS)*100, 1), 1, format="f"), "%"), xpd=NA)
} # end k


###############################################################################
# Changes in observation bias part-way through time series 
###############################################################################

#------------------------------------------------------------------------------
# Example of scenarios
#------------------------------------------------------------------------------
par(mar=c(4,5,2,5), oma=c(0,0,0,0))
plot(1:50, rep(-0.4, 50), "l", ylim=c(-1.6, 0), ylab=expression(paste("Observation bias (", bar(delta), ")")), las=1, xaxs="i", lwd=3, col=grey(0.8), xlab="Year in simulation")
mtext(side=4, line=3.5, "Matching Expansion Factor III")
axis(side=4, at=log(1/c(1:5)), labels=1:5, las=1)

lines(c(1:50), rep(c(-0.4, -1.6), each=25))
lines(c(1:50), rep(c(-0.4, -0.7), each=25), lty=2)
lines(c(1:50), rep(c(-0.4, 0), each=25), lty=3)

#------------------------------------------------------------------------------

delistedOut <- out_obsBiasChange2

quartz(width = 6, height= 4.5, pointsize = 10)
# layout(matrix(c(1,2,3,3), nrow=2, byrow=2))
par(mfrow=c(2,2), mar=c(3, 4, 1, 1), oma=c(2,2,2,0))

#------------------------------------------------------------------------------
# Ppn Wrong
#------------------------------------------------------------------------------

for(j in 1:2){
	
	bp <- barplot(t(delistedOut$ppn[[j]]), col=c(statusCols, grey(0.8), 1), las=1, ylab = "", border=NA, xaxs="i", yaxs="i", space=0.1)
	
	abline(h=0)
	
	uBP <- par('usr')
	
	# axis(side=1, at=bp, labels=FALSE, tck=-0.02)
	axis(side=1, at=bp, labels=obs_bias2)
	
	if(j == 1) mtext(side=2, line = 4, "Proportion of simulations")
	mtext(side=3, line=2, c("Stock-recruitment", "Historical spawners")[j])
	mtext(side=3, line=1, adj=0, c("a)", "b)")[j])
}
#------------------------------------------------------------------------------
# RBs
#------------------------------------------------------------------------------
jit <- 0.1
ylims <- range(c(delistedOut$RB$Sgen1[, '25th'], delistedOut$RB$Sgen1[, '75th'], delistedOut$RB$Smsy[, '25th'], delistedOut$RB$Smsy[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th']))

plotCI(bp - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim = ylims, xaxt="n", xlim=c(0, 4.5))

u <- par('usr')

abline(h=0)

plotCI(bp - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=0.8, pch=25, pt.bg="white", add=TRUE)

plotCI(bp + jit, delistedOut$RB$Smsy[, 'median'], li = delistedOut$RB$Smsy[, '25th'], ui=delistedOut$RB$Smsy[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=0.8, pch=24, pt.bg="white")

plotCI(bp, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=0.8, pch=21, pt.bg="white")

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Relative bias")

axis(side=1, at=bp, labels=obs_bias2)

mtext(side=3, line=0.5, adj=0, "c)")

# legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN1])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), pt.cex = 0.8, xpd=NA, ncol=3, bg="white")

#--------------------------
# Historical spawners
ylims <- range(c(delistedOut$RB$S25[, '25th'], delistedOut$RB$S25[, '75th'], delistedOut$RB$S50[, '25th'], delistedOut$RB$S50[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th']))

plotCI(bp - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim =ylims, xaxt="n", xlim=c(0, 4.5))

abline(h=0)

u <- par('usr')

plotCI(bp - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=0.8, pch=25, pt.bg="white", add=TRUE)
plotCI(bp + jit, delistedOut$RB$S50[, 'median'], li = delistedOut$RB$S50[, '25th'], ui=delistedOut$RB$S50[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=0.8, pch=24, pt.bg="white")
plotCI(bp, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=0.8, pch=21, pt.bg="white")

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)

axis(side=1, at=bp, labels=obs_bias2)

mtext(side=3, line=0.5, adj=0, "d)")

# legend("topright", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[25])), expression(italic(S[AVG])), expression(paste(italic(S[50])))), pt.cex = 0.8, xpd=NA, bty="n")

mtext(side=1, outer=TRUE, "Observation bias post change", line=0)


#####################################################
####################################################
# Supplemental
indCol <- c(ind = "#475A83", nonInd = "#C2B642")
indCol2 <- c(ind = "#475A8380", nonInd = "#C2B64280")

# Ricker parameter calculations
 # What does the distribution look like at low productivity given re-draw crtieria?
y <- rnorm(10^6, quantile(datLoc$prod, 0.025), sd(datLoc$prod))
while(length(which(y<0.4)) > 0){
	y[which(y<0.4)] <- rnorm(length(which(y<0.4)), quantile(datLoc$prod, 0.025), sd(datLoc$prod))
}
dy <- density(y)

x <- seq(0, 3, 0.1)
dx <- dnorm(x, mean(datLoc$prod), sd(datLoc$prod))

quartz(width = 6.3, height = 2.8, pointsize=10)
par(mfrow=c(1,2), mar=c(4,4,2,1))
# Histogram of productivity over all rivers
hist(datLoc$prod, breaks = seq(0, 3, 0.1), main="", xlab=expression(paste("Productivity (", italic(a), ")", sep="")), las=1, col=grey(0.8), border=grey(0.8))
abline(v = mean(datLoc$prod), col=statusCols['g'], lwd=2)
abline(v = quantile(datLoc$prod, 0.025), col=statusCols['r'], lwd=2 )
abline(v = 0.4, lty=3)
lines(x, dx/max(dx)*20, col=statusCols['g'], lwd=1.5)
lines(dy$x, dy$y/max(dy$y)*20, col=statusCols['r'], lwd=1.5)

mtext(side=3, line=0.5, adj=0, "a) Productivity")


# Histogram of log Smax = log(1/b)
hist(log(1/datLoc$densDep), main="", xlab=expression(paste("log ", S[MAX], " (", log(1/italic(b)), ")", sep="")), las=1, , breaks=seq(3, 13, 0.5), col=grey(0.8), border=grey(0.8))
hist(log(1/datLoc$densDep[datLoc$indicator=="Y"]), col=indCol2['ind'], border=indCol['ind'], breaks=seq(3, 13, 0.5), add=TRUE)
hist(log(1/datLoc$densDep[datLoc$indicator=="N"]), col=indCol2['nonInd'], border=indCol['nonInd'], add=TRUE, breaks=seq(3, 13, 0.5))
abline(v = c(mean(log(1/datLoc$densDep[datLoc$indicator=="Y"])), mean(log(1/datLoc$densDep[datLoc$indicator=="N"]))), lwd=1.5, col=indCol)

# log Smax for indicator and non-indicator
avgSmax <- c(mean(log(1/datLoc$densDep[datLoc$indicator=="Y"])), mean(log(1/datLoc$densDep[datLoc$indicator=="N"])))
sdSmax <- c(sd(log(1/datLoc$densDep[datLoc$indicator=="Y"])), sd(log(1/datLoc$densDep[datLoc$indicator=="N"])))

x <- seq(0, 14, 0.1)
lines(x, dnorm(x, avgSmax[1], sdSmax[1])*length(which(datLoc$indicator=="Y")), col=indCol['ind'], lwd=1.5)
lines(x, dnorm(x, avgSmax[2], sdSmax[2])*length(which(datLoc$indicator=="N")), col=indCol['nonInd'], lwd=1.5)
# legend(8, 35, fill=indCol, c("indicator", "non-indicator"), bty="n", xpd=NA)
text(11, 30, "indicator", font=2, col=indCol['ind'])
text(11, 27, "non-indicator", font=2, col=indCol['nonInd'])

mtext(side=3, line=0.5, adj=0, "b) Capacity")
