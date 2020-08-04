###############################################################################
# This code plots the specific figures presented in 
# Peacock et al. 2020, Canadian Journal of Fisheries and Aquatic Sciences
#
# Evaluating the consequences of common assumptions in run reconstructions 
# on Pacific-salmon biological status assessments
#
# Corresponding author: Stephanie J. Peacock <stephanie.j.peacock@gmail.com>
###############################################################################

library(gplots)
statusCols <- c(g = "#8EB687", a = "#DFD98D", r = "#B66A64")

###############################################################################
# Base case simulations
###############################################################################
library(mvtnorm)
library(here)
library(corpcor) # make.positive.definite function for Sigma
library(data.table) # for shift function
library(gsl) # for lambert_W0 function
library(doParallel) # for parallelizing function application over parameters or MCMC iterations


source("model/populationSubmodFns.R")
source("model/obsSubmodFns.R")
source("model/expansionFactors.R")
source("model/benchmarkFns.R")
source("model/performanceFns.R")
source("model/reconstrSimulator.R")
source("runSims/runSensitivity.R")

# Load base parameter values
simPar_all <- read.csv(here("data/baseSimPar.csv"), stringsAsFactors = F)

simPar_base <- list(simPar_all[simPar_all$scenario == "baseGreen",], simPar_all[simPar_all$scenario == "baseAmber",], simPar_all[simPar_all$scenario == "baseRed",])

# for(i in 1:3) simPar_base[[i]]$correlPop <- 0.2

set.seed(8293) #original seed
# set.seed(9823)
out_base <- runSensitivity(parList = simPar_base, nSim = 4000, nCores = 3)[[1]]
# saveRDS(object = out_base, file = "workspaces/out_base.rds")
# readRDS(file = "workspaces/out_base.rds")
# Figure 6: cases2plot <- c(1, 3)
# Figure S7: cases2plot <- c(2)

if(length(cases2plot) == 1){
	# quartz(width = 6.3, height= 3, pointsize = 10)
	pdf(file = "FigS7.pdf", width = 6.3, height= 3, pointsize = 10)
	I <- cases2plot
	layout(matrix(c(2,1,1,2,1,1,7,3,3,5,4,4,5,4,4,8,6,6), nrow=3))
} else if (length(cases2plot) == 2){
	# quartz(width = 6.3, height = 5.5, pointsize=10)
	pdf(file = "Fig6.pdf", width = 6.3, height= 5.5, pointsize = 10)
	I <- cases2plot
	layout(matrix(c(rep(c(2,1,1,8,7,7), 2), 13,3,3,15,9,9, rep(c(5,4,4,11,10,10),2), 14,6,6,16,12,12), nrow=6))
} else {
	quartz(width = 6, height= 7, pointsize = 9)
	par(oma=c(0,3.5,2,0))
	I <- cases2plot
	layout(matrix(c(rep(c(2,1,1,8,7,7,14,13,13), 2), 19,3,3,21,9,9,23,15,15, rep(c(5,4,4,11,10,10,17,16,16),2), 20,6,6,22,12,12,24,18,18), nrow=9))
}
 par(oma=c(0,4,2,0))

 lc <- 0
for(i in I){
for(k in 1:2){
	lc <- lc+1
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
	
	if(k==1){ mtext(side=2, line = 6, c("High productivity\nHarvest Control Rule", "Low productivity\nModerate harvest", "Low productivity\nHigh harvest")[i], cex=0.8, font=2)
	}
	
	par(mar=c(0,5,4,0))
	TR <- apply(statusDiffMat, 2, sum)/length(statusDiff)
	bp <- barplot(TR, col=statusCols[c('r', 'a', 'g')], yaxt="n", names.arg = FALSE)
	text(bp, TR, pos=3, paste(formatC(round(TR*100, 1), 1, format="f"), "%"), xpd=NA)
	
	if(length(cases2plot) == 1){
		mtext(side = 3, adj = 0, line=2, c("a) Spawner-recruitment", "b) Percentile")[k], cex=0.8)
	} else {
		mtext(side = 3, adj = 0, line=2, paste(letters[lc], ")", sep=""), cex=0.8)
		if(i==1) mtext(side=3, line=3.5, c("Spawner-recruitment", "Percentile")[k], cex=0.8, font = 2)
	}
	
	par(mar=c(5,0,0,4))
	OBS <- apply(statusDiffMat, 1, sum)/length(statusDiff)
	bp <- barplot(rev(OBS), col=rev(statusCols[c('r', 'a', 'g')]), yaxt="n", horiz=TRUE, xaxt="n", names.arg=FALSE)
	text(rev(OBS), bp, pos=4, paste(formatC(round(rev(OBS)*100, 1), 1, format="f"), "%"), xpd=NA)
	
	} # end k
	
} # end i

 dev.off()

#------------------------------------------------------------------------------
# Why misclassifications? RB in benchmarks
#------------------------------------------------------------------------------

delistedOut <- delistSensitivity(out_base)

cases2plot <- c(1:3)

if(length(cases2plot) == 2) pdf(file = "Fig7.pdf", width = 6, height= 3, pointsize = 9)
if(length(cases2plot) == 3) pdf(file = "FigS8.pdf", width = 6, height= 3, pointsize = 9)
par(mfrow=c(1,2), mar=c(3, 4, 1, 1), oma=c(2,3.5,2,2))

# RBs

if(length(cases2plot) == 3){
	x <- c(1:3)
	jit <- 0.1
} else if(length(cases2plot) == 2){
	x <- c(1.25, 2.75)
	jit <- 0.15
}
	ptCex <- 1
	

plotCI(x, delistedOut$RB$Sgen1[cases2plot, 'median'], "n", li = delistedOut$RB$Sgen1[cases2plot, '25th'], ui= delistedOut$RB$Sgen1[cases2plot, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", xlim=c(0.5, 3.5))

axis(side=1, at=x, labels=FALSE)
u <- par('usr')
if(length(cases2plot)==3){
	text(x, rep(u[3] - 0.1*(u[4]-u[3]), 3), c("high", "low", "low ")[cases2plot], xpd=NA)
	text(x, rep(u[3] - 0.2*(u[4]-u[3]), 3), c("HCR", "moderate", "high")[cases2plot], xpd=NA)
	text(0, u[3] - 0.1*(u[4]-u[3]), "productivity:", xpd=NA, adj=1, font=2)
	text(0, u[3] - 0.2*(u[4]-u[3]), "target harvest:", xpd=NA, adj=1, font=2)
	
} else {
	text(x, rep(u[3] - 0.1*(u[4]-u[3]), 3), c("high prod.", "low prod.", "low prod.")[cases2plot], xpd=NA)
	text(x, rep(u[3] - 0.2*(u[4]-u[3]), 3), c("HCR", "moderate harvest", "high harvest")[cases2plot], xpd=NA)
}

abline(h=0)

plotCI(x - jit, delistedOut$RB$Sgen1[cases2plot, 'median'], "n", li = delistedOut$RB$Sgen1[cases2plot, '25th'], ui= delistedOut$RB$Sgen1[cases2plot, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)

plotCI(x + jit, delistedOut$RB$Smsy[cases2plot, 'median'], li = delistedOut$RB$Smsy[cases2plot, '25th'], ui=delistedOut$RB$Smsy[cases2plot, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])

plotCI(x, delistedOut$RB$avgS[cases2plot, 'median'], li = delistedOut$RB$avgS[cases2plot, '25th'], ui= delistedOut$RB$avgS[cases2plot, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Relative bias")

legend("topright", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), bty="n", pt.cex = ptCex, pt.bg = c(statusCols['r'], 1, statusCols['g']))

mtext(side=3, line=0.5, adj=0, "a) Spawner-recruitment")
mtext(side=1, "Case", line=3.5)

# Percentile

yAll <- c(delistedOut$RB$S25[cases2plot, '25th'], delistedOut$RB$S25[cases2plot, '75th'], delistedOut$RB$S55[cases2plot, '25th'], delistedOut$RB$S50[cases2plot, '75th'], delistedOut$RB$avgS[cases2plot, '25th'], delistedOut$RB$avgS[cases2plot, '75th'])

plotCI(x, delistedOut$RB$S25[cases2plot, 'median'], "n", li = delistedOut$RB$S25[cases2plot, '25th'], ui= delistedOut$RB$S25[cases2plot, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", ylim=range(yAll), xlim=c(0.5, 3.5))

axis(side=1, at=x, labels=FALSE)
u <- par('usr')
if(length(cases2plot)==3){
	text(x, rep(u[3] - 0.1*(u[4]-u[3]), 3), c("high", "low", "low ")[cases2plot], xpd=NA)
	text(x, rep(u[3] - 0.2*(u[4]-u[3]), 3), c("HCR", "moderate", "high")[cases2plot], xpd=NA)
	
	} else {
		text(x, rep(u[3] - 0.1*(u[4]-u[3]), 3), c("high prod.", "low prod.", "low prod.")[cases2plot], xpd=NA)
		text(x, rep(u[3] - 0.2*(u[4]-u[3]), 3), c("HCR", "moderate harvest", "high harvest")[cases2plot], xpd=NA)
	}
abline(h=0)

plotCI(x - jit, delistedOut$RB$S25[cases2plot, 'median'], "n", li = delistedOut$RB$S25[cases2plot, '25th'], ui= delistedOut$RB$S25[cases2plot, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)

plotCI(x + jit, delistedOut$RB$S50[cases2plot, 'median'], li = delistedOut$RB$S50[cases2plot, '25th'], ui=delistedOut$RB$S50[cases2plot, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])

plotCI(x, delistedOut$RB$avgS[cases2plot, 'median'], li = delistedOut$RB$avgS[cases2plot, '25th'], ui= delistedOut$RB$avgS[cases2plot, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)

mtext(side=3, line=0.5, adj=0, "b) Percentile")
mtext(side=1, "Case", line=3.5)

legend("topright", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[25])), expression(italic(S[AVG])), expression(paste(italic(S[50])))), bty="n", pt.cex = ptCex, pt.bg = c(statusCols['r'], 1, statusCols['g']))

dev.off()

#------------------------------------------------------------------------------
# Comparing to 100% monitoring coverage of indicator, non-indicator, and both
#------------------------------------------------------------------------------

sensitivityPar <- c(1:4)
delistedOut.all <- list(); length(delistedOut.all) <- 3
delistedOut.all[[1]] <- readRDS("workspaces/base100Mon_delisted_baseGreen.rds")
delistedOut.all[[2]] <- readRDS("workspaces/base100Mon_delisted_baseAmber.rds")
delistedOut.all[[3]] <- readRDS("workspaces/base100Mon_delisted_baseRed.rds")

pdf(file = "FigS9.pdf", width = 6, height= 6.5, pointsize = 10)
# quartz(width = 6, height= 6.5, pointsize = 10)
# layout(matrix(c(1,2,3,3), nrow=2, byrow=2))
par(mfrow=c(3,2), mar=c(3, 4, 1, 1), oma=c(0,6,5,1))

x <- c(1:4)
jit <- 0.1
ptCex <- 1

for(i in 1:3){
	delistedOut <- delistedOut.all[[i]]
	
	# Spawner-recruitment
	yAll <- c(delistedOut$RB$Sgen1[, '25th'], delistedOut$RB$Sgen1[, '75th'], delistedOut$RB$Smsy[, '25th'], delistedOut$RB$Smsy[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th'])


	plotCI(x - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", xlim=c(0.5, 4.5), ylim=range(yAll))

	u <- par('usr')
	axis(side=1, at=x, labels=FALSE)
	text(x, u[3] - 0.15*(u[4]-u[3]), c("Default", "100%\nind.", "100%\nnon-ind", "100%\nboth"), xpd=NA)

	abline(h=0)
	
	plotCI(x - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)
	
	plotCI(x + jit, delistedOut$RB$Smsy[, 'median'], li = delistedOut$RB$Smsy[, '25th'], ui=delistedOut$RB$Smsy[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])
	
	plotCI(x, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)
	
	A <- axis(side=2, labels = FALSE)
	axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
	mtext(side=2, line=4, "Relative bias")
	
	if(i==1) legend(3, u[4] + 0.22*(u[4]-u[3]), pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), pt.cex = ptCex, pt.bg=c(statusCols['r'], 1, statusCols['g']), xpd=NA, bg="white")
	
	if(i==1) mtext(side=3, line=4, font = 2, "Spawner-recruitment")
	mtext(side = 2, line = 6, c("High-productivity\nHCR", "Low-productivity\nModerate-harvest", "Low-productivity\nHigh-harvest")[i], font = 2)
	
	# Percentile
	yAll <- c(delistedOut$RB$S25[, '25th'], delistedOut$RB$S25[, '75th'], delistedOut$RB$S50[, '25th'], delistedOut$RB$S50[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th'])
	
	plotCI(x - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", ylim=range(yAll), xlim=c(0.5, 4.5))
	
	u <- par('usr')
	axis(side=1, at=x, labels=FALSE)
	text(x, u[3] - 0.15*(u[4]-u[3]), c("Default", "100%\nind.", "100%\nnon-ind", "100%\nboth"), xpd=NA)
	
	
	abline(h=0)
	
	plotCI(x - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)
	
	plotCI(x + jit, delistedOut$RB$S50[, 'median'], li = delistedOut$RB$S50[, '25th'], ui=delistedOut$RB$S50[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])
	
	plotCI(x, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)
	
	A <- axis(side=2, labels = FALSE)
	axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
	
	if(i==1) mtext(side=3, line=4, font = 2, "Percentile")
	
	if(i==1) legend(4.1, u[4] + 0.22*(u[4]-u[3]), pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[S25])), expression(italic(S[AVG])), expression(paste(italic(S[50])))), pt.cex = 0.8, pt.bg=c(statusCols['r'], 1, statusCols['g']), xpd=NA, bg="white")

} # end i

dev.off()

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
baseValue <- 1
ptCex <- 1

for(baseCaseNum in 1:3){
	
if(baseCaseNum == 1){
	delistedOut <- readRDS("workspaces/nPop_delisted_baseGreen.rds")
	pdf(file = "FigS10.pdf", width = 6, height= 4.5, pointsize = 10)
} else if(baseCaseNum == 2){
	delistedOut <- readRDS("workspaces/nPop_delisted_baseAmber.rds")
	pdf(file = "FigS11.pdf", width = 6, height= 4.5, pointsize = 10)
} else if(baseCaseNum == 3){
	delistedOut <- readRDS("workspaces/nPop_delisted_baseRed.rds")
	pdf(file = "FigS12.pdf", width = 6, height= 4.5, pointsize = 10)
}
# quartz(width = 6, height= 4.5, pointsize = 10)
par(mfrow=c(2,2), mar=c(3, 4, 1, 1), oma=c(3,2,2,0))

# top: proportion wrong
for(j in 1:2){ # for SR and HS (percentile)
	bp <- barplot(t(delistedOut$ppn[[j]]), 
								col=c(statusCols, grey(0.8), 1), las=1, ylab = "", 
								space=c(0.5, 0.5, 0.1, 0.5, 0.1),
								border=NA, yaxs="i")
	abline(h = 0)
	
	# text(bp, 0.02, srt=90, rep(c("Base", "Decline"), 4), adj=0, col="white")
	
	text(bp, rep(-0.1, 5), LETTERS[1:5], xpd=NA)
	
	points(bp[1], 1.05, pch=8, xpd=NA) # Base case
	
	if(j == 1) mtext(side=2, line = 4, "Proportion of simulations")
	
	mtext(side=3, line=2, c("Spawner-recruitment", "Percentile")[j])
	mtext(side=3, line=0.5, adj=0, c("a)", "b)")[j])
}

uBP <- par('usr')

# RBs

x <- bp
jit <- 0.1

# Spawner-recruitment

yAll <- c(delistedOut$RB$Sgen1[, '25th'], delistedOut$RB$Sgen1[, '75th'], delistedOut$RB$Smsy[, '25th'], delistedOut$RB$Smsy[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th'])


plotCI(x - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", xlim=uBP[1:2], ylim=range(yAll))

axis(side=1, at=x, LETTERS[1:5])
u <- par('usr')
polygon(x=x[1] + c(-2, 2, 2, -2)*jit, y = rep(u[3:4], each=2), col="#00000020", border=NA)
abline(h=0)

plotCI(x - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=c("white", rep(c("white", statusCols['r']), 2)), add=TRUE)

plotCI(x + jit, delistedOut$RB$Smsy[, 'median'], li = delistedOut$RB$Smsy[, '25th'], ui=delistedOut$RB$Smsy[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=c("white", rep(c("white", statusCols['g']), 2)))

plotCI(x, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=c("white", rep(c("white", 1), 2)))

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
mtext(side=2, line=4, "Relative bias")

# legend("topright", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), bty="n", pt.cex = 0.8)

mtext(side=3, line=0.5, adj=0, "c)")

u <- par('usr')
y <- u[3] - 0.25*(u[4] - u[3])
segments(x0=x[1]-jit, x1=x[1]+jit, y0=y, y1=y, xpd=NA); text(x[1], y, pos=1, "Default", xpd=NA)
segments(x0=x[2] - jit, x1=x[3] + jit, y0=y, y1=y, xpd=NA); text(mean(x[2:3]), y, pos=1, "Small CUs", xpd=NA)
segments(x0=x[4] - jit, x1=x[5] + jit, y0=y, y1=y, xpd=NA); text(mean(x[4:5]), y, pos=1, "Large CUs", xpd=NA)

# Percentile

yAll <- c(delistedOut$RB$S25[, '25th'], delistedOut$RB$S25[, '75th'], delistedOut$RB$S50[, '25th'], delistedOut$RB$S50[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th'])

plotCI(x - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", ylim=range(yAll), xlim=uBP[1:2])

axis(side=1, at=x, LETTERS[1:5])
u <- par('usr')
polygon(x=x[1] + c(-2, 2, 2, -2)*jit, y = rep(u[3:4], each=2), col="#00000020", border=NA)
abline(h=0)

plotCI(x - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=c("white", rep(c("white", statusCols['r']), 2)), add=TRUE)

plotCI(x + jit, delistedOut$RB$S50[, 'median'], li = delistedOut$RB$S50[, '25th'], ui=delistedOut$RB$S50[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=c("white", rep(c("white", statusCols['g']), 2)))

plotCI(x, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=c("white", rep(c("white", 1), 2)))

A <- axis(side=2, labels = FALSE)
axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)

mtext(side=3, line=0.5, adj=0, "d)")

u <- par('usr')
y <- u[3] - 0.25*(u[4] - u[3])
segments(x0=x[1]-jit, x1=x[1]+jit, y0=y, y1=y, xpd=NA); text(x[1], y, pos=1, "Default", xpd=NA)
segments(x0=x[2] - jit, x1=x[3] + jit, y0=y, y1=y, xpd=NA); text(mean(x[2:3]), y, pos=1, "Small CUs", xpd=NA)
segments(x0=x[4] - jit, x1=x[5] + jit, y0=y, y1=y, xpd=NA); text(mean(x[4:5]), y, pos=1, "Large CUs", xpd=NA)

mtext(side=1, outer=TRUE, "Scenario for number of indicator & non-indicator streams", line=1)

dev.off()

} # end all base cases



###############################################################################
# Capacity X Real monitoring scenarios
###############################################################################

plot_capacityXmon <- function(){
	par(mfrow=c(2,2), mar=c(4, 4, 1, 1), oma=c(1,2,3,0))
	
	# Ppn Wrong
	
	
	for(j in 1:2){
		
		bp <- barplot(t(delistedOut$ppn[[j]][includeIndices, ]), col=c(statusCols, grey(0.8), 1), las=1, ylab = "", border=NA, xaxs="i", yaxs="i", space=rep(c(0.6, rep(0.1, 3)), length(scenarios)))
		
		abline(h=0)
		
		uBP <- par('usr')
		monLab <- c(mean(bp[2:3]), mean(bp[6:7]), mean(bp[10:11]), mean(bp[14:15]))
		# axis(side = 3, at = c(mean(bp[4:5]), mean(bp[8:9]), mean(bp[12:13])), labels=FALSE)
		text(monLab, rep(uBP[4], 4), pos=3, xpd=NA, monLabels)
		
		# axis(side=1, at=bp, labels=FALSE, tck=-0.02)
		axis(side=1, at=bp, labels=rep(redHab, length(scenarios)), tck=-0.02, las=lasLab, mgp=c(3, 0.5, -0.1))
		
		if(j == 1) mtext(side=2, line = 4, "Proportion of simulations")
		mtext(side=3, line=3, c("Spawner-recruitment", "Percentile")[j])
		mtext(side=3, line=1.5, adj=0, c("a)", "b)")[j])
	}
	
	# RBs
	
	ptCex <- 1
	
	ylims <- range(c(delistedOut$RB$Sgen1[includeIndices, '25th'], delistedOut$RB$Sgen1[includeIndices, '75th'], delistedOut$RB$Smsy[includeIndices, '25th'], delistedOut$RB$Smsy[includeIndices, '75th'], delistedOut$RB$avgS[includeIndices, '25th'], delistedOut$RB$avgS[includeIndices, '75th']))
	# ylims <- c(-0.2, 1.5)
	plotCI(bp - jit, delistedOut$RB$Sgen1[includeIndices, 'median'], "n", li = delistedOut$RB$Sgen1[includeIndices, '25th'], ui= delistedOut$RB$Sgen1[includeIndices, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim = ylims, xaxt="n", xlim=range(c(bp-jit, bp+jit)))
	
	u <- par('usr')
	# mtext(side=3, "Monitoring scenario", line=1.5)
	axis(side = 3, at = c(mean(bp[4:5]), mean(bp[8:9]), mean(bp[12:13])), labels=FALSE)
	text(monLab, rep(u[4], 4), pos=3, xpd=NA, monLabels)
	
	abline(h=0)
	
	plotCI(bp - jit, delistedOut$RB$Sgen1[includeIndices, 'median'], "n", li = delistedOut$RB$Sgen1[includeIndices, '25th'], ui= delistedOut$RB$Sgen1[includeIndices, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)
	
	plotCI(bp + jit, delistedOut$RB$Smsy[includeIndices, 'median'], li = delistedOut$RB$Smsy[includeIndices, '25th'], ui=delistedOut$RB$Smsy[includeIndices, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])
	
	plotCI(bp, delistedOut$RB$avgS[includeIndices, 'median'], li = delistedOut$RB$avgS[includeIndices, '25th'], ui= delistedOut$RB$avgS[includeIndices, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)
	
	abline(v=c(mean(bp[4:5]), mean(bp[8:9]), mean(bp[12:13])))
	
	A <- axis(side=2, labels = FALSE)
	axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
	mtext(side=2, line=4, "Relative bias")
	
	axis(side=1, at=bp, labels=rep(redHab, length(scenarios)), tck=-0.02, mgp=c(3, 0.5, 0), las=lasLab)
	
	mtext(side=3, line=1.5, adj=0, "c)")
	
	# legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), pt.cex = ptCex, xpd=NA, ncol=3, bg="white", pt.bg = c(statusCols['r'], 1, statusCols['g']))
	
	# Percentile
	
	ylims <- range(c(delistedOut$RB$S25[includeIndices, '25th'], delistedOut$RB$S25[includeIndices, '75th'], delistedOut$RB$S50[includeIndices, '25th'], delistedOut$RB$S50[includeIndices, '75th'], delistedOut$RB$avgS[includeIndices, '25th'], delistedOut$RB$avgS[includeIndices, '75th']))
	# ylims <- c(-0.2, 4.5)
	plotCI(bp - jit, delistedOut$RB$S25[includeIndices, 'median'], "n", li = delistedOut$RB$S25[includeIndices, '25th'], ui= delistedOut$RB$S25[includeIndices, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim =ylims, xaxt="n", xlim=range(c(bp-jit, bp+jit)))
	
	abline(h=0)
	
	u <- par('usr')
	# mtext(side=3, "Monitoring scenario", line=1.5)
	axis(side = 3, at = c(mean(bp[4:5]), mean(bp[8:9]), mean(bp[12:13])), labels=FALSE)
	text(monLab, rep(u[4], 4), pos=3, xpd=NA, monLabels)
	
	
	plotCI(bp - jit, delistedOut$RB$S25[includeIndices, 'median'], "n", li = delistedOut$RB$S25[includeIndices, '25th'], ui= delistedOut$RB$S25[includeIndices, '75th'], gap=0, sfrac = 0, col=statusCols['r'], pch=25, pt.bg=statusCols['r'], add=TRUE, cex = ptCex)
	plotCI(bp + jit, delistedOut$RB$S50[includeIndices, 'median'], li = delistedOut$RB$S50[includeIndices, '25th'], ui=delistedOut$RB$S50[includeIndices, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])
	plotCI(bp, delistedOut$RB$avgS[includeIndices, 'median'], li = delistedOut$RB$avgS[includeIndices, '25th'], ui= delistedOut$RB$avgS[includeIndices, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)
	
	abline(v=c(mean(bp[4:5]), mean(bp[8:9]), mean(bp[12:13])))
	
	A <- axis(side=2, labels = FALSE)
	axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
	
	axis(side=1, at=bp, labels=rep(redHab, length(scenarios)), tck=-0.02, las=lasLab, mgp=c(3, 0.5, 0))
	
	mtext(side=3, line=1.5, adj=0, "d)")
	
	# legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[25])), expression(italic(S[AVG])), expression(paste(italic(S[50])))), pt.cex = ptCex, xpd=NA, ncol = 3, bg="white", pt.bg = c(statusCols['r'], 1, statusCols['g']))
	
	mtext(side=1, outer=TRUE, "Percentage of subpopulations with severe declines in capacity", line=-1)

	} # end function plot_capacityXmon

#------------------------------------------------------------------------------

redHab <- c(seq(0, 50, 25), 100)
declCap <- redHab
jit <- 0.2

for(baseCaseNum in 1:3){
	if(baseCaseNum == 1){
		delistedOut <- readRDS("workspaces/capacityXmon_delisted_baseGreen.rds")
		
		# Include all four monitoring scenarios or just the first two?
		for(s in 1:2){
			if(s == 1){
				scenarios <- c(1:2)
			} else {
				scenarios <- c(1:4)
			}
		
			includeIndices <- ((scenarios[1]-1)*4+1):((scenarios[1]-1)*4+4)
			
			if(length(scenarios) > 1){
				for(i in 2:length(scenarios)){
					includeIndices <- c(includeIndices, ((scenarios[i]-1)*4+1):((scenarios[i]-1)*4+4))
				}
			}
			
			if(length(scenarios) > 2){
				monLabels <- LETTERS[1:4]
				lasLab <- 2 
			} else{
				monLabels <- c("No change in monitoring", "Observed decline in monitoring", "Chum", "Recent")
				lasLab <- 1
			}
			
			if(length(scenarios) > 2){
			pdf(file = "FigS13.pdf", width = 8, height= 5, pointsize = 10)
			plot_capacityXmon()
			dev.off()
		} else {
			pdf(file = "Fig8.pdf", width = 8, height = 5, pointsize = 10)
			plot_capacityXmon()
			dev.off()
		}
		} # end plot2 vs plot4
	# end baseCase 1
	}	else if(baseCaseNum == 2){
		delistedOut <- readRDS("workspaces/capacityXmon_delisted_baseAmber.rds")
		pdf(file = "FigS14.pdf", width = 8, height= 5, pointsize = 10)
		plot_capacityXmon()
		dev.off()
	} else if(baseCaseNum == 3){
		delistedOut <- readRDS("workspaces/capacityXmon_delisted_baseRed.rds")
		pdf(file = "FigS15.pdf", width = 8, height= 5, pointsize = 10)
		plot_capacityXmon()
		dev.off()
	}

} # end for baseCaseNum in 1:3

#------------------------------------------------------------------------------
# Different correlation among subpopulations
scenarios <- c(1:2)
includeIndices <- ((scenarios[1]-1)*4+1):((scenarios[1]-1)*4+4)
for(i in 2:length(scenarios)){
		includeIndices <- c(includeIndices, ((scenarios[i]-1)*4+1):((scenarios[i]-1)*4+4))
}
monLabels <- c("No change in monitoring", "Observed decline in monitoring", "Chum", "Recent")
lasLab <- 1

delistedOut <- readRDS("workspaces/capacityXmon_delisted_baseGreen_correlPop0.rds")
pdf(file = "FigS16.pdf", width = 8, height= 5, pointsize = 10)
plot_capacityXmon()
dev.off()

delistedOut <- readRDS("workspaces/capacityXmon_delisted_baseGreen_correlPop9.rds")
pdf(file = "FigS17.pdf", width = 8, height= 5, pointsize = 10)
plot_capacityXmon()
dev.off()

###############################################################################
# Observation bias
###############################################################################

plot_obsBias <- function(){
	par(mfrow=c(2,2), mar=c(3, 4, 1, 1), oma=c(1,2,2,0))

	# Ppn Wrong
	
	for(j in 1:2){
		
		bp <- barplot(t(delistedOut$ppn[[j]]), col=c(statusCols, grey(0.8), 1), las=1, ylab = "", space=0.1, border=NA, xaxs="i", yaxs="i")
		axis(side = 1 , at = predict(lm(bp ~ sensitivityPar), newdata = data.frame(sensitivityPar = seq(-1.6, 0, 0.4))), labels=formatC(seq(-1.6, 0, 0.4), 1, format="f"))
		
		points(predict(lm(bp ~ sensitivityPar), newdata = data.frame(sensitivityPar = baseValue)), 1.05, pch=8, xpd=NA)
		if(j == 1) mtext(side=2, line = 4, "Proportion of simulations")
		mtext(side=3, line=2, c("Spawner-recruitment", "Percentile")[j])
		mtext(side=3, line=0.5, adj=0, c("a)", "b)")[j])
		}
	
	
	# RBs
	ptCex <- 1
	jit <- 0.02
	
	# Spawner-recruitment
	yAll <- range(c(delistedOut$RB$Sgen1[, '75th'], delistedOut$RB$Sgen1[, '25th'], delistedOut$RB$Smsy[, '25th'], delistedOut$RB$Smsy[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th']))
	
	plotCI(sensitivityPar - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xlim=c(-1.7, 0.1), xaxt="n", ylim=yAll)
	axis(side=1, at=seq(-1.6, 0, 0.4))
	
	abline(v = baseValue, lwd=16, col=grey(0.8))
	abline(h=0)
	
	plotCI(sensitivityPar - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)
	
	plotCI(sensitivityPar + jit, delistedOut$RB$Smsy[, 'median'], li = delistedOut$RB$Smsy[, '25th'], ui=delistedOut$RB$Smsy[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])
	
	plotCI(sensitivityPar, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)
	
	A <- axis(side=2, labels = FALSE)
	axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
	mtext(side=2, line=4, "Relative bias")
	
	mtext(side=3, line=0.5, adj=0, "c)")
	legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), bty="n", pt.cex = ptCex, pt.bg = c(statusCols['r'], 1, statusCols['g']))
	
	# Percentile
	
	yAll <- range(c(delistedOut$RB$S25[, '75th'], delistedOut$RB$S25[, '25th'], delistedOut$RB$S50[, '25th'], delistedOut$RB$S50[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th']))
	
	plotCI(sensitivityPar - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], yaxt="n", col=NA, ylab="", xlab="", ylim=yAll, xlim=c(-1.7, 0.1), xaxt="n")
	axis(side=1, at=seq(-1.6, 0, 0.4))
	
	abline(v = baseValue, lwd=16, col=grey(0.8))
	abline(h=0)
	
	plotCI(sensitivityPar - jit, delistedOut$RB$S25[, 'median'], "n", li = delistedOut$RB$S25[, '25th'], ui= delistedOut$RB$S25[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)
	plotCI(sensitivityPar + jit, delistedOut$RB$S50[, 'median'], li = delistedOut$RB$S50[, '25th'], ui=delistedOut$RB$S50[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])
	plotCI(sensitivityPar, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)
	
	A <- axis(side=2, labels = FALSE)
	axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
	
	mtext(side=3, line=0.5, adj=0, "d)")
	legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[25])), expression(italic(S[AVG])), expression(italic(S[50]))), bty="n", pt.cex = ptCex, pt.bg = c(statusCols['r'], 1, statusCols['g']))
	
	mtext(side=1, outer=TRUE, expression(paste("Observation bias of spawners (", bar(delta), ")")), line=0)

} # end function plot_obsBias

#------------------------------------------------------------------------------

obsBias <- seq(-1.6, 0, 0.2)
sensitivityPar <- obsBias
baseValue <- -0.4
# quartz(width = 6, height= 4.5, pointsize = 10)

for(baseCaseNum in 1:3){
	if(baseCaseNum == 1){
		delistedOut <- readRDS("workspaces/obsBias_delisted_baseGreen.rds")
		pdf(file = "Fig9.pdf", width = 6, height= 4.5, pointsize = 10)
		plot_obsBias()
		dev.off()
	} else if(baseCaseNum == 2){
		delistedOut <- readRDS("workspaces/obsBias_delisted_baseAmber.rds")
		pdf(file = "FigS18.pdf", width = 6, height= 4.5, pointsize = 10)
		plot_obsBias()
		dev.off()
	} else if(baseCaseNum == 3){
		delistedOut <- readRDS("workspaces/obsBias_delisted_baseRed.rds")
		pdf(file = "FigS19.pdf", width = 6, height= 4.5, pointsize = 10)
		plot_obsBias()
		dev.off()
	}
		
} # end baseCaseNum	

###############################################################################
# Catch bias
###############################################################################


#------------------------------------------------------------------------------
# Don't need to include Percentile here because catch bias does not 
# affect spawner estimates, only SR relationship.

plot_catchBias <- function(){
	par(mfrow=c(2,1), mar=c(3, 4, 1, 1), oma=c(1,1,2,0))
# Ppn Wrong


	j <- 1
	# for(j in 1:2){
		
		bp <- barplot(t(delistedOut$ppn[[j]]), col=c(statusCols, grey(0.8), 1), las=1, ylab = "", space=0.1, border=NA, xaxs="i", yaxs="i")
		
		axis(side = 1, at = bp[seq(1, length(bp), 2)], labels=formatC(seq(-1, 1, 0.4), 1, format="f"))
		
		points(bp[findInterval(baseValue, sensitivityPar)], 1.05, pch=8, xpd=NA)
		if(j == 1) mtext(side=2, line = 4, "Proportion of simulations")
		mtext(side=3, line=2, c("Spawner-recruitment", "Percentile")[j])
		mtext(side=3, line=0.5, adj=0, c("a)", "b)")[j])
	# }
		
	
	# RBs
		jit <- 0.015
		
		
		ptCex <- 1
		
		# Spawner-recruitment
		yAll <- range(c(delistedOut$RB$Sgen1[, '75th'], delistedOut$RB$Sgen1[, '25th'], delistedOut$RB$Smsy[, '25th'], delistedOut$RB$Smsy[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th']))
		
		
	plotCI(sensitivityPar, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", ylim=yAll)
	axis(side=1, at = seq(-1, 1, 0.4))
	abline(v = baseValue, lwd=10, col=grey(0.8))
	abline(h=0)
	
	plotCI(sensitivityPar - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)
	plotCI(sensitivityPar + jit, delistedOut$RB$Smsy[, 'median'], li = delistedOut$RB$Smsy[, '25th'], ui=delistedOut$RB$Smsy[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])
	plotCI(sensitivityPar, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)
	
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
	# legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), bty="n", pt.cex = 0.8)
	
	# # Percentile
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
} # end plot_catch Bias 

#------------------------------------------------------------------------------

catchBias <- seq(-1, 1, 0.2)
sensitivityPar <- catchBias
baseValue <- 0
# quartz(width = 3.2, height= 4.5, pointsize = 10)
# layout(matrix(c(1,2,3,3), nrow=2, byrow=2))

for(baseCaseNum in 1:3){
	if(baseCaseNum == 1){
		delistedOut <- readRDS("workspaces/catchBias_delisted_baseGreen.rds")
		pdf(file = "Fig10.pdf", width = 3.2, height= 4.5, pointsize = 10)
		plot_catchBias()
		dev.off()
	} else if(baseCaseNum == 2){
		delistedOut <- readRDS("workspaces/catchBias_delisted_baseAmber.rds")
		pdf(file = "FigS22.pdf", width = 3.2, height= 4.5, pointsize = 10)
		plot_catchBias()
		dev.off()
	} else if(baseCaseNum == 3){
		delistedOut <- readRDS("workspaces/catchBias_delisted_baseRed.rds")
		pdf(file = "FigS23.pdf", width = 3.2, height= 4.5, pointsize = 10)
		plot_catchBias()
		dev.off()
	}
	
} # end baseCaseNum	

#------------------------------------------------------------------------------

# What is the change in misclassification rates over change in bias?
mis <- data.frame(
	TrueStatus = rep(c("green", "amber", "red"), each = length(catchBias)),
	CatchBias = rep(catchBias, 3),
	PercentChange = exp(rep(catchBias, 3)),
	green = rep(NA, 3*length(catchBias)),
	amber = rep(NA, 3*length(catchBias)),
	red = rep(NA, 3*length(catchBias)),
	pessimistic = rep(NA, 3*length(catchBias)),
	optimistic = rep(NA, 3*length(catchBias)),
	totalMis = rep(NA, 3*length(catchBias)),
	propMis = rep(NA, 3*length(catchBias)))

delistedOut <- readRDS("workspaces/catchBias_delisted_baseGreen.rds")
mis[mis$TrueStatus=="green", 4:8] <- delistedOut$ppn$SR
delistedOut <- readRDS("workspaces/catchBias_delisted_baseAmber.rds")
mis[mis$TrueStatus=="amber", 4:8] <- delistedOut$ppn$SR
delistedOut <- readRDS("workspaces/catchBias_delisted_baseRed.rds")
mis[mis$TrueStatus=="red", 4:8] <- delistedOut$ppn$SR

mis$totalMis <- mis$pessimistic + mis$optimistic

for(i in 1:3){
	I <- which(mis$TrueStatus == c("green", "amber", "red")[i])
	mis$propMis[I] <- (mis$totalMis[I] - (mis$totalMis[I])[which(mis$CatchBias[I] == 0)]) / (mis$totalMis[I])[which(mis$CatchBias[I] == 0)]
}

mis$col <- as.numeric(mis$TrueStatus)
mis$col[mis$TrueStatus == "amber"] <- statusCols['a']
mis$col[mis$TrueStatus == "green"] <- statusCols['g']
mis$col[mis$TrueStatus == "red"] <- statusCols['r']

pdf(file = "FigS24.pdf", width = 3.2, height = 4.5, pointsize = 9)
# quartz(width = 3.2, height= 4.5, pointsize = 9)
# layout(matrix(c(1,2,3,3), nrow=2, byrow=2))
par(mfrow=c(2,1), mar=c(3, 5, 2, 1), oma=c(1,0,2,0))

plot(mis$PercentChange, mis$totalMis, "n", ylab="Rate of\nmisclassification", las=1, ylim=c(0.2, 0.6))
for(i in 1:3){
	I <- which(mis$TrueStatus == c("green", "amber", "red")[i])
	points(mis$PercentChange[I], mis$totalMis[I], "o", pch=21, bg="white", col=mis$col[I])
}
axis(side=3, at=mis$PercentChange, labels=mis$CatchBias)
mtext(side = 3, expression(paste("Observation bias in catch (", bar(chi), ")")), line = 2.5)
abline(v = 1)

plot(mis$PercentChange, mis$propMis, "n", col=mis$col, ylab="Change in rate of\nmisclassification", las=1)
for(i in 1:3){
	I <- which(mis$TrueStatus == c("green", "amber", "red")[i])
	points(mis$PercentChange[I], mis$propMis[I], "o", pch=21, bg="white", col=mis$col[I])
}
axis(side=3, at=mis$PercentChange, labels=mis$CatchBias)
abline(v = 1)
abline(h = 0)
text(-0.4, 0.4, "b)", xpd = NA)
text(-0.4, 1.32, "a)", xpd = NA)


mtext(side = 1, expression(paste("Percent observation bias (", italic(e)^{bar(chi)}, ")")), line = 2.5)
abline(h = c(0.1), lty=3)

dev.off()

mis[round(mis$CatchBias, 1) == round(0.4,1), ]
###############################################################################
# Changes in observation bias part-way through time series 
###############################################################################
obs_bias2 <- c(-1.6, -0.7, -0.4,  0)

#------------------------------------------------------------------------------
# Example of scenarios
#------------------------------------------------------------------------------
pdf(file = "FigS20.pdf", width = 4.5, height = 3, pointsize = 10)
par(mar=c(4,5,2,5), oma=c(0,0,0,0))
plot(1:50, rep(-0.4, 50), "l", ylim=c(-1.6, 0), ylab=expression(paste("Observation bias (", bar(delta), ")")), las=1, xaxs="i", lwd=3, col=grey(0.8), xlab="Year in simulation")
mtext(side=4, line=3.5, "Matching Expansion Factor III")
axis(side=4, at=log(1/c(1:5)), labels=1:5, las=1)

lines(c(1:50), rep(c(-0.4, -1.6), each=25))
lines(c(1:50), rep(c(-0.4, -0.7), each=25), lty=2)
lines(c(1:50), rep(c(-0.4, 0), each=25), lty=3)
dev.off()
#------------------------------------------------------------------------------

delistedOut <- readRDS("workspaces/obsBiasChange_delisted_baseGreen.rds")
pdf(file = "FigS21.pdf", width = 6, height = 4.5, pointsize = 10)
# quartz(width = 6, height= 4.5, pointsize = 10)
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
	mtext(side=3, line=2, c("Spawner-recruitment", "Percentile")[j])
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

# legend("topleft", pch = c(25, 21, 24), col=c(statusCols['r'], 1, statusCols['g']), c(expression(italic(S[GEN])), expression(italic(S[AVG])), expression(paste("80%", italic(S[MSY])))), pt.cex = 0.8, xpd=NA, ncol=3, bg="white")

#--------------------------
# Percentile
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

dev.off()
###############################################################################
# Supplemental
###############################################################################

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


###############################################################################
# Inter-annual variability in age-at-return
###############################################################################

# Don't need to include Percentile here because does not 
# affect spawner estimates, only SR relationship.

delistedOut.all <- list(); length(delistedOut.all) <- 3
delistedOut.all[[1]] <- readRDS("workspaces/ageErr_delisted_baseGreen.rds")
delistedOut.all[[2]] <- readRDS("workspaces/ageErr_delisted_baseAmber.rds")
delistedOut.all[[3]] <- readRDS("workspaces/ageErr_delisted_baseRed.rds")
# ageErr <- seq(0.1, 0.9, 0.1)
ageErr <- seq(0.2, 1.6, 0.2)

sensitivityPar <- ageErr
baseValue <- 0.8

#------------------------------------------------------------------------------
# Plot
pdf(file = "FigS25.pdf", width = 6.3, height= 3.6, pointsize = 10)
# quartz(width = 6.3, height= 3.6, pointsize = 10)
par(mfcol=c(2,3), mar=c(3, 4, 1, 1), oma=c(1,1,3,0))

for(s in 1:3){ # For each base case
	delistedOut <- delistedOut.all[[s]]
	
	# Ppn Wrong
	j <- 1
	bp <- barplot(t(delistedOut$ppn[[j]]), col=c(statusCols, grey(0.8), 1), las=1, ylab = "", space=0.1, border=NA, xaxs="i", yaxs="i")
	
	axis(side = 1, at = bp[seq(1, length(bp), 2)], labels=formatC(ageErr[seq(1, length(ageErr),2)], 1, format="f"))
	
	points(bp[findInterval(baseValue, sensitivityPar)], 1.05, pch=8, xpd=NA)
	if(s == 1) mtext(side=2, line = 4, "Proportion of simulations", cex=0.8)
	mtext(side=3, line=0.5, adj=0, c("a)", "c)", "e)")[s], cex=0.8)
	
	mtext(side = 3, line = 1.5, c("High productivity\nHarvest Control Rule", "Low productivity\nModerate harvest", "Low productivity\nHigh harvest")[s], cex=0.8, font=2)
	
	# RBs
	jit <- 0.015
	ptCex <- 1
	
	# Spawner-recruitment
	yAll <- range(c(delistedOut$RB$Sgen1[, '75th'], delistedOut$RB$Sgen1[, '25th'], delistedOut$RB$Smsy[, '25th'], delistedOut$RB$Smsy[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th']))
	
	
	plotCI(sensitivityPar, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", ylim=yAll, xlim=c(0.15, 1.65))
	axis(side=1, at = ageErr[seq(1, length(ageErr),2)])
	abline(v = baseValue, lwd=10, col=grey(0.8))
	abline(h=0)
	
	plotCI(sensitivityPar - jit, delistedOut$RB$Sgen1[, 'median'], "n", li = delistedOut$RB$Sgen1[, '25th'], ui= delistedOut$RB$Sgen1[, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)
	plotCI(sensitivityPar + jit, delistedOut$RB$Smsy[, 'median'], li = delistedOut$RB$Smsy[, '25th'], ui=delistedOut$RB$Smsy[, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])
	plotCI(sensitivityPar, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)
	
	A <- axis(side=2, labels = FALSE, at=c(0, 0.4, 0.8, 1.2))
	axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
	if(s==1) mtext(side=2, line=4, "Relative bias", cex=0.8)
	
	mtext(side=3, line=0.5, adj=0, c("b)", "d)", "f)")[s], cex=0.8) 
	if(s==2) mtext(side=1, expression(paste("Interannual variability in age-at-maturity (", bar(omega), ")")), line=3)
	
}
dev.off()
#-----------
# What are the numbers in the text?
# Increase in mislcassifications from omega = 0.2 to default value of 0.8:
case <- 1
cbind(ageErr, apply(delistedOut.all[[case]]$ppn$SR[,4:5], 1, sum))

#-----------
# How variable is this really?
ppnAge <- list(
	ppnAgeErr(ppnAgeVec = c(0, 0.23, 0.64, 0.13, 0), omega = 0.2, nYears = 1000), 
	ppnAgeErr(ppnAgeVec = c(0, 0.23, 0.64, 0.13, 0), omega = 0.8, nYears = 1000),
	ppnAgeErr(ppnAgeVec = c(0, 0.23, 0.64, 0.13, 0), omega = 1.6, nYears = 1000) 
)

quartz(width = 3.2, height = 5, pointsize = 10)
par(mfrow=c(3,1), mar=c(4,4,2,1))
for(i in 1:3){
	hist(ppnAge[[i]][, 3], xlim=c(0, 1), col="#FF000050", border=NA, main=substitute(paste(bar(omega)==o), list(o = c(0.2, 0.8, 1.6)[i])), las=1, xlab="", breaks=seq(0, 1, 0.02), ylim=c(0, 250), cex.lab=1.2, cex.axis=1.2)
	hist(ppnAge[[i]][, 2], xlim=c(0, 1), col="#00FF0050", border=NA, add=TRUE, breaks=seq(0, 1, 0.02))
	hist(ppnAge[[i]][, 4], xlim=c(0, 1), col="#0000FF50", border=NA, add=TRUE, breaks=seq(0, 1, 0.02))
	abline(v = c(0.23, 0.64, 0.13), col=c(3,2,4))
	if(i==1) legend('topright', lwd=1, col=c(3,2,4), legend=c("Age 3", "Age 4", "Age 5"), bty="n", cex=1.2)
	mtext(side=3, line=0.5, adj=0, paste(letters[i], ")", sep=""))
}
mtext(side=1, "Proportion of recruits", line=3)


###############################################################################
# Variability in recruitment deviates (sigma_u)
###############################################################################

delistedOut.all <- list(); length(delistedOut.all) <- 3
delistedOut.all[[1]] <- readRDS("workspaces/sigma_u_delisted_baseGreen.rds")
delistedOut.all[[2]] <- readRDS("workspaces/sigma_u_delisted_baseAmber.rds")
delistedOut.all[[3]] <- readRDS("workspaces/sigma_u_delisted_baseRed.rds")

sigma_u <- seq(0.5, 1.5, 0.1)

sensitivityPar <- sigma_u
baseValue <- 1.13

X <- seq(0.5, 1.5, 0.2)
#------------------------------------------------------------------------------
# Plot
# quartz(width = 4.5, height= 7, pointsize = 9)
pdf(file = "FigS26.pdf", width = 4.5, height= 7, pointsize = 9)

par(mfcol=c(6,2), mar=c(3, 3, 1, 1), oma=c(1,3,3,4))

for(j in 1:2){ # for SR and HS
	for(s in 1:3){ # For each base case
	delistedOut <- delistedOut.all[[s]]
	
	# Ppn Wrong #-----------------------------------------------------------------
		bp <- barplot(t(delistedOut$ppn[[j]]), col=c(statusCols, grey(0.8), 1), las=1, ylab = "", space=0.1, border=NA, xaxs="i", yaxs="i", yaxt="n")
	axis(side = 2, at = seq(0, 1, 0.5), las = 1)
	axis(side = 1, at = bp[seq(1, length(bp), 2)], labels=X)
	
	points(bp[findInterval(baseValue, sensitivityPar)], 1.05, pch=8, xpd=NA)
	if(j == 1) mtext(side=2, line = 3.5, "Proportion of\nsimulations", cex=0.8)
	
	if(j==1) mtext(side=3, line=0.5, adj=0, c("a)", "e)", "i)")[s], cex=0.8)
	if(j==2) mtext(side=3, line=0.5, adj=0, c("b)", "f)", "j)")[s], cex=0.8)
	
	if(j==1 & s ==1) mtext(side = 3, line = 2, "Spawner-recruitment", font = 2, cex = 0.8)
	if(j==2 & s ==1) mtext(side = 3, line = 2, "Percentile", font = 2, cex = 0.8)
	
	# mtext(side = 3, line = 1.5, c("High productivity\nHarvest Control Rule", "Low productivity\nModerate harvest", "Low productivity\nHigh harvest")[s], cex=0.8, font=2)
	if(j == 2 & s == 1) text(14, -0.5, srt=-90,"High productivity\nHarvest Control Rule", font = 2, xpd= NA, cex = 1.2)
	if(j == 2 & s == 2) text(14, -0.5, srt=-90,"Low productivity\nModerate harvest", font = 2, xpd= NA, cex = 1.2)
	if(j == 2 & s == 3) text(14, -0.5, srt=-90,"Low productivity\nHigh harvest", font = 2, xpd= NA, cex = 1.2)
	
	# RBs #-----------------------------------------------------------------------
	jit <- 0.015
	ptCex <- 1
	
	# Spawner-recruitment
	if(j==1) 
	if(j==2) yAll <- range(c(delistedOut$RB$S25[, '75th'], delistedOut$RB$S25[, '25th'], delistedOut$RB$S50[, '25th'], delistedOut$RB$S50[, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th']))
	
	J1 <- (j-1)*2 + 2 # Index for lower
	J2 <- (j-1)*2 + 3 # Index for upper
	
	yAll <- range(c(delistedOut$RB[[J1]][, '75th'], delistedOut$RB[[J1]][, '25th'], delistedOut$RB[[J2]][, '25th'], delistedOut$RB[[J2]][, '75th'], delistedOut$RB$avgS[, '25th'], delistedOut$RB$avgS[, '75th']))
	
	plotCI(sensitivityPar, delistedOut$RB[[J1]][, 'median'], "n", li = delistedOut$RB[[J1]][, '25th'], ui= delistedOut$RB[[J1]][, '75th'], yaxt="n", col=NA, ylab="", xlab="", xaxt="n", ylim=yAll)
	axis(side=1, at = sigma_u[seq(1, length(sigma_u),2)])
	abline(v = baseValue, lwd=10, col=grey(0.8))
	abline(h=0)
	
	plotCI(sensitivityPar - jit, delistedOut$RB[[J1]][, 'median'], "n", li = delistedOut$RB[[J1]][, '25th'], ui= delistedOut$RB[[J1]][, '75th'], gap=0, sfrac = 0, col=statusCols['r'], cex=ptCex, pch=25, pt.bg=statusCols['r'], add=TRUE)
	plotCI(sensitivityPar + jit, delistedOut$RB[[J2]][, 'median'], li = delistedOut$RB[[J2]][, '25th'], ui=delistedOut$RB[[J2]][, '75th'], gap=0, sfrac = 0, add=TRUE, col=statusCols['g'], cex=ptCex, pch=24, pt.bg=statusCols['g'])
	plotCI(sensitivityPar, delistedOut$RB$avgS[, 'median'], li = delistedOut$RB$avgS[, '25th'], ui= delistedOut$RB$avgS[, '75th'], gap=0, sfrac = 0, add=TRUE, cex=ptCex, pch=21, pt.bg=1)
	
	A <- axis(side=2, labels = FALSE, at=pretty(yAll, n=4))
	axis(side=2, at = A, labels = paste(formatC(round(A*100, 1), 0, format="f"), "%", sep=""), las=1)
	if(j == 1) mtext(side=2, line=4, "Relative bias", cex=0.8)
	
	if(j==1) mtext(side=3, line=0.5, adj=0, c("c)", "g)", "k)")[s], cex=0.8) 
	if(j==2) mtext(side=3, line=0.5, adj=0, c("d)", "h)", "l)")[s], cex=0.8) 
	}
}

mtext(outer=TRUE, side=1, expression(paste("Standard deviation in recruitment deviates (", sigma[upsilon], ")")), cex=0.8)

dev.off()

#------------------------------------------------------------------------------
# How is true status (regardless of how it's assessed) affected by sigma_u?
trueStatus_sigma_u <- readRDS("workspaces/sigma_u_trueStatus.rds")

pdf(file = "FigS27.pdf", width = 3, height = 5.5, pointsize = 9)

par(mfrow = c(3,1), mar = c(4,4,2,1), oma = c(2,2,1,0))

for(baseCaseNum in 1:3){
bp <- barplot(t(trueStatus_sigma_u[[baseCaseNum]][, j, ]), col=c(statusCols), las=1, ylab = "", space=0.1, border=NA, xaxs="i", yaxs="i", yaxt="n")
axis(side = 2, at = seq(0, 1, 0.5), las = 1)
axis(side = 1, at = bp[seq(1, length(bp), 2)], labels=sigma_u[seq(1, length(sigma_u), 2)])

points(bp[findInterval(baseValue, sensitivityPar)], 1.05, pch=8, xpd=NA)

mtext(side=3, line=1, adj=0, c("a) High productivity/HCR", "b) Low productivity/Moderate harvest", "c) Low productivity/High harvest")[baseCaseNum], cex=0.8)
}

mtext(side=2, outer=TRUE, "Proportion of simulations with green, amber, and red true status", cex=0.8)
mtext(outer=TRUE, side=1, expression(paste("Standard deviation in recruitment deviates (", sigma[upsilon], ")")), cex=0.8)

dev.off()
