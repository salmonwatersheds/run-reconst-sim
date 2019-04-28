#' Plot the matrix of true vs observed status assessments
#'
#'
#' @param statusDiff A numeric vector with length nSim giving the status
#' difference code returned by perfStatus 
#' @return Plot
#'
#' @examples
#' #
#' 

plotStatusDiff <- function(statusDiff){
	
	statusCols <- c(g = "#8EB687", a = "#DFD98D", r = "#B66A64")
	# Create matrix for number of simulations with true status R, A, G (columns) 
	# and observed status R, A, G (rows)
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
	par(mar=c(4, 5, 2, 1))
	plot(1:4, 1:4, "n", xaxt = "n", yaxt = "n", xlab = "True status", ylab = "")
	for(s in 1:2){
		axis(side = s, at = c(1:4), labels=FALSE)
	}
	axis(side = 1, at = seq(1.5, 3.5, 1), c("Red", "Amber", "Green"), tck = 0, las=1)
	axis(side = 2, at = rev(seq(1.5, 3.5, 1)), c("Red", "Amber", "Green"), tck = 0, las=1)
	
	mtext(side = 2, line = 4, "Observed status")
	for(i in 1:3){
		for(j in 1:3){
			if(i == j) polygon(x = c(j, j+1, j+1, j), y = 4 - c(i , i , i-1, i-1), border = NA, col = statusCols[4-i])
			if(j < i){ # Risky misclassifications
				fnt <- 2 
				fntcol <- 1
			}else if (j > i){ # Cautious misclassifications
				fnt <- 1
				fntcol <- grey(0.6)
			} else {
				fnt <- 1
				fntcol <- 1
			}
			text(j + 0.5, 4 - i + 0.5, paste(formatC(round(statusDiffMat[i,j]/length(statusDiff)*100, 1), 1, format="f"), "%"), font=fnt, col=fntcol)
			# text(j + 0.5, 4 - i + 0.5, statusDiffMat[i,j])
		}
	}
	
} # end function
