#______________________________________________________________________________
#' Assess performance of observed versus true status 
#' 
#' This function compares the observed benchmarks and status to the true  
#' benchmarksand status, and returns the difference between true and observed  
#' benchmarks, as well as a code indicating whether the resulting status
#' assignment was right or not (and how so).
#' 
#' @param trueStatus A list outcome from the assessPop function applied to
#' the true data.
#' @param obsStatus A list outcome from the assessPop function applied to
#' the observed data.
#' @return Returns a list, with the first element equal to a 2x2 matrix of the 
#' mean percent error ((observed - true)/true) for the upper and lower benchmarks 
#' (columns) and both SR and percentile metrics (rows) and second element equal
#' to a status code (see details).
#' @details The status code indicates the difference between true and observed
#' status outcomes for the two metrics. For each metric, the code is a numeric
#' from 1 - 9 that indicates:
#' 1: correct green status
#' 2: correct amber status
#' 3: correct red status
#' 4: assessed as red, when really amber (cautious)
#' 5: assessed as red, when really green (very cautious)
#' 6: assessed as amber, when really red (risky)
#' 7: assessed as amber, when really green (cautious)
#' 8: assessed as green, when really red (very risky)
#' 9: assessed as green, when really amber (risky)
#' 
#' @examples 
#' 

perfStatus <- function(trueStatus, obsStatus) {
	
	benchMPE <- cbind(
		(obsStatus$lowerBenchmark - trueStatus$lowerBenchmark)/trueStatus$lowerBenchmark, 
		(obsStatus$upperBenchmark - trueStatus$upperBenchmark)/trueStatus$upperBenchmark)
	colnames(benchMPE) <- c("lower", "upper")
	rownames(benchMPE) <- c("SR", "perc")
	
	statusDiff <- numeric(2)
	for (i in 1:2){ # For each metric
		if (trueStatus$statusNum[i] == obsStatus$statusNum[i]){
			statusDiff[i] <- obsStatus$statusNum[i]
		} else if (trueStatus$statusNum[i] == 2 & obsStatus$statusNum[i] == 3){
			statusDiff[i] <- 4
		} else if (trueStatus$statusNum[i] == 1 & obsStatus$statusNum[i] == 3){
			statusDiff[i] <- 5
		} else if (trueStatus$statusNum[i] == 3 & obsStatus$statusNum[i] == 2){
			statusDiff[i] <- 6
		} else if (trueStatus$statusNum[i] == 1 & obsStatus$statusNum[i] == 2){
			statusDiff[i] <- 7
		} else if (trueStatus$statusNum[i] == 3 & obsStatus$statusNum[i] == 1){
			statusDiff[i] <- 8
		} else if (trueStatus$statusNum[i] == 2 & obsStatus$statusNum[i] == 1){
			statusDiff[i] <- 9
		} else {
			stop("Status assessments not defined!")
		}
	}
	
	return(list(MPE = benchMPE, statusDiff = statusDiff))
	
}
