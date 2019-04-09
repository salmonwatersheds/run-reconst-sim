#______________________________________________________________________________
#' Make list of parameter sets to be run in a sensitivity analysis 
#' 
#' This function takes a numeric vector of values for a given parameter in
#' `simPar` and returns a list of parameter sets for input into
#' `runSensitivity()`.
#' 
#' @param basePar A named vector of the base parameter values to be applied.
#' @param sensName The character name of the parameter that is being varied,
#' matching exactly as the name appears in `basePar`
#' @param sensValues A vector of values for the parameter named in sensName
#' over which the sensitivity analysis is to be run.
#' @return Returns a list, where each element is a unique parameter combination
#' to be run in a sensitivity analysis. 
#' 
#' @details 
#' 
#' @examples 
#' 

makeParList <- function(basePar = simPar, sensName = "obs_bias", sensValues = log(1/seq(1, 5, 0.2))){
	
	# Check that parameter to vary is correctly named in basePar
	if(is.element(sensName, names(basePar)) == FALSE){
		stop("sensName not in basePar")
	}
	
	# Create list of parameters for each value of sensValues
	parList <- list(); length(parList) <- length(sensValues)
	for(i in 1:length(sensValues)){
		parList[[i]] <- basePar
		parList[[i]][which(names(basePar) == sensName)] <- sensValues[i]
	}
	
	return(parList)
}

#______________________________________________________________________________
#' Run a sensitivitity analysis over a list of parameter sets 
#' 
#' This function takes a list of parameters and runs multiple MC simulations
#' of `reconstrSim()` for each parameter set (i.e., element) in that list.
#' 
#' @param parList A list of parameter sets, where each element is a named
#' vector of all parameters being fed into `reconstrSim()`. Can be the output
#' of `makeParList()` above.
#' @param nSim The number of MC simulations to be run for each parameter set
#' @param nCores THe number of cores over which to run the MC simulations in
#' parallel.
#' @return Returns a list, with the first element equal to the overall output 
#' of the sensitivity analysis (see Details) and the second element equal to
#' the total time the simulation took in minutes.
#' 
#' @details The overall output is given as a list, where each element is the
#' performace output for the corresponding parameter set in the given `parList`.
#' The performace is itself a list, where:
#' 1. the first element is the relative bias (RB, a.k.a. mean percent error 
#' MPE) in the current spawner abundance, 
#' 2. the second element is the RB in the stock-recruit benchmarks and the 
#' (mis)classification of status (as a numebr from 1-8; see documentation for
#' `perfStatus()`, each given as columns in a matrix with rows equal to `nSim`,
#' 3. the third element is the RB and (mis)classification according to the
#' historical spawners metric.
#' 
#' @examples 
#' 

runSensitivity <- function(parList, nSim, nCores){
	
	# Parallelize over multiple cores
	registerDoParallel(cores = nCores)
	
	ptime <- system.time({ 
		
	out_overall <- foreach (j = 1:length(parList)) %dopar% {
		
		performance <- list(
			avgS = rep(NA, nSim),
			SR = matrix(NA, nrow = nSim, ncol = 3),
			HS = matrix(NA, nrow = nSim, ncol = 3))
		
		for(i in 1:nSim){
			out <- reconstrSim(parList[[j]])
			
			performance[[1]][i] <- out$performance$currentMPE
			for(b in 2:3){ # SR benchmarks (b = 2), or HS benchmarks (b = 3)
				performance[[b]][i, 1:2] <- out$performance$benchMPE[b-1, ]
				performance[[b]][i, 3] <- out$performance$statusDiff[b-1]
			}
			
		} # end nSim
		
		performance
	} # end over all parList
	
})[3]
	
	return(list(out_overall, time = ptime/60))
}


#______________________________________________________________________________
#' Delist output of `runSensitivity` 
#' 
#' This function ...
#' 
#' @param 
#' @param 
#' @return 
#' 
#' @details 
#' 
#' @examples 
#' 

delistSensitivity <- function(out){
	
	nSim <- length(out[[1]]$avgS)
	n <- length(out)
	
	RB <- list(
		avgS =  matrix(NA, nrow = n, ncol = 4, dimnames=list(NULL, c("mean", "median", "25th", "75th"))),
		Sgen1 = matrix(NA, nrow = n, ncol = 4, dimnames=list(NULL, c("mean", "median", "25th", "75th"))), 
		Smsy = matrix(NA, nrow = n, ncol = 4, dimnames=list(NULL, c("mean", "median", "25th", "75th"))), 
		S25 = matrix(NA, nrow = n, ncol = 4, dimnames=list(NULL, c("mean", "median", "25th", "75th"))),
		S50 = matrix(NA, nrow = n, ncol = 4, dimnames=list(NULL, c("mean", "median", "25th", "75th")))
	)
	
	for (i in 1:n) {
		
		RB$avgS[i, 'mean'] <- mean(out[[i]]$avgS)
		RB$avgS[i, 'median'] <- quantile(out[[i]]$avgS, 0.5, na.rm=TRUE)
		RB$avgS[i, '25th'] <- quantile(out[[i]]$avgS, 0.25, na.rm=TRUE)
		RB$avgS[i, '75th'] <- quantile(out[[i]]$avgS, 0.75, na.rm=TRUE)
		
		RB$S25[i, 'mean'] <- mean(out[[i]]$HS[,1], na.rm=TRUE)
		RB$S25[i, 'median'] <- quantile(out[[i]]$HS[,1], 0.5, na.rm=TRUE)
		RB$S25[i, '25th'] <- quantile(out[[i]]$HS[,1], 0.25, na.rm=TRUE)
		RB$S25[i, '75th'] <- quantile(out[[i]]$HS[,1], 0.75, na.rm=TRUE)
		
		RB$S50[i, 'mean'] <- mean(out[[i]]$HS[,2], na.rm=TRUE)
		RB$S50[i, 'median'] <- quantile(out[[i]]$HS[,2], 0.5, na.rm=TRUE)
		RB$S50[i, '25th'] <- quantile(out[[i]]$HS[,2], 0.25, na.rm=TRUE)
		RB$S50[i, '75th'] <- quantile(out[[i]]$HS[,2], 0.75, na.rm=TRUE)
		
		RB$Sgen1[i, 'mean'] <- mean(out[[i]]$SR[,1], na.rm=TRUE)
		RB$Sgen1[i, 'median'] <- quantile(out[[i]]$SR[,1], 0.5, na.rm=TRUE)
		RB$Sgen1[i, '25th'] <- quantile(out[[i]]$SR[,1], 0.25, na.rm=TRUE)
		RB$Sgen1[i, '75th'] <- quantile(out[[i]]$SR[,1], 0.75, na.rm=TRUE)
		
		RB$Smsy[i, 'mean'] <- mean(out[[i]]$SR[,2], na.rm=TRUE)
		RB$Smsy[i, 'median'] <- quantile(out[[i]]$SR[,2], 0.5, na.rm=TRUE)
		RB$Smsy[i, '25th'] <- quantile(out[[i]]$SR[,2], 0.25, na.rm=TRUE)
		RB$Smsy[i, '75th'] <- quantile(out[[i]]$SR[,2], 0.75, na.rm=TRUE)
	}
	
	#----------------------------------------------------
	# Proportion of simulations with wrong status
	
	ppn <- list(
		SR = matrix(NA, nrow = n, ncol = 5, dimnames = list(NULL, c("green", 'amber', 'red', 'pessimistic', 'optimistic'))), 
		HS = matrix(NA, nrow = n, ncol = 5, dimnames = list(NULL, c("green", 'amber', 'red', 'pessimistic', 'optimistic')))
	)
	for (i in 1:n) {
		for(j in 1:2) {
			ppn[[j]][i, 1] <- sum(out[[i]][[j+1]][, 3] == 1)/nSim
			ppn[[j]][i, 2] <- sum(out[[i]][[j+1]][, 3] == 2)/nSim
			ppn[[j]][i, 3] <- sum(out[[i]][[j+1]][, 3] == 3)/nSim
			ppn[[j]][i, 4] <- sum(is.element(out[[i]][[j+1]][, 3], c(4,5,7)))/nSim
			ppn[[j]][i, 5] <- sum(is.element(out[[i]][[j+1]][, 3], c(6,8,9)))/nSim
			}}
	
	return(list(RB = RB, ppn = ppn))
}
	