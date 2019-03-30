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


