runSensitivity <- function(
	basePar = simPar, 
	sensitivity = list("obs_bias", log(1/seq(1, 5, 0.2))), 
	nSim = 6000,
	nCores = 8,
	multiPar = FALSE
	){
	
	if(multiPar == FALSE & length(basePar[[1]]) > 1){
		stop("Error: mutliple base parameter sets supplied but multiPar == FALSE")
	}
	
	if(multiPar == TRUE & is.null(sensitivity) == FALSE){
		stop("Looking for multiple base parameter sets (multiPar = TRUE) but single-value sensitivity supplied.")
	}
	
	if(multiPar == FALSE & is.null(sensitivity)){
		stop("Looking for single-value sensitivity supplied (multiPar = FALSE) but sensitivity is NULL.")
	}
	
	if(multiPar == FALSE){
		if(is.element(sensitivity[[1]], names(basePar)) == FALSE){
			stop("Error: Sensitivity parameter not named in basePar.")
		}}
	
	# Parallelize over multiple cores
	registerDoParallel(cores = nCores)
	
	# If supplied a list of all  parameters instead of 
	# range of values for a single parameter
	if(multiPar == TRUE){
		n.i <- length(basePar)
	} else {
		n.i <- length(sensitivity[[2]])
	}
	
	ptime <- system.time({ 
		
		out_overall <- foreach (i = 1:n.i) %dopar% {
			
			if(multiPar == TRUE){
				simPar.i <- basePar[[i]]
			} else {
				simPar.i <- basePar
				simPar.i[which(names(basePar) == sensitivity[[1]])] <- sensitivity[[2]][i]
			}
			
			dummy <- list(
				avgS = rep(NA, nSim),
				SR = matrix(NA, nrow = nSim, ncol = 3),
				HS = matrix(NA, nrow = nSim, ncol = 3))
			
			performance <- list(all = dummy, a = dummy, b=dummy, true.data = dummy)
			
			for(i in 1:nSim){
				out <- reconstrSim(simPar.i)
				
				for(j in 1:4){ #1 = trueParams, 2 = perfect obs (no obs_bias), 3 = complete coverage, 4 = trueData
					# avgS (b = 1)
					performance[[j]][[1]][i] <- out$performance[[j]]$currentMPE
					for(b in 2:3){ # SR benchmarks (b = 2), or HS benchmarks (b = 3)
						performance[[j]][[b]][i, 1:2] <- out$performance[[j]]$benchMPE[b-1, ]
						performance[[j]][[b]][i, 3] <- out$performance[[j]]$statusDiff[b-1]
					}}
			}
			
			performance
		}
		
	})[3]
	
	return(list(out_overall, time = ptime/60))
}


