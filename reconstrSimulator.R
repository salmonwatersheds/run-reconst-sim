#' Recovery simulator
#'
#' Closed-loop simulation tool to assess management procedures and inform
#' Pacific salmon rebuilding strategies. Based on C. Holt's south coast chum 
#' model (CSAS 2018) and the samSim package (C. Freshwater) . Function runs
#' a single MCMC trial, and can therefore be run in parallel using dopar 
#' if many MCMC runs are desired. 
#' 
#' ** Dec 20, 2018: May make more sense to run MCMC trials within function, and parallelize
#' across parameters/scenarios. Will consider this later. **
#' 
#' The model generates salmon stock-recruit data including variation in age
#' structure, survey design, and variable exploitation rules. OM uses Ricker
#' formulation \code{R=S(exp(a-bS))}. 
#' 
#' @importFrom here here
#' @importFrom dplyr group_by summarise
#' @param simPar Vector of parameters for the simulation from .csv file
#' @param cuCustomCorrMat Vector of parameters for the simulation from .csv file

#' @return TO BE DEFINED.
#' @export

#Temporary inputs
# here <- here::here
simPar <- read.csv(here("data/baseSimPar.csv"), stringsAsFactors = F)
# cuCustomCorrMat <- read.csv(here("data/baseCorrMatrix.csv"), stringsAsFactors=F)

set.seed(987)
a <- rnorm(simPar$nPop, simPar$a_mean, simPar$sigma_a)

recoverySim <- function(simPar, ricker_a = a, cuCustomCorrMat=NULL,
												dirName, uniqueProd=TRUE, seed = NULL) {
	
	# If a seed for simulation is provided, then set the seed for the simulation here
	if(length(seed) == 1){ 
		set.seed(seed)
	}
	
	#-----------------------------------------------------------------------------
	# Setup
	#-----------------------------------------------------------------------------
	
	# Unpack parameters referred to often
	nPop <- simPar$nPop
	nYears <- (simPar$gen + 2) + simPar$simYears
	ages <- 2:6
	
	#_____
	# Calculate proportion at age  
	ppnAge <- ppnAgeErr(
		ppnAgeVec = c(simPar$meanRec2, simPar$meanRec3, simPar$meanRec4, 
									simPar$meanRec5, simPar$meanRec6), 
		omega = simPar$ageErr, 
		nYears = nYears) 
	
	if(length(ages) != dim(ppnAge)[2]) { 
		stop("Length of 'ages' doesn't match 'ppnAgeVec'.") 
		}
	
	#_____
	# Calculate true harvest rates (assumed constant across subpopulations)
	# Equation (F19) from Holt et al. (2018 CSAS)
	# ** Target harvest rates are not changing through time for now, but
	# this is a second priority to vary. **
	harvestRate <- simPar$targetHarvest + qnorm(runif(nYears, 0.0001, 0.9999), 0, simPar$sigma_harvest)
	# If harvest rate is < 0 or > 1, resample as in Holt et al. (2018)
	while (length(which(harvestRate > 1 | harvestRate < 0)) > 0) {
		harvestRate[which(harvestRate > 1 | harvestRate < 0)] <- simPar$targetHarvest + qnorm(runif(length(which(harvestRate > 1 | harvestRate < 0)), 0.0001, 0.9999), 0, simPar$sigma_harvest)
	}
	
	#_____
	# Calculate corMat: correlation matrix for recruitment residuals
	# referred to as "Sigma" in equations in Appendix F of Holt et al. (2018 CSAS)
	# Set covariance among subpopulations in recruitment deviations
	if (simPar$corrMat == FALSE) { 
		correlPop <- simPar$correlPop 
	} else if (simPar$corrMat == TRUE) { # Replace uniform correlation w/ custom matrix
		if (nrow(cuCustomCorrMat) != nPop) {
			stop("Custom correlation matrix does not match number of subpopulations.")
		} else {
			correlPop <- as.matrix(cuCustomCorrMat)
		}
	}
	
	sigMat <- matrix(as.numeric(simPar$sigma_u), nrow = 1, ncol = nPop) 
	covMat <- t(sigMat) %*% sigMat # Calculate shared variance
	corMat <- covMat * correlPop # Correct based on correlation
	diag(corMat) <- as.numeric(simPar$sigma_u^2) # Add variance
	
	#-----------------------------------------------------------------------------
	# Population dynamics submodel
	#-----------------------------------------------------------------------------
	
	# Set up matrices to store output
	phi <- matrix(NA, nrow = nYears + 1, ncol = nPop) # recruitment deviations
	recruitsBY <- matrix(NA, nrow = nYears, ncol = nPop) # recruits by brood year
	recruitsRY <- matrix(NA, nrow = nYears, ncol = nPop) # recruits by return year
	spawners <- matrix(NA, nrow = nYears, ncol = nPop) # spawners in each year
	
	#_____
	# LOOP 1: Initialize subpopulations over a generation with constant age structure  
	# and no harvest
	#_____
	
	# Initialize spawners at 20% of Seq (Holt 2018 CSAS) 
	# and recruitment error for first year
	spawners[1:(simPar$gen + 2), ] <- 0.2 * a / simPar$b
	phi[1,]<-rmvnorm(1, rep(0, nPop), sigma = corMat)
	
	# Loop over first 7 years for chum (par$gen + 2)
	for (y in 1:(simPar$gen + 2)){ #first obsLag period necessary to generate recBY
		dum <- rickerModel(S = spawners[y, ], a = a, b = rep(simPar$b, nPop), 
											 error = rmvnorm(1, rep(0, nPop), sigma = corMat),
											 rho = simPar$rho,
											 phi_last = phi[y, ])
		recruitsBY[y, ] <- apply(rbind(dum[[1]], simPar$recCap), 2, min)
		phi[y+1, ] <- dum[[2]]
	}
	
	#_____
	# LOOP 2: Simulate population dynamics to nYears
	#_____
	
	for(y in (simPar$gen + 3):nYears){
		
		# Sum return from recruitment t-2 to t-6 years ago * ppn age at return
		recruitsRY[y, ] <- ppnAge[cbind(y - ages, 1:simPar$gen)] %*% recruitsBY[y - ages,]
		
		spawners[y, ] <- harvestRate[y] * recruitsRY[y, ]
		
		# Apply Ricker model to calculate recruits
		dum <- rickerModel(S = spawners[y, ], a = a, b = rep(simPar$b, nPop), 
											 error = rmvnorm(1, rep(0, nPop), sigma = corMat),
											 rho = simPar$rho,
											 phi_last = phi[y, ])
		
		recruitsBY[y, ] <- apply(rbind(dum[[1]], simPar$recCap), 2, min)
		phi[y+1, ] <- dum[[2]]
		
	} # end nYears y
	
	#-----------------------------------------------------------------------------
	# Observation submodel
	#-----------------------------------------------------------------------------
	
	#-----------------------------------------------------------------------------
	# Assessment submodel
	#-----------------------------------------------------------------------------
	
	return()
	
} # end recoverySim function