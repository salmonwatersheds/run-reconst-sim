#' Recovery simulator
#'
#' Closed-loop simulation tool to assess management procedures and inform
#' Pacific salmon rebuilding strategies. Based on C. Holt's south coast chum 
#' model (CSAS 2018) and the samSim package (C. Freshwater) . Function runs
#' a single MCMC trial, and can therefore be run in parallel using dopar 
#' if many MCMC runs are desired. 
#' 
#' ** Dec 20, 2018: May make more sense to run MCMC trials within function, and 
#' parallelize across parameters/scenarios. Will consider this later. **
#' 
#' The model generates salmon stock-recruit data including variation in age
#' structure, survey design, and variable exploitation rules. OM uses Ricker
#' formulation \code{R=S(exp(a-bS))}. 
#' 
#' @importFrom here here
#' @importFrom dplyr group_by summarise
#' @param simPar Vector of parameters for the simulation from .csv file, including:
#'   * Simulation basics *
#'   scenario: names of given scenario corresponding to these parameters
#'   nameOM: ??
#'   nameMP: ??
#'   species: salmon species being simulated (always chum for runReconst project)
#'   simYears: number of years to simulate and use for status assessment (not 
#'   including initialization for age structure and recruitment calculation)
#'   nIndicator: number of indicator streams within the hypothetical CU
#'   nNonIndicator: number of non-indicator streams within the hypothetical CU
#'   nPop: total number of streams or sub-populations (nIndicator + nNonIndicator)
#'   * Population submodel *
#'   a_mean: average productivity parameter (a) among sub-populations within CU
#'   sigma_a: standard deviation in a among sub-populations within CU
#'   b: density dependence parameter (same for all subpopulations)
#'   rho: temporal autocorrelation in Ricker residuals (numeric 0-1)
#'   sigma_u: standard deviation in residuals (same for all subpopulations;
#'   sigma_upsilon)
#'   correlPop: spatial autocorrelation among subpopulations (covariance; 0-1)
#'   corrMat: is a custom correlation matrix for recruitment deviations (Sigma)
#'   provided? If FALSE, then corrMat calculated from sigma_u and correlPop
#'   recCap: recruitment cap
#'   gen: number of different ages that salmon can return at
#'   meanRec2: average proportion of fish that return as age 2
#'   meanRec3: average proportion of fish that return as age 3
#'   meanRec4: average proportion of fish that return as age 4
#'   meanRec5: average proportion of fish that return as age 5
#'   meanRec6: average proportion of fish that return as age 6
#'   ageErr: magnitude of natural interannual variation in age-at-return by 
#'   brood year, referred to as \bar{\omega} by Holt et al. (2018 CSAS)
#'   * Harvest submodel *
#'   harvContRule: harvest control rule (always fixedER for runReconst project)
#'   targetHarvest: the target harvest rate for the CU
#'   sigma_harvest: the standard deviation in differences between target and
#'   realized harvest rates
#'   * Observation submodel *
#'   sigma_obs: standard deviation in observation error of spawners
#'   ppnSampled: proportion of indicator streams and proportion of non-indicator 
#'   streams that are monitored in a given year
#'   extinctThresh: number of spawners below which the subpopulation is 
#'   considered extinct
#'   
#' @param ricker_a Numeric vector of productivity parameters for each 
#' subpopulation within the CU, so can be held constant across MCMC trials.
#' Drawn from \code{rnorm(simPar$nPop, simPar$a_mean, simPar$sigma_a)}.  
#'
#' @param cuCustomCorrMat If corrMat == FALSE, this is applied as a custom 
#' matrix of spatial autocorrelation among subpopulations. Allows flexibility
#' for correlPop to differ among different pairs of subpopulations.
#'
#' @param seed Random number seed for the given MCMC trial. Can set to ensure
#' results of simulation are reproducible, but should be different for each
#' MCMC trial!!
#' 
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

	# true aggregate (i.e., summed across subpopulations) catch each year based on 
	# realized harvest rate and true number of returns to each subpopulation
	trueCatch <- numeric(nYears) 
	
	#_____
	# LOOP 1: Initialize subpopulations over a generation with constant age structure  
	# and no harvest
	#_____
	
	# Initialize spawners at 20% of Seq (Holt 2018 CSAS) 
	# and recruitment error for first year
	spawners[1:(simPar$gen + 2), ] <- 0.2 * a / simPar$b
	phi[1,]<-0
	
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
		
		spawners[y, ] <- (1 - harvestRate[y]) * recruitsRY[y, ]
		trueCatch[y] <- sum(harvestRate[y] * recruitsRY[y, ])
		
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
	
	#_____
	# Add lognormal observation error to all spawners
	obsSpawners <- spawners * matrix(exp(
		qnorm(runif(nYears*nPop, 0.0001, 0.9999), 
					0, # -simPar$sigma_obs^2 / 2, # Change this to zero???
					simPar$sigma_obs)), 
		nrow = nYears, ncol = nPop)
	
	#_____
	# Add error to observed catch for CU
	# ** If adding bias in catch, this is where you'd put it **
	obsCatch <- trueCatch * exp(qnorm(runif(nYears, 0.0001, 0.9999), 
					0, #-simPar$sigma_obs^2 / 2, # Change this to zero???
					simPar$sigma_catch))
	
	#_____
	# Calculate observed average age-at-return
	# ** This may change to subsampling from true ages (and adding obs error?) **
	obsPpnAge <- ppnAgeErr(
		ppnAgeVec = c(simPar$meanRec2, simPar$meanRec3, simPar$meanRec4, 
									simPar$meanRec5, simPar$meanRec6), 
		omega = simPar$obsAgeErr, 
		nYears = 1) 
	
	#_____
	# Apply monitoring design (random sampling of propSampled each year)
	
	sampled <- samplingDesign(ppnSampled = simPar$ppnSampled, nPop, nYears,
														ppnChange = simPar$ppnChange, 
														samplingDeclStart = simPar$samplingDeclStart,
														samplingDeclEnd = simPar$samplingDeclEnd,
														gen = simPar$gen)
	
	#_____
	
	
	#-----------------------------------------------------------------------------
	# Assessment submodel
	#-----------------------------------------------------------------------------
	
	return()
	
} # end recoverySim function