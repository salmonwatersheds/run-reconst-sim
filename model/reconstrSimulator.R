#' Run-reconstruction simulator
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


reconstrSim <- function(simPar, a = a, cuCustomCorrMat=NULL,
												dirName = NULL, uniqueProd=TRUE, seed = NULL) {
	
	# If a seed for simulation is provided, then set the seed for the simulation here
	if(length(seed) == 1){ 
		set.seed(seed)
	}
	
	#-----------------------------------------------------------------------------
	# Setup
	#-----------------------------------------------------------------------------
	
	# Unpack parameters referred to often
	nPop <- simPar$nPop
	simYears <- simPar$simYears
	nYears <- (simPar$gen + 2) + simYears
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
	# Calculate covMat: covariance matrix for recruitment residuals
	# See https://en.wikipedia.org/wiki/Covariance_and_correlation
	# cov(X,Y) = cor(X,Y) * sd(X) * sd(Y)
	# cov(X,Y) referred to as "Sigma" in equations in Appendix F of Holt et al. (2018 CSAS)
	# but referred to as corMat in samSim package
	
	# a) Set *correlation* among subpopulations in recruitment deviations (correlPop)
	if (simPar$corrMat == FALSE) { 
		correlPop <- simPar$correlPop 
	} else if (simPar$corrMat == TRUE) { # Replace uniform correlation w/ custom matrix
		if (nrow(cuCustomCorrMat) != nPop) {
			stop("Custom correlation matrix does not match number of subpopulations.")
		} else {
			correlPop <- as.matrix(cuCustomCorrMat)
		}
	}
	
	# b) Create matrix of variances sd(X) * sd(Y)
	# Note: we assume sd(X) = sd(Y)
	sigMat <- matrix(as.numeric(simPar$sigma_u), nrow = 1, ncol = nPop) 
	varMat <- t(sigMat) %*% sigMat # Calculate shared variance
	
	# c) Multiply by correlation to get covariance
	covMat <- correlPop * varMat # Correct based on correlation
	
	# d) Fix diagonals to be true variance within populations
	diag(covMat) <- as.numeric(simPar$sigma_u^2) # Add variance
	
	# Check "sigma" for rmvnorm() is positive semi-definite:
	if(is.positive.definite(covMat) == FALSE){
		warning("Covariance matrix for residuals is not positive definite. Adjusted to conform using make.positive.definite()")
		covMat <- make.positive.definite(covMat)
	}
	
	
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
											 error = rmvnorm(1, rep(0, nPop), sigma = covMat),
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
											 error = rmvnorm(1, rep(0, nPop), sigma = covMat),
											 rho = simPar$rho,
											 phi_last = phi[y, ])
		
		recruitsBY[y, ] <- apply(rbind(dum[[1]], simPar$recCap), 2, min)
		phi[y+1, ] <- dum[[2]]
		
	} # end nYears y
	
	#-----------------------------------------------------------------------------
	# Observation submodel
	#-----------------------------------------------------------------------------
	# Observations only made over all years including initialization
	
	#_____
	# Add lognormal observation error to all spawners
	obsSpawners <- spawners[(simPar$gen + 3):nYears, ] * matrix(exp(
		qnorm(runif(simYears*nPop, 0.0001, 0.9999), 
					0, # -simPar$sigma_obs^2 / 2, # Change this to zero???
					simPar$sigma_obs)), 
		nrow = simYears, ncol = nPop)
	
	#_____
	# Apply monitoring design (random sampling of propSampled each year)
	z <- samplingDesign(ppnSampled = simPar$ppnSampled, nPop, simYears,
														ppnChange = simPar$ppnChange, 
														samplingDeclStart = simPar$samplingDeclStart,
														samplingDeclEnd = simPar$samplingDeclEnd,
														gen = simPar$gen)
	
	sampledSpawners <- z * obsSpawners # 0 = not monitored
	
	#_____
	# Add error to observed catch for CU
	# ** If adding bias in catch, this is where you'd put it **
	obsCatch <- trueCatch[(simPar$gen + 3):nYears] * exp(qnorm(runif(simYears, 0.0001, 0.9999), 
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
	
	#-----------------------------------------------------------------------------
	# Assessment submodel
	#-----------------------------------------------------------------------------
	
	#_____
	# Expansion Factors to calculate total spawners
	
	# Expansion Factor I to account for indicator streams not monitored
	dumExp1 <- ExpFactor1(sampledSpawners = sampledSpawners[, 1:simPar$nIndicator])
	spawnersExp1 <- dumExp1[[1]] * apply(sampledSpawners[, 1:simPar$nIndicator], 1, sum)
	
	# Expansion Factor II to account for non-indicator streams
	dumExp2 <- ExpFactor2(
		spawnersInd = sampledSpawners[, 1:simPar$nIndicator], 
		spawnersNonInd = sampledSpawners[, (simPar$nIndicator + 1):simPar$nPop])
	spawnersExp2 <- dumExp2[[1]] * spawnersExp1
	
	# Expansion Factor II to account for observer efficiency
	# spawnersExp3 <- simPar$ExpFactor3 * spawnersExp2
	spawnersExp3 <- 1 * spawnersExp2
	
	#_____
	# Reconstructing recruitment
	# Observed returns based on total spawners and obsCatch
	obsReturn <- obsCatch + spawnersExp3
	
	# Observed recruits by brood year
	# Note: included the obsPpnAge!=0 so that NAs aren't produced when we're 
	# missing, e.g., age 6 returns
	recruitsShifted <- matrix(unlist(shift(x = obsReturn, n = ages[obsPpnAge!=0], type = "lead")), ncol = length(ages[obsPpnAge!=0]))
	obsRecruitsBY <-  recruitsShifted %*% obsPpnAge[obsPpnAge!=0]

	# Sanity check: Does the fancy matrix jiggery pokery work?
	# obsReturn[4] * obsPpnAge[2] + obsReturn[5] * obsPpnAge[3] + obsReturn[6] * obsPpnAge[4]
	# obsRecruitsBY[1]
	
	#_____
	# Benchmarks: observed
	obsData <- data.frame(S = spawnersExp3, R = obsRecruitsBY)
	obsStatus <- assessPop(SR.pairs = obsData, gen = simPar$gen)
	
	#_____
	# Benchmarks: true
	trueData <- data.frame(S = apply(spawners[(simPar$gen + 3):nYears, ], 1, sum), R = apply(recruitsBY[(simPar$gen + 3):nYears, ], 1, sum))	
	trueStatus <- assessPop(SR.pairs = trueData, gen = simPar$gen)
	
	#-----------------------------------------------------------------------------
	# Performance
	#-----------------------------------------------------------------------------
	
	P <- perfStatus(trueStatus, obsStatus)
	
	#-----------------------------------------------------------------------------
	# END
	#-----------------------------------------------------------------------------
	
	return(list(
		performance = P, 
		status = list(obsStatus, trueStatus),
		data = list(obsData, trueData)
	))
	
} # end recoverySim function
