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


reconstrSim <- function(simPar, cuCustomCorrMat=NULL, seed = NULL, returnObsCases = FALSE) {
	
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
	# Draw productivity parameters
	a <- rnorm(nPop, simPar$a_mean, simPar$sigma_a)
	
	#_____
	# Determine density dependence parameter and change over time
	# Draw Smax - spawner abundance that yields maximum recruitment (above which there's neg dens dep)
	Smax.initial <- c(
		rlnorm(simPar$nIndicator, simPar$logSmax_ind_mean, simPar$logSmax_ind_sd), # Indicator streams
		rlnorm(simPar$nNonIndicator, simPar$logSmax_nonInd_mean, simPar$logSmax_nonInd_sd) #Non-indicator
		)
	
	# Check
	if(length(Smax.initial) != nPop) stop("Number of indicator and non-indicator streams does not match nPop.")
	
	# Create matrix of annual capacity, assuming no changes over time
	Smax <- matrix(rep(Smax.initial, each = nYears), nrow = nYears, ncol = nPop, byrow = FALSE)
	
	# If a decline in Smax is to be included, alter the capacity of affected subpopulations
	if(simPar$greenHab < 100){
		if(sum(c(simPar$greenHab, simPar$amberHab, simPar$redHab)) != 100){
			stop("Habitat status categories sum to > 100%")
		}
		dum <- rmultinom(n = nPop, size = 1, prob = c(simPar$greenHab, simPar$amberHab, simPar$redHab) / 100)
		habitatCategory <- t(dum) %*% matrix(c(1:3), nrow=3)
		SmaxDecl <- matrix(NA, nrow = nYears, ncol = nPop)
		for(j in 1:nPop){
			if(habitatCategory[j] == 1){
				SmaxDecl[, j] <- 1
			} else if (habitatCategory[j] == 2){
				SmaxDecl[, j] <- c(rep(1, (simPar$gen + 2)), seq(1, runif(1, 0.50, 0.75), length.out = simYears))
			} else if (habitatCategory[j] == 3){
				SmaxDecl[, j] <- c(rep(1, (simPar$gen + 2)), seq(1, runif(1, 0.25, 0.50), length.out = simYears))
			}
		} # end nPop
		# Calculate Smax incorporating decline
		Smax <- Smax * SmaxDecl
	} # end if
	
	# Calculate density dependence parameter as the inverse of spawner abundance for max. recruitment
	b <- 1/Smax
	
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
	# If there is a fixed target harvest rate, can calculate
	# realized harvest rates (assumed constant across subpopulations)
	# outside of population dynamics loop
	if(simPar$harvContRule == "noError"){
		harvestRate <- rep(simPar$targetHarvest, nYears)
	} else if(simPar$harvContRule == "fixedER"){
		harvestRate <- realizedHarvestRate(
			targetHarvest = simPar$targetHarvest, 
			sigmaHarvest = simPar$sigma_harvest,
			nYears = nYears, 
			errorType = "beta")
		targetHarvest <- c(rep(0, nYears - simPar$simYears), rep(simPar$targetHarvest, simPar$simYears))
	} else {
		targetHarvest <- c(rep(0, nYears - simPar$simYears), rep(NA, simPar$simYears))
		harvestRate <- numeric(nYears)
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
		stop("Covariance matrix for residuals is not positive definite.")
		# covMat <- make.positive.definite(covMat)
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
	# Changed to 20% Smax to avoid negative numbers if a < 0
	# and recruitment error for first year
	spawners[1:(simPar$gen + 2), ] <- 0.2 * 1 / b[1:(simPar$gen + 2), ]
	phi[1,] <- 0
	
	# Loop over first 7 years for chum (par$gen + 2)
	for (y in 1:(simPar$gen + 2)){ #first obsLag period necessary to generate recBY
		dum <- rickerModel(S = spawners[y, ],
											 a = a, # constant productivity
											 b = b[y, ], # time-varying capacity
											 error = rmvnorm(1, rep(0, nPop), sigma = covMat), #-simPar$sigma_u^2 / 2 # No lognormal bias correction
											 rho = simPar$rho,
											 phi_last = phi[y, ],
											 recCap = simPar$recCap, 
											 extinctThresh = 0)#simPar$extinctThresh)
		recruitsBY[y, ] <- dum[[1]]
		phi[y+1, ] <- dum[[2]]
	}
	
	#_____
	# LOOP 2: Simulate population dynamics to nYears
	#_____
	
	for(y in (simPar$gen + 3):nYears){
		
		# Sum return from recruitment t-2 to t-6 years ago * ppn age at return
		recruitsRY[y, ] <- ppnAge[cbind(y - ages, 1:simPar$gen)] %*% recruitsBY[y - ages,]
		
		# If using variable harvest rate, calculate targetHarvest based on true total return
		if(simPar$harvContRule == "variableER"){
			targetHarvest[y] <- round(simPar$maxHarvest * (1 - exp(simPar$d * (simPar$m - sum(recruitsRY[y, ])))), 4)
			harvestRate[y] <- realizedHarvestRate(
				targetHarvest = targetHarvest[y], 
				sigmaHarvest = simPar$sigma_harvest,
				nYears = 1, 
				errorType = "beta")
		}
		
		spawners[y, ] <- (1 - harvestRate[y]) * recruitsRY[y, ]
		trueCatch[y] <- sum(harvestRate[y] * recruitsRY[y, ])
		
		# Apply Ricker model to calculate recruits
		dum <- rickerModel(S = spawners[y, ], 
											 a = a, 
											 b = b[y, ], 
											 error = rmvnorm(1, rep(0, nPop), sigma = covMat), #-simPar$sigma_u^2 / 2 # No lognormal bias correction
											 rho = simPar$rho,
											 phi_last = phi[y, ],
											 recCap = simPar$recCap, 
											 extinctThresh = 0)#simPar$extinctThresh)
		
		recruitsBY[y, ] <- dum[[1]]
		phi[y+1, ] <- dum[[2]]
		
	} # end nYears y
	
	#-----------------------------------------------------------------------------
	# Observation submodel
	#-----------------------------------------------------------------------------
	# Observations only made over all years including initialization
	
	#_____
	# Add lognormal observation error to all spawners
	if(is.na(simPar$obs_biasChangeYear) == FALSE){
		obs_bias <- c(rep(simPar$obs_bias, simPar$obs_biasChangeYear), rep(simPar$obs_bias2, simYears - simPar$obs_biasChangeYear))
	} else {
		obs_bias <- rep(simPar$obs_bias, simYears)
	}
	
	obsSpawners <- spawners[(simPar$gen + 3):nYears, ] * matrix(exp(
	qnorm(runif(simYears*nPop, 0.0001, 0.9999), 
				rep(obs_bias, nPop), #- simPar$sigma_obs^2 / 2, # No lognormal bias correction
				simPar$sigma_obs)), 
	nrow = simYears, ncol = nPop, byrow = FALSE)
	
	#_____
	# Apply monitoring design (random sampling of propSampled each year)
	# Different for indicator (first) and non-indicator (second)
	z <- cbind(
			samplingDesign(ppnSampled = simPar$ppnSampled_ind, 
										 nPop = simPar$nIndicator, 
										 simYears = simYears,
										 ppnChange = simPar$ppnChange_ind,
										 samplingDeclStart = simPar$samplingDeclStart_ind,
										 samplingDeclEnd = simPar$samplingDeclEnd_ind,
										 gen = simPar$gen),
			samplingDesign(ppnSampled = simPar$ppnSampled_nonInd, 
										 nPop = simPar$nNonIndicator, 
										 simYears = simYears,
										 ppnChange = simPar$ppnChange_nonInd,
										 samplingDeclStart = simPar$samplingDeclStart_nonInd,
										 samplingDeclEnd = simPar$samplingDeclEnd_nonInd,
										 gen = simPar$gen))
	
	# Constraint: at least one indicator stream has to be monitored each year,
	# otherwise you get a value for ExpFactor1 of Inf
	nIndicatorMonitored <- apply(z[,1:simPar$nIndicator] == 1, 1, sum)
	if(length(which(nIndicatorMonitored == 0)) > 0){
		# Randomly choose one indicator stream to be monitored for each year that has none monitored
		z[cbind(which(nIndicatorMonitored == 0), sample(1:simPar$nIndicator, size = length(which(nIndicatorMonitored == 0)), replace = TRUE))] <- 1
	}
	
	
	# #---
	# # Check: proportion of all indicator and non-indicator streams with at least one year of data
	# # in each decade
	# decade.counts <- matrix(NA, nrow = 5, ncol = nPop)
	# for(i in 1:nPop) decade.counts[ ,i] <- tapply(z[,i], rep(1:5, each=10), sum)
	# rbind(apply(decade.counts[, 1:simPar$nIndicator] > 0, 1, sum)/simPar$nIndicator, apply(decade.counts[, (simPar$nIndicator+1):nPop] > 0, 1, sum)/simPar$nNonIndicator)
	# #---
	
	sampledSpawners <- z * obsSpawners # 0 = not monitored
	
	#_____
	# Add error to observed catch for CU
	obsCatch <- trueCatch[(simPar$gen + 3):nYears] * exp(qnorm(runif(simYears, 0.0001, 0.9999), 
					simPar$catch_bias,# - simPar$sigma_obs^2 / 2, # No lognormal bias correction
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
	spawnersExp3 <- simPar$ExpFactor3 * spawnersExp2
	
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
	obsStatus <- assessPop(SR.pairs = obsData, gen = 4)
	
	#_____
	# Benchmarks: true
	trueData <- data.frame(
		S = apply(spawners[(simPar$gen + 3):nYears, ], 1, sum), 
		R = apply(recruitsBY[(simPar$gen + 3):nYears, ], 1, sum))	
	
	# trueStatus.data <- assessPop(SR.pairs = trueData, gen = simPar$gen)
	# Q: In assessing true status under a decline in capacity, set a nominal period 
	# or assess based on current SR parameters?
	# A: Use nominal period that reflects initial capacity parameters to avoid shifting
	# baselines.
	trueStatus <- assessTruePop(SR.pairs = trueData, SR.params = cbind(a, b[1, ]), gen = 4)
	
	#-----------------------------------------------------------------------------
	# Performance
	#-----------------------------------------------------------------------------
	# Base case: use true status derived from parameters for SR metric
	# trueStatus <- trueStatus.params 
	
	P <- perfStatus(trueStatus, obsStatus)
	
	# ****************************************************************************
	if(returnObsCases == TRUE){
		# Two categories of partial application of the observation submodel:
		# a) Perfect observation but incomplete coverage (separates out effect
		#		 of Expansion Factors I and II)
		# Include incomplete monitoring of spawners observed without a bias
		# to separate out effects of Expansion Factor I and II without III
		sampledSpawners_noBias <- z * spawners[(simPar$gen + 3):nYears, ]
		
		dumExp1a <- ExpFactor1(sampledSpawners = sampledSpawners_noBias[, 1:simPar$nIndicator])
			spawnersExp1a <- dumExp1a[[1]] * apply(sampledSpawners_noBias[, 1:simPar$nIndicator], 1, sum)
			
			dumExp2a <- ExpFactor2(
				spawnersInd = sampledSpawners_noBias[, 1:simPar$nIndicator], 
				spawnersNonInd = sampledSpawners_noBias[, (simPar$nIndicator + 1):simPar$nPop])
			spawnersExp2a <- dumExp2a[[1]] * spawnersExp1a
			
			spawnersExp3a <- 1 * spawnersExp2a # Expansion Factor III = 1 because no obs_bias
			
			obsReturna <- obsCatch + spawnersExp3a
			
			recruitsShifteda <- matrix(unlist(shift(x = obsReturna, n = ages[obsPpnAge!=0], type = "lead")), ncol = length(ages[obsPpnAge!=0]))
			
			obsRecruitsBYa <-  recruitsShifteda %*% obsPpnAge[obsPpnAge!=0]
			
			obsDataa <- data.frame(S = spawnersExp3a, R = obsRecruitsBYa)
			obsStatusa <- assessPop(SR.pairs = obsDataa, gen = simPar$gen)
			
	# b) Imperfect observation but complete coverage (separates out effect
	#		 of Expansion Factor III)
			spawnersExp3b <- simPar$ExpFactor3 * apply(obsSpawners, 1, sum)
			
			obsReturnb <- obsCatch + spawnersExp3b
			
			recruitsShiftedb <- matrix(unlist(shift(x = obsReturnb, n = ages[obsPpnAge!=0], type = "lead")), ncol = length(ages[obsPpnAge!=0]))
			
			obsRecruitsBYb <-  recruitsShiftedb %*% obsPpnAge[obsPpnAge!=0]
			
			obsDatab <- data.frame(S = spawnersExp3b, R = obsRecruitsBYb)
			obsStatusb <- assessPop(SR.pairs = obsDatab, gen = simPar$gen)
			
			# Performance 
			Pa <- perfStatus(trueStatus, obsStatusa)
			Pb <- perfStatus(trueStatus, obsStatusb)
			
	} # end ObsCases
	# ****************************************************************************
	
	#-----------------------------------------------------------------------------
	# END
	#-----------------------------------------------------------------------------
	
	if(returnObsCases == TRUE){
		return(list(
			performance = list(
				trueParams = P, 
				caseA = Pa, 
				caseB = Pb),
			status = list(
				true = trueStatus, 
				obs = obsStatus, 
				obsA = obsStatusa, 
				obsB = obsStatusb),
			data = list(true = trueData, obs = obsData),
			RickerPar = list(a = a, b = b, Smax = Smax),
			trueHarvest = data.frame(
				targetHarvest = targetHarvest[(simPar$gen + 3):nYears], 
				realizedHarvest = harvestRate[(simPar$gen + 3):nYears], 
				trueCatch = trueCatch[(simPar$gen + 3):nYears],
				obsCatch = obsCatch[(simPar$gen + 3):nYears])
		))
	} else {
		return(list(
			performance = P,
			status = list(
				true = trueStatus, 
				obs = obsStatus),
			data = list(true = trueData, obs = obsData),
			RickerPar = list(a = a, b = b, Smax = Smax),
			trueHarvest = data.frame(
				targetHarvest = targetHarvest[(simPar$gen + 3):nYears], 
				realizedHarvest = harvestRate[(simPar$gen + 3):nYears], 
				trueCatch = trueCatch[(simPar$gen + 3):nYears],
				obsCatch = obsCatch[(simPar$gen + 3):nYears])
		))
	}
	
} # end recoverySim function
