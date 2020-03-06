#' Generate recruit abundance with Ricker model
#'
#' This function calculates recruitment from Ricker curve with AR(1) process
#' (according to Peterman et al. 2003; modified to take more recent parameter-
#' ization). Uses parameters from arima.mle (a, -b, sig, tau in log space) with
#' multivariate normally distributed errors. Note that by default
#' utminus1 and tau are zero, resulting in a standard Ricker model.
#'
#' @param S A numeric vector of spawner abundances.
#' @param a A numeric vector of alpha values, i.e. productivity at low spawner
#' abundance.
#' @param b A numeric vector of beta values, i.e. density dependence para-
#' meter.
#' @param error A numeric vector of recruitment errors (upsilon), typically generated
#' using \code{rmvnorm()} and relevant process variance estimates (sigma).
#' @param tau A numeric vector of temporal autocorrelation values, i.e. AR1 coefficient.
#' outside of model using multivariate normal (or equivalent) distribution.
#' @param phi_last A numeric vector representing recruitment deviations (phi) from
#' previous brood year (t-1).
#' @param recCap The recruitment cap. If R > recCap, then set R = recCap.
#' @param extinctThresh The extinction threshold. If S < extinctThresh, set R = 0
#' @return Returns a list of R, a numeric representing recruit abundance, and
#' \code{ut}, the recruitment deviation for year t, which is used to generate 
#' subsequent process error.
#'
#' @examples
#' #Spawner and recruit values represent millions of fish, stock-recruit
#' parameters approximate those of Fraser River sockeye salmon Chilko CU.
#'
#' #without autoregressive error
#' rickerModel(S = 1.1, a = 1.8, b = 1.2, error = 0.3)
#'
#' #with autoregressive error
#' rickerModel(S = 1.1, a = 1.8, b = 1.2, error = 0.3, tau = 0.2,
#' phi_last = 0.7)
#' 
#' # For 10 subpopulations
#' nPop <- 10
#' Sigma <- matrix(rep(0.01, nPop^2), nPop, nPop)
#' diag(Sigma) <- 0.1
#' rickerModel(S = runif(nPop, 0.8, 1.5), a = rnorm(nPop, 1, 1), 
#' b = rep(1, nPop), error = rmvnorm(1, rep(0, nPop), Sigma), tau = 0.2,
#' phi_last = rmvnorm(1, rep(0, nPop), Sigma))
#'
#' @export

rickerModel <- function(S, a, b, error, tau = 0, phi_last = 0, recCap = NULL, extinctThresh = 0) {
	
	
	# err <- utminus1 * tau + error
	# phi <- tau * phi_last + error # WRONG! As pointed out by Reviewer #2, this needs to be corrected to account for effect of autocorrelation on the variance
	phi <- tau * phi_last + error
	
	# if (a >= 0) {
	# 	if (b != 0 & S > 0) {
			R <- S * exp(a - b * S) * exp(phi)
	# 		err.next <- log(R / S) - (a - b * S) + error
	# 	}
	# 	if (b == 0 & S > 0) {
	# 		R <- S * exp(err)
	# 		err.next <- log(R / S) - 0
	# 	}
	# }
	# if (a < 0 & S > 0) {
	# 	R <- S * exp(a) * exp(error)
	# 	err.next <- log(R / S) - 0 + error
	# }
	# if (S == 0) {
	# 	R <- 0
	# 	err.next <- err
	# }
	# return(list(R, err.next))
	
	# If any subpopulations are below extinction threshold,
	# then set R = 0 for those subpopulations		
	if(sum(S <= extinctThresh) > 0){ 
		R[which(S <= extinctThresh)] <- 0
	}

	if(length(recCap) > 0){ #If recCap is not NULL
		R[which(R > recCap)] <- recCap
	}

	return(list(R, phi))
}

#______________________________________________________________________________

#' Generate log-normal error associated with proportions data
#'
#' This function generates proporations at age with multivariate logistic error
#' (Schnute and Richards 1995, eqns.S.9 and S.10). 
#'
#' @param ppnAgeVec A numeric vector of mean proportions of individuals returning 
#' at a given age. \code{nAges = length(ppnAgeVec)}.
#' @param omega A numeric specifying the parameter that controls interannual 
#' variability in proportions of fish returning at each age.
#' @param nYears A numeric vector of length 1 giving the number of years to 
#' generate random ages for.
#' @return Returns a numeric matrix, \code{p}, representing proportions for 
#' each class, with number of rows equal to \code{nYears} and number of columns 
#' equal to \code{nAges} 
#'
#' @examples
#' ppnAgeErr(ppnAgeVec = c(0.2, 0.4, 0.3, 0.1), omega = 0.8, nYear = 1)
#'
#' @export

ppnAgeErr <- function(ppnAgeVec, omega, nYears) {
	nAges <- length(ppnAgeVec)
	
	#NAs produced when dividing 0 by 0 for ppn of recruits at age; replace w/ 0s
	ppnAgeVec[is.na(ppnAgeVec)] <- 0
	
	# Calculate matrix of random normal deviates, epsilon
	epsilon <- matrix(qnorm(runif(nYears*nAges, 0.0001, 0.9999)), nrow = nYears, ncol = nAges)
	
	# Dummy variable in order to ensure proportions sum to one.
	p.dum <- matrix(rep(ppnAgeVec, each=nYears), nrow = nYears, ncol = nAges) * exp(omega*epsilon)
	
	p <- p.dum/apply(p.dum, 1, sum)
	return(p)
}


#______________________________________________________________________________
#' Calclate realized harvest rate given target (mean) and standard deviation
#'
#' This function generates harvest rates bound between zero and one,
#' incorporating either beta or normal error around a target harvest rate, h'.
#'
#' @param targetHarvest A single value or vector (for temporal changes in 
#' harvest rates) giving the target harvest rate
#' @param sigmaHarvest A single value for the standard deviation in error
#' around the target harvest rate
#' @param n The number of realized harvest rates to return. If 
#' targetHarvest is a vector, then n must equal length(targetHarvest). May be 
#' the number of years or number of populations.
#' @param errorType Must be one of "beta" or "normal" specifying the error
#' distribution for the realized harvest rates.  Default is beta.
#' 
#' @return Returns a numeric vector of realized harvest rates for each n
#'
#' @examples
#'
#' @export

realizedHarvestRate <- function(targetHarvest, sigmaHarvest, n = length(targetHarvest), errorType = "beta") {
	
	# Checks
	if(errorType != "beta" & errorType != "normal"){
		stop("Unknown error distribution. Must be beta or normal.")
	}
	
	if(length(targetHarvest) > 1 & length(targetHarvest) != n){
		stop("If length(targetHarvest)>1, then must equal n")
	}
	
	#-----------------------------
	# BETA error
	if(errorType == "beta"){
		
		beta1 <- (targetHarvest^2 - targetHarvest^3 - sigmaHarvest^2*targetHarvest)/(sigmaHarvest^2)
		beta2 <- (targetHarvest * (1 - targetHarvest)^2 - sigmaHarvest^2*(1 - targetHarvest))/(sigmaHarvest^2)
		
		harvestRate <- rbeta(n = n, shape1 = beta1, shape2 = beta2)
		
		# BETA CHECK #1:
		# Beta distribution requires shape parameters > 0. If this is violated, then
		# the above harvestRate will have NaNs. Check and fix:
		if(length(which(beta1 <= 0 & beta2 <= 0)) > 0){
			
			# For beta1, beta2 > 0, we must satisfy
			# targetHarvest * (1 - targetHarvest) > sigmaHarvest^2
			# Calculate threshold targetHarvest for given sigmaHarvest:
			tH <- seq(0, 1, 0.0001)
			unacceptable <- tH[which(tH*(1-tH) <= sigmaHarvest^2)]
			unacceptableInd <- which(round(diff(unacceptable), 3) > 0.0001)
			thresh <- unacceptable[c(unacceptableInd, unacceptableInd+1)]
			# The acceptable range is just inside the unacceptable
			thresh <- c(thresh[1] + 0.0001, thresh[2] - 0.0001)
			
			# For each violation of the above rule
			for(j in which(beta1 <= 0 & beta2 <=0)){ 
				
				# 1) If targetHarvest[j] is closer to zero than one:
				if(abs(diff(c(targetHarvest[j], 0))) < abs(diff(c(targetHarvest[j], 1)))){
					
					# a) If target harvest is closer to zero than threshold, 
					if(abs(diff(c(targetHarvest[j], 0))) <= abs(diff(c(targetHarvest[j], thresh[1])))){
						harvestRate[j] <- 0 # Choose realized harvest rate = 0
					} else { # b) If target harvest rate is closer to threshold, draw from that
						harvestRate[j] <- rbeta(n = 1, shape1 = (thresh[1]^2 - thresh[1]^3 - sigmaHarvest^2*thresh[1])/(sigmaHarvest^2), shape2 = (thresh[1] * (1 - thresh[1])^2 - sigmaHarvest^2*(1 - thresh[1]))/(sigmaHarvest^2))
					}
					
			# 2) If targetHarvest[j] is closer to one than zero:		
				} else {
					# a) If target harvest is closer to one than threshold, 
					if(abs(diff(c(targetHarvest[j], 1))) <= abs(diff(c(targetHarvest[j], thresh[2])))){
						harvestRate[j] <- 1 # Choose realized harvest rate = 0
					} else { # b) If target harvest rate is closer to threshold, draw from that
						harvestRate[j] <- rbeta(n = 1, shape1 = (thresh[2]^2 - thresh[2]^3 - sigmaHarvest^2*thresh[2])/(sigmaHarvest^2), shape2 = (thresh[2] * (1 - thresh[2])^2 - sigmaHarvest^2*(1 - thresh[2]))/(sigmaHarvest^2))
					}
			}
			}
		} # END BETA CHECK #1
		
		# BETA CHECK #2
		# The beta distribution is bimodal under certain conditions. If the realized 
		# harvest rate is much different from the target, then re-draw
		bimodalInd <- which(targetHarvest^2 *(1 - targetHarvest)/(1 + targetHarvest) < sigmaHarvest ^ 2 & targetHarvest * (1-targetHarvest)^2 / (2-targetHarvest)	< sigmaHarvest^2)
		
		# If there are instances where targetHarvest indicates we're in the bimodal zone
		if(length(bimodalInd) > 0){
			
			# If we are in the bimodal range and the harvestRate is > 0.5 different
			# from the targetHarvest, then assume it's a error to to bimodality (wrongEnd)
			wrongEnd <- which(abs(harvestRate[bimodalInd]-targetHarvest[bimodalInd]) > 0.5)
			
			while(length(wrongEnd) > 0){
				# Redraw offending bimodal values
				harvestRate[bimodalInd[wrongEnd]] <- rbeta(n = length(wrongEnd), shape1 = beta1[bimodalInd[wrongEnd]], shape2 = beta2[bimodalInd[wrongEnd]])

				# Re-calculate if there are any that come from opposite ends
				wrongEnd <- which(abs(harvestRate[bimodalInd]-targetHarvest[bimodalInd]) > 0.5)
			} 
		} # END BETA CHECK #2
			
	# Finally, check that there are no NaNs from something else I missed!
		if(sum(is.na(harvestRate)) > 0){
			warning
			(paste("Realized harvest still contains NAs for targetHarvest = ", targetHarvest[is.na(harvestRate)], "\n Imputed with harvest = 0"))
			harvestRate[is.na(harvestRate)] <- 0
		}
		 
	#-----------------------------
	# NORMAL error
	} else if (errorType == "normal"){
		
		# Equation (F19) from Holt et al. (2018 CSAS)
		harvestRate <- targetHarvest + qnorm(runif(n, 0.0001, 0.9999), 0, sigmaHarvest)
		
		# If harvest rate is < 0 or > 1, resample as in Holt et al. (2018)
		while (length(which(harvestRate > 1 | harvestRate < 0)) > 0) {
			if(length(targetHarvest) == 1){
				harvestRate[which(harvestRate > 1 | harvestRate < 0)] <- targetHarvest + qnorm(runif(length(which(harvestRate > 1 | harvestRate < 0)), 0.0001, 0.9999), 0, sigmaHarvest)
			} else if(length(targetHarvest) > 1){
				harvestRate[which(harvestRate > 1 | harvestRate < 0)] <- targetHarvest[which(harvestRate > 1 | harvestRate < 0)] + qnorm(runif(length(which(harvestRate > 1 | harvestRate < 0)), 0.0001, 0.9999), 0, sigmaHarvest)
			} # end if
		} #end while
	} # end normal error

	return(harvestRate)
	}

