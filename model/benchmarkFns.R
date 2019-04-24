#' Calculate upper stock-recruitment benchmark, Smsy
#'
#' This function calculates Smsy, or the spawner abundance projected to 
#' maintain long-term maximum sustainable yield from a population with 
#' Ricker dynamics. It applied the explicit solution for Smsy given by
#' Scheuerell (2016), PeerJ, DOI 10.7717/peerj.1623
#' This function uses the lambert_W0 function from the `gsl` library.
#'
#' @param a A numeric value giving Ricker parameter a (or log alpha), estimated
#' from the observed spawner-recruitment data.
#' @param b A numeric value giving Ricker parameter b (strength of density 
#' dependence; units spawners^(-1)), estimated from the observed spawner-
#' recruitment data.
#' @return Returns the value of Smsy.
#'
#' @examples
#' #

calcSmsy <- function(a, b) {
	
	Smsy = (1 - lambert_W0(exp(1 - a))) / b
	
	return(as.numeric(Smsy))
}

#______________________________________________________________________________
#' Optimization routine for calculation of Sgen
#' 
#' This function calculates the likelihood of residuals between the projected
#' recruits from an estimated Sgen.hat and a known value of Smsy, to be used in
#' the optimization of Sgen in the calcSgen function. This function is based on
#' the Sgen.optim function in the samSim package (Freshwater et al. 2018).
#' 
#' @param Sgen.hat A proposed numeric value for Sgen
#' @param theta A numeric vector containing the estimated parameters a, b, and
#' sigma from the Ricker function fit to data
#' @param Smsy The calculated value of Smsy based on theta, from the calcSmsy
#' function.
#' @return Returns the negative log likelihood of the residuals.

Sgen.optim <- function (Sgen.hat, theta, Smsy) {
	# # Add warning and adjustment for non-zero spawner abundances
	# if(any(Sgen.hat < 0.00001)){
	# 	Sgen.hat[Sgen.hat < 0.00001] <- 0.0001
	# 	print(c("Sgen.hat abundance must be > 0. Negative values replaced w/ small positive"))
	# }
	# 
	# if(any(Smsy < 0.00001)){
	# 	Smsy[Smsy < 0.00001] <- 0.0001
	# 	print(c("Smsy must be > 0. Negative values replaced w/ small positive"))
	# }
	
	a <- theta[1]
	b <- theta[2]
	sig <- exp(theta[3])
	
	# Compute projected recruits based on Sgen.hat
	Smsy.hat <- Sgen.hat * exp(a - b * Sgen.hat) 
	
	# Calculate residuals and negative log likelihood
	epsilon <- log(Smsy) - log(Smsy.hat)
	nloglike <- - sum(dnorm(x = epsilon, mean = 0, sd = sig, log = TRUE))	
	
	return(nloglike)
}
	
#______________________________________________________________________________
#' Calculate lower stock-recruitment benchmark, Sgen: 
#' 
#' This function calculates Sgen1, or the spawner abundances that would result 
#' in recovery to Smsy within one generation. 
#' 
#' @param Sgen.hat A proposed numeric value for Sgen
#' @param theta A numeric vector containing the estimated parameters a, b, and
#' sigma from the Ricker function fit to data
#' @param Smsy The calculated value of Smsy based on theta, from the calcSmsy
#' function.
#' @return Returns the numeric value of Sgen1.
#' 
#' @examples 
#' 

calcSgen <- function(Sgen.hat, theta, Smsy) {
	
	fit <- optimize(f = Sgen.optim, interval = c(0, Smsy), theta = theta, Smsy = Smsy)
	
	# Give warning if Sgen1 is at upper or lower bound
	if(round(fit$minimum) == 0) warning("Sgen1 at lower bound of zero")
	
	if(round(fit$minimum) == round(Smsy)){
		warning("Lower benchmark greater than upper benchmark (Sgen1 > Smsy). Set to NA.")
		return(Sgen1 = NA)
	} else {
		return(Sgen1 = as.numeric(fit$minimum))
	}
}
	
#______________________________________________________________________________
#' Return status category given upper and lower benchmarks for a given metric
#' 
#' This function returns the status category as "green", "amber", or "red"
#' given upper and lower benchmarks and the current value for a given 
#' metric. If the upper and lower benchmarks are vectors of length 3, giving 
#' the mean, upper, and lower 95% credible intervals for the benchmarks, then
#' the function will return the probability of falling within each status
#' category.
#' 
#' @param current The current value of the metric (e.g., spawner abundance
#' over the most recent generation)
#' @param lower The lower benchmark, either a single numeric value or a vector
#' of length 3, giving the mean, upper, and lower 95% credible intervals for 
#' the lower benchmark
#' @param upper The upper benchmark, either a single numeric value or a vector
#' of length 3, giving the mean, upper, and lower 95% credible intervals for 
#' the upper benchmark
#' @return Returns a list containing the status (green, amber, or red), status 
#' number number (1 = green, 2 = amber, 3 = red), current, lower, and upper 
#' benchmarks, and (if confidence intervals on benchmarks are given) the 
#' probability of being in each status category.
#' 
#' @examples 
#' 

assessMetric <- function(current, lower, upper) {
	
	if(is.na(current) == TRUE) warning("Current abundance NA")
	if(is.na(lower) == TRUE) warning("Lower benchmark NA")
	if(is.na(upper) == TRUE) warning("Upper benchmark NA")
	
	if(is.na(lower) == TRUE & is.na(upper) == TRUE){
		status <- list(NA, 0, current, lower, upper)
	} else {
		if(length(lower) != length(upper)) stop("Lower and upper benchmarks must be of the same length")
		
		# If just a single benchmark is given (no CI):
		if(length(lower) == 1){
			# 1) current >= upper -> green
			if (is.na(upper) == FALSE & current >= upper){
				status <- list("green", 1, current, lower, upper)
			# 2) current <= lower -> red
			} else if (is.na(lower) == FALSE & current <= lower){
				status <- list("red", 3, current, lower, upper)
			# 3) no lower -> current <= upper = lower -> red
			} else if (is.na(lower) == TRUE & current <= upper){
				status <- list("red", 3, current, lower, upper)
			# 4) 
			} else {
				status <- list("amber", 2, current, lower, upper)
			}
		} # end if length(lower) == 1
		
		if(length(lower) == 3){
			stop("Probabilities for status not yet implemented for this function.")
		}}
		
	return(status)
}

#______________________________________________________________________________
#' Return population assessment based on spawner abundance using both
#' stock-recruitment and historic spawners metrics 
#' 
#' This function returns the status category as "green", "amber", or "red",
#' current geometric-mean spawner abundance, and upper and lower benchmarks 
#' for both stock-recruitment (SR) and historic spawners (HS) benchmarks.
#' 
#' @param SR.pairs A data frame containing S (total spawner abundance to the CU)
#' and R (recruits to the CU by brood year corresponding to the number of 
#' spawners in the S)
#' @param gen The length of a generation (for calculating geometric-mean 
#' spawners over most recent generation)
#' @return Returns a list containing the status for SR and HS metrics
#' (green, amber, or red), corresponding status number (1 = green, 2 = amber, 
#' 3 = red), current gemetric-mean spawner abundance, lower and upper 
#' benchmarks for both metrics.
#' 
#' @examples 
#' 

assessPop <- function(SR.pairs, gen) {
	
	S <- SR.pairs$S
	R <- SR.pairs$R
	
	# Calculate geometric mean spawners over most recent generation:
	AvgS <- prod(tail(S, gen))^(1/gen)
	
	# Stock-recruit benchmarks
	
	# Reconstructed run (rr) to estimate Ricker a and b parameters:
	rr <- data.frame(x = S, y = log(R/S)) 
	# plot(rr$x, rr$y)
	
	fit <- lm(y ~ x, data = rr)
	theta <- c(a = as.numeric(fit$coefficients[1]),
						 b = - as.numeric(fit$coefficients[2]),
						 sig = as.numeric(summary(fit)$sigma))
	
	# Upper SR benchmark (Smsy)
	Smsy <- calcSmsy(a = theta[1], b = theta[2])
	
	# Lower SR benchmark (Sgen1)
	# Start optimization at 20% of Smsy
	Sgen1 <- calcSgen(Sgen.hat = 0.2*Smsy, theta = theta, Smsy = Smsy)
	
	# Changed upper benchmark to 80% Smsy
	statusSR <- assessMetric(current = AvgS, lower = Sgen1, upper = 0.8*Smsy)
	
	# historic spawners benchmarks
	
	# Upper BM changed to 50th percentile
	upperP <- quantile(S, 0.5)
	
	# Based on Table 6 of Holt et al. (2018), if alpha < 2.5, then
	# upper and lower benchmarks collapse to S50
	if(exp(theta['a']) < 2.5){
		lowerP <- NA 
	} else {
		lowerP <- quantile(S, 0.25)
	}
	statusHS <- assessMetric(current = AvgS, lower = lowerP, upper = upperP)
	
	return(list(
		status = c(SR = statusSR[[1]], HS = statusHS[[1]]),
		statusNum = c(SR = statusSR[[2]], HS = statusHS[[2]]),
		current = AvgS,
		lowerBenchmark = c(SR = Sgen1, HS = lowerP),
		upperBenchmark = c(SR = Smsy, HS = upperP)
		))
	
}

#______________________________________________________________________________
#' Return population assessment based on spawner abundance given
#' true underlying parameters for stock-recruitment relationship 
#' 
#' This function returns the status category as "green", "amber", or "red",
#' current geometric-mean spawner abundance, and upper and lower benchmarks 
#' for both stock-recruitment (SR) and historic spawners (HS) benchmarks.
#' 
#' @param SR.pairs A data frame containing S (total spawner abundance to the CU)
#' and R (recruits to the CU by brood year corresponding to the number of 
#' spawners in the S)
#' @param SR.params A matrix with two columns, one for the productivity 
#' parameter (a) and one for the density-dependence parameter (b), where rows
#' are the different sub-populations within the CU
#' @param gen The length of a generation (for calculating geometric-mean 
#' spawners over most recent generation)
#' @return Returns a list containing the status for SR and HS metrics
#' (green, amber, or red), corresponding status number (1 = green, 2 = amber, 
#' 3 = red), current gemetric-mean spawner abundance, lower and upper 
#' benchmarks for both metrics.
#' 
#' @examples 
#' 

assessTruePop <- function(SR.pairs, SR.params, gen) {
	
	S <- SR.pairs$S
	R <- SR.pairs$R
	nPop <- dim(SR.params)[1]
	
	# Calculate geometric mean spawners over most recent generation:
	AvgS <- prod(tail(S, gen))^(1/gen)
	
	# Stock-recruit benchmarks
	
	# Upper SR benchmark (Smsy)
	Smsy <- pmax(0, calcSmsy(a = SR.params[,1], b = SR.params[,2]))
	
	# Check
	if(sum(is.na(Smsy)) > 0){
		stop("Error: NAs in Smsy for true benchmarks calculated from actual parameters.")
	}
	
	# Lower SR benchmark (Sgen1)
	# Start optimization at 20% of Smsy
	Sgen1 <- rep(NA, nPop)
	for (i in which(Smsy > 0)) { # For each subpopulation, optimize Sgen1
		Sgen1[i] <- calcSgen(
			Sgen.hat = 0.2*Smsy[i], 
			theta = c(
				a = as.numeric(SR.params[i,1]), 
				b = as.numeric(SR.params[i,2]), 
				sig = 10^-10), 
			Smsy = Smsy[i])
	}
	
	# # Check
	# if(sum(is.na(Sgen1)) > 0){
	# 	stop("Error: NAs in Sgen1 for true benchmarks calculated from actual parameters.")
	# }
	
	# In assessing true status, remove NAs (equivalent to setting to zero) for lower or
	# upper benchmarks? Most common that low Smsy -> Sgen = NA and lower benchmark is NA
	statusSR <- assessMetric(current = AvgS, lower = sum(Sgen1, na.rm=TRUE), upper = 0.8*sum(Smsy))
	
	# According to Holt et al. (2018) there is no "true" historic spawners status
	# The only "true" benchmarks are based on underlying SR parameters
	# # historic spawners benchmarks
	# 
	# lowerP <- quantile(S, 0.25)
	# upperP <- quantile(S, 0.75)
	# statusHS <- assessMetric(current = AvgS, lower = lowerP, upper = upperP)
	
	# Return status SR for comparison to observed SR and observed HS
	return(list(
		status = c(SR = statusSR[[1]], HS = statusSR[[1]]),
		statusNum = c(SR = statusSR[[2]], HS = statusSR[[2]]),
		current = AvgS,
		lowerBenchmark = c(SR = sum(Sgen1), HS = sum(Sgen1)),
		upperBenchmark = c(SR = sum(Smsy), HS = sum(Smsy))
	))
	
}
