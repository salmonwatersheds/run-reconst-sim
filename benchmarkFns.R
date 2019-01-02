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
	# Add warning and adjustment for non-zero spawner abundances
	# if(any(s < 0.00001)){
	# 	s[s < 0.00001] <- 0.0001
	# 	# print(c("s abundance must be > 0. Negative values replaced w/ small positive"))
	# }
	
	# if(any(s.msy < 0.00001)){
	# 	s.msy[s.msy < 0.00001] <- 0.0001
	# 	# print(c("s.msy must be > 0. Negative values replaced w/ small positive"))
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
#' This function calculates Sgen, or the spawner abundances that would result 
#' in recovery to Smsy within one generation. 
#' 
#' @param 
#' @param 
#' @return
#' 
#' @examples 
#' 

calcSgen <- function(Sgen.hat, theta, Smsy) {
	
	fit <- optimize(f = Sgen.optim, interval = c(0, ((theta[1] / theta[2]) * (0.5 - 0.07 * theta[1]))))
}
	

solver.sgen <- function(s, theta, s.msy) 
{
	#fit=optimize(f=fn.sgen,interval=c(0,S.msy[m]), theta=theta, s.msy=s.msy)	
	fit = optimize(f = fn.sgen,interval = c(0, ((theta[1] / theta[2]) * (0.5 - 0.07 * theta[1]))), 
								 theta = theta, s.msy = s.msy)         
	#if(fit$convergence!=0)print("Non-convergence for benchmark S.lower.3")
	return(list(fit = fit$minimum))
}