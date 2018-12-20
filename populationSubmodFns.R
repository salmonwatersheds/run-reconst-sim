#' Generate recruit abundance with Ricker model
#'
#' This function calculates recruitment from Ricker curve with AR(1) process
#' (according to Peterman et al. 2003; modified to take more recent parameter-
#' ization). Uses parameters from arima.mle (a, -b, sig, rho in log space) with
#' multivariate normally distributed errors. Note that by default
#' utminus1 and rho are zero, resulting in a standard Ricker model.
#'
#' @param S A numeric vector of spawner abundances.
#' @param a A numeric vector of alpha values, i.e. productivity at low spawner
#' abundance.
#' @param b A numeric vector of beta values, i.e. density dependence para-
#' meter.
#' @param error A numeric vector of recruitment errors (upsilon), typically generated
#' using \code{rmvnorm()} and relevant process variance estimates (sigma).
#' @param rho A numeric vector of rho values, i.e. AR1 coefficient.
#' outside of model using multivariate normal (or equivalent) distribution.
#' @param phi_last A numeric vector representing recruitment deviations (phi) from
#' previous brood year (t-1).
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
#' rickerModel(S = 1.1, a = 1.8, b = 1.2, error = 0.3, rho = 0.2,
#' phi_last = 0.7)
#' 
#' # For 10 subpopulations
#' nPop <- 10
#' Sigma <- matrix(rep(0.01, nPop^2), nPop, nPop)
#' diag(Sigma) <- 0.1
#' rickerModel(S = runif(nPop, 0.8, 1.5), a = rnorm(nPop, 1, 1), 
#' b = rep(1, nPop), error = rmvnorm(1, rep(0, nPop), Sigma), rho = 0.2,
#' phi_last = rmvnorm(1, rep(0, nPop), Sigma))
#'
#' @export

rickerModel <- function(S, a, b, error, rho = 0, phi_last = 0) {
	
	# err <- utminus1 * rho + error
	phi <- rho * phi_last + error
	
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
