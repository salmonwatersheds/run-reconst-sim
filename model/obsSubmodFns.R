#' Sampling matrix for 'monitoring' spawners
#'
#' This function computes the sampling design over multiple years and streams
#' given either constant or changing monitoring coverage (proportion of streams
#' that are monitored each year). In the case of declining coverage, a start 
#' and end year for the decline can be specified. 
#' Note: in it's current form this assumed the same proportion of indicator and 
#' non-indicator streams are monitored each year, but this could be revised 
#' if desired.
#' 
#'
#' @param ppnSampled The initial proportion of streams that are monitored.
#' If there is no change, then this is the constant proportion among years.
#' @param nPop Number of subpopulations (indicator + non-indicator)
#' @param simYears Number of years to return monitoring coverage.
#' @param ppnChange The change in the proportion of streams that are monitored.
#' This can be positive or negative, but care should be taken that the 
#' \code{ppnSampled + ppnChange} is not less than zero or greater than one.
#' Otherwise, the function will return an error.
#' @param samplingDeclStart The year that the first change in sampling effort
#' began.
#' @param samplingDeclEnd The last year that had a change in sampling effort.
#' @param gen The number of different ages that fish return at. Used to 
#' distinguish the initialization years when calculating the rows in which 
#' \code{samplingDeclStart} and \code{samplingDeclEnd} should be applied.
#' @return Returns a matrix of zeroes and ones, with the number of rows equal 
#' to \code{simYears} and the number of columns equal to \code{nPop}. 
#' Zeroes indicate that subpopulation was not monitored in that year, wheras
#' ones indicate that the subpopulation was monitored in the given year.
#'
#' @examples
#' Example based on historic monitoring effort and recent declines reported for 
#' chum CUs on BC's central coast from the Pacific Salmon Explorer.
#' 
#' samplingDesign(ppnSampled = 0.85, nPop = 35, simYears = 57, 
#' ppnChange = -0.18, samplingDeclStart = 40, samplingDeclEnd = 50, gen = 5)
#'
#' @export

samplingDesign <- function(ppnSampled, nPop, simYears, ppnChange = 0, 
													 samplingDeclStart = NULL, 
													 samplingDeclEnd = samplingDeclStart, 
													 gen = NULL){
		
	if (ppnChange == 0) { # If ppnSampled is constant
		sampled <- matrix(sample(x = c(0,1), size = nPop*simYears, replace = TRUE, 
														 prob = c(1-ppnSampled, ppnSampled)),
											nrow = simYears, ncol = nPop)
	} else {
	
		# Calculate ppnSampled for each year
		ppnSampled.change <- rep(ppnSampled, simYears)
		
		# If decline happens over time, interpolate change:
		if ((samplingDeclEnd - samplingDeclStart) > 0) {
			ppnSampled.change[samplingDeclStart : samplingDeclEnd] <- 
				ppnSampled + approx(
					x = c(0, ppnChange), 
					n = samplingDeclEnd - samplingDeclStart + 2)$y[2:(samplingDeclEnd - samplingDeclStart + 2)]
		}
		
		# After change remain at new ppnSampled:
		ppnSampled.change[samplingDeclEnd : simYears] <- ppnSampled + ppnChange
		
		# If any ppnSampled is less than 0 or greater than 1, error!
		if(length(which(ppnSampled.change > 1 | ppnSampled.change < 0)) > 0) stop("Proportion of streams observed outside of (0,1).")
		
		# Compute matrix of sampling design (0 = not sampled) for each year
		# and subpopulation
		sampled <- matrix(NA, nrow = simYears, ncol = nPop)
		for(y in 1:simYears){
			sampled[y, ] <- sample(x = c(0, 1), size = nPop, replace = TRUE,
														 prob = c(1 - ppnSampled.change[y], ppnSampled.change[y]))
		} # end simYears
	} # end else (decline in monitoring)

	return(sampled)

} # end function