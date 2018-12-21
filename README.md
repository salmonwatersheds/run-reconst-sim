Evaluating the consequences of assumptions in run-reconstructions when assessing the biological status of salmon
================
Stephanie Peacock
2018-12-21

Model overview
--------------

This repo contains code that runs a stochastic simulation model of salmon population dynamics, allowing control over various biological and management factors that may influence the accuracy of status assessments. This model will be based on previous studies by Carrie Holt and colleagues (e.g., [Holt and Folkes 2015](http://dx.doi.org/10.1016/j.fishres.2015.01.002), Holt et al. 2018) that developed a simulation model comprised of sub-models for salmon population dynamics, observation of spawners, assessment, harvest, and performance (Fig. 1).

![Fig. 1. Schematic of the simulation model, including sub-models for population dynamics, harvest, observation, assessment, and performance. The entire process will be repeated for different autocorrelation in residuals among sub-populations, inter-annual variability in age-at-return, bias in harvest, and observation errors (see research questions, above). Adapted from [Holt et al. (2016)](https://www.psc.org/fund-project/adapting-benchmarks/).](model.png)

Model equations
---------------

### Population sub-model

The population dynamics of multiple sub-populations, *j*, designated as indicator or non-indicator streams, are simulated within a single hypothetical CU following a Ricker type stock-recruitment relationship, with parameters based on observations from central coast chum CUs ([Connors et al. 2018](https://salmonwatersheds.ca/library/lib_442/)). The true population dynamics are simulated in a loop over `nYears`, with calculations done in vectors across subpopulations.

The first loop is an initialization that calculates the recruits by brood year `recruitsBY`, *R*′<sub>*y*, *j*</sub>, for years 1 to `gen + 2`, where `gen` is the number of different ages fish can return at. These values are needed in order to calculate the first recruits by return year `recruitsRY`, which includes recruits that return as age 3, 4, or 5 year olds in the case of chum salmon (the model is flexible to incorporate the possible of 2 and 6 year olds returning too). For each year in this initialization, we assumed that the number of spawners was equal to , 20% of equilibrium spawner abundance, *S*<sub>*j*</sub><sup>\*</sup> = *a*<sub>*j*</sub>/*b* for subpopulation *j* (Holt et al. 2018 CSAS).

The second loop simulates the true population dynamics from year `gen + 3` to `nYears`. The number of salmon returning for a given year *y* and subpopulation *j*, *R*<sub>*y*, *j*</sub>, is calculated as:

*R*<sub>*y*,  *j*</sub> = *R*′<sub>*y* − 3,  *j*</sub> *p*<sub>*y* − 3,  3</sub> + *R*′<sub>*y* − 4,  *j*</sub> *p*<sub>*y* − 4,  4</sub> + *R*′<sub>*y* − 5,  *j*</sub> *p*<sub>*y* − 5,  5</sub>

where *p*<sub>*y*, *g*</sub> is the proportion of recruits from brood year *y* returning as *g* year olds. Note that we assume that the proportion of recruits returning at a given age is the same among subpopulations, but incorporate interannual variability as in Holt et al. (2018 CSAS):

$$ p\_{y,g} = \\frac{\\bar{p}\_g \\! \\exp ({\\bar{\\omega} \\, \\varepsilon\_{y,g}})}{\\sum\_{G = 3}^{5} \\bar{p}\_G \\! \\exp ({\\bar{\\omega} \\, \\varepsilon\_{y,G}})} $$

In the model, this is coded as:

    recruitsRY[y, ] <- ppnAge[cbind(y - ages, 1:simPar$gen)] %*% recruitsBY[y - ages,]

where the `ppnAge` matrix incorporating natural interannual variability is calculated prior to the population dynamics loop using the `ppnAgeErr` function.

The number of spawners returning in year *y* is the number of returning salmon *R*<sub>*y*, *j*</sub> times 1 − *h*<sub>*y*</sub>, where *h*<sub>*y*</sub> is the realized harvest rate for year *y*. We assume a predetermined target harvest rate, *h*′<sub>*y*</sub>, but incorporate normally distirbuted error around that:

*h*<sub>*y*</sub> ∼ *N*(*h*′<sub>*y*</sub>,  *σ*<sub>*h*</sub>)

with the constraint that 0 &lt; *h*<sub>*y*</sub> &lt; 1. If *h*<sub>*y*</sub> ≤ 0 or *h*<sub>*y*</sub> ≥ 1, we drew another random *h*<sub>*y*</sub> until this condition was satisdied. As with the `ppnAge`, we assume the same harvest rate among subpopulatins and the vector of harvest rates `harvestRate` is calculated prior to the population dynamics loop and then applied within.

The number of fish returning to spawn in year *y* and subpopulation *j* is calculated as:

*S*<sub>*y*, *j*</sub> = *R*<sub>*y*, *j*</sub> (1 − *h*<sub>*y*</sub>)

or in R code:

    spawners[y, ] <- harvestRate[y] * recruitsRY[y, ]

For this project, we ignore potential straying of returning adults among subpopulations that was included in some previous analyses (e.g., Holt et al. 2018 and [Peacock & Holt 2012](http://www.nrcresearchpress.com/doi/full/10.1139/f2012-004)).

Finally, to calculate the number of recruits from brood year *y* and subpopulation *j*, we apply the Ricker model:

*R*′<sub>*y*, *j*</sub> = *S*<sub>*y*, *j*</sub> exp(*a*<sub>*j*</sub> − *b* *S*<sub>*y*, *j*</sub>) exp(*ϕ*<sub>*y*, *j*</sub>)

where *a*<sub>*j*</sub> is the log recruits per spawner at low spawner abundance (i.e., productivity), which is assumed to be normally distributed among subpopulations with some mean $\\bar{a}$ and variance *σ*<sub>*a*</sub><sup>2</sup>, density-dependence parameter *b*, which is assumed to be the same among all subpopulations, and *ϕ*<sub>*y*, *j*</sub> are the recruitment deviations for year *y* and subpopulation *j*. We incorporated temporal autocorrelation in recruitment residuals:

*ϕ*<sub>*y*, *j*</sub> = *ρ* *ϕ*<sub>*y* − 1, *j*</sub> + *υ*<sub>*y*, *j*</sub>

where *ρ* is the temporal autocorrelation coefficient and *υ*<sub>*y*, *j*</sub> is drawn from a multivariate normal distribution with mean $ - \_^2 / 2$, to make the arithmetic mean of the lognormally distributed recruitment variation equal to 1 (as in Holt and Peterman 2008), and variance-covariance matrix:

$$ \\Sigma\_{j \\times j} = \\left\[ \\begin{array} 
 \\sigma\_\\upsilon^2 & \\rho\_\\upsilon \\, \\sigma\_\\upsilon^2 & \\ldots & \\rho\_\\upsilon \\, \\sigma\_\\upsilon^2 \\\\
\\rho\_\\upsilon \\, \\sigma\_\\upsilon^2 & \\sigma\_\\upsilon^2 & \\ldots & \\rho\_\\upsilon \\, \\sigma\_\\upsilon^2 \\\\
\\vdots & \\vdots & \\ddots & \\vdots \\\\
\\rho\_\\upsilon \\, \\sigma\_\\upsilon^2 & \\rho\_\\upsilon \\, \\sigma\_\\upsilon^2 & \\ldots & \\sigma\_\\upsilon^2 \\\\ 
\\end{array} \\right\]\_{j \\times j} $$

**Question: I think *ρ*<sub>*υ*</sub> (autocorrelation among subpopulations) is different from *ρ* (temporal autocorrelation), but this is not clear in Holt et al. (2018 CSAS) equations (F7). Check.** Here, *σ*<sub>*υ*</sub> is the standard deviation in residuals without autocorrelation (Rickery 1975, Holt and Bradford 2011) and ***ρ*<sub>*υ*</sub> is the "spatial" autocorrelation among subpopulations.**

The calculation of recruits is done in the model using the `rickerModel` function, and the resulting recruits and recruitment deviations are stored in matrices `recruitsBY` and `phi`, respectively.

### Observation sub-model
