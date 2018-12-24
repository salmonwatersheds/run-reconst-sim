Evaluating the consequences of assumptions in run-reconstructions when assessing the biological status of salmon
================
Stephanie Peacock
2018-12-24

Model overview
--------------

This repo contains code that runs a stochastic simulation model of salmon population dynamics, allowing control over various biological and management factors that may influence the accuracy of status assessments. This model will be based on previous studies by Carrie Holt and colleagues (e.g., [Holt and Folkes 2015](http://dx.doi.org/10.1016/j.fishres.2015.01.002), Holt et al. 2018) that developed a simulation model comprised of sub-models for salmon population dynamics, observation of spawners, assessment, harvest, and performance (Fig. 1).

![Fig. 1. Schematic of the simulation model, including sub-models for population dynamics, harvest, observation, assessment, and performance. The entire process will be repeated for different autocorrelation in residuals among sub-populations, inter-annual variability in age-at-return, bias in harvest, and observation errors (see research questions, above). Adapted from [Holt et al. (2016)](https://www.psc.org/fund-project/adapting-benchmarks/).](model.png)

Model equations
---------------

### Population sub-model

The population dynamics of multiple sub-populations, *j*, designated as indicator or non-indicator streams, are simulated within a single hypothetical CU following a Ricker type stock-recruitment relationship, with parameters based on observations from central coast chum CUs ([Connors et al. 2018](https://salmonwatersheds.ca/library/lib_442/)). The true population dynamics are simulated in a loop over `nYears`, with calculations done in vectors across subpopulations.

The first loop is an initialization that calculates the recruits by brood year `recruitsBY`, *R*′<sub>*y*, *j*</sub>, for years 1 to `gen + 2`, where `gen` is the number of different ages fish can return at. These values are needed in order to calculate the first recruits by return year `recruitsRY`, which includes recruits that return as age 3, 4, or 5 year olds in the case of chum salmon (the model is flexible to incorporate the possible of 2 and 6 year olds returning too). For each year in this initialization, we assumed that the number of spawners was equal to 20% of equilibrium spawner abundance, *S*<sub>*j*</sub><sup>\*</sup> = *a*<sub>*j*</sub>/*b* for subpopulation *j* (Holt et al. 2018 CSAS).

The second loop simulates the true population dynamics from year `gen + 3` to `nYears`. The number of salmon in **return year *t*** and subpopulation *j*, *R*<sub>*t*, *j*</sub>, is calculated as:

*R*<sub>*t*,  *j*</sub> = *R*′<sub>*t* − 3,  *j*</sub> *p*<sub>*t* − 3,  3</sub> + *R*′<sub>*t* − 4,  *j*</sub> *p*<sub>*t* − 4,  4</sub> + *R*′<sub>*t* − 5,  *j*</sub> *p*<sub>*t* − 5,  5</sub>

where *p*<sub>*y*, *g*</sub> is the proportion of recruits from **brood year *y*** returning as *g* year olds. Note that we assume that the proportion of recruits returning at a given age is the same among subpopulations, but incorporate interannual variability as in Holt et al. (2018 CSAS):

$$ p\_{y,g} = \\frac{\\bar{p}\_g \\! \\exp ({\\bar{\\omega} \\, \\varepsilon\_{y,g}})}{\\sum\_{G = 3}^{5} \\bar{p}\_G \\! \\exp ({\\bar{\\omega} \\, \\varepsilon\_{y,G}})} $$

In the model, this is coded as:

    recruitsRY[y, ] <- ppnAge[cbind(y - ages, 1:simPar$gen)] %*% recruitsBY[y - ages,]

where the `ppnAge` matrix incorporating natural interannual variability is calculated prior to the population dynamics loop using the `ppnAgeErr` function.

The number of spawners **returning in year *t*** is the number of returning salmon *R*<sub>*t*, *j*</sub> times 1 − *h*<sub>*t*</sub>, where *h*<sub>*t*</sub> is the realized harvest rate for year *t*. We assume a predetermined target harvest rate, *h*′<sub>*t*</sub>, but incorporate normally distributed error around that:

*h*<sub>*t*</sub> ∼ *N*(*h*′<sub>*t*</sub>,  *σ*<sub>*h*</sub>)

with the constraint that 0 &lt; *h*<sub>*t*</sub> &lt; 1. If *h*<sub>*t*</sub> ≤ 0 or *h*<sub>*t*</sub> ≥ 1, we drew another random *h*<sub>*t*</sub> until this condition was satisdied. As with the `ppnAge`, we assume the same harvest rate among subpopulatins and the vector of harvest rates `harvestRate` is calculated prior to the population dynamics loop and then applied within.

The number of fish **returning to spawn in year *t*** and subpopulation *j* is calculated as:

*S*<sub>*t*, *j*</sub> = (1 − *h*<sub>*t*</sub>) *R*<sub>*t*, *j*</sub>

or in R code:

    spawners[y, ] <- harvestRate[y] * recruitsRY[y, ]

The true total catch of fish that would have returned to streams within the CU is calculated as:

*C*<sub>*t*</sub> = (1 − *h*<sub>*t*</sub>) ∑<sub>*j*</sub>*R*<sub>*t*, *j*</sub>

For this project, we ignore potential straying of returning adults among subpopulations that was included in some previous analyses (e.g., Holt et al. 2018 and [Peacock & Holt 2012](http://www.nrcresearchpress.com/doi/full/10.1139/f2012-004)).

Finally, to calculate the number of recruits from **brood year *y*** and subpopulation *j*, we apply the Ricker model:

*R*′<sub>*y*, *j*</sub> = *S*<sub>*y*, *j*</sub> exp(*a*<sub>*j*</sub> − *b* *S*<sub>*y*, *j*</sub>) exp(*ϕ*<sub>*y*, *j*</sub>)

where *a*<sub>*j*</sub> is the log recruits per spawner at low spawner abundance (i.e., productivity), which is assumed to be normally distributed among subpopulations with some mean $\\bar{a}$ and variance *σ*<sub>*a*</sub><sup>2</sup>, density-dependence parameter *b*, which is assumed to be the same among all subpopulations, and *ϕ*<sub>*y*, *j*</sub> are the recruitment deviations for year *y* and subpopulation *j*. We incorporated temporal autocorrelation in recruitment residuals:

*ϕ*<sub>*y*, *j*</sub> = *ρ* *ϕ*<sub>*y* − 1, *j*</sub> + *υ*<sub>*y*, *j*</sub>

where *ρ* is the temporal autocorrelation coefficient and *υ*<sub>*y*, *j*</sub> is drawn from a multivariate normal distribution with means zero (\*\*right? I think Cam mentioned they are no longer using the $ - \_^2 / 2 $?\*\*) and variance-covariance matrix:

$$ \\Sigma\_{j \\times j} = \\left\[ \\begin{array} 
 \\sigma\_\\upsilon^2 & \\rho\_\\upsilon \\, \\sigma\_\\upsilon^2 & \\ldots & \\rho\_\\upsilon \\, \\sigma\_\\upsilon^2 \\\\
\\rho\_\\upsilon \\, \\sigma\_\\upsilon^2 & \\sigma\_\\upsilon^2 & \\ldots & \\rho\_\\upsilon \\, \\sigma\_\\upsilon^2 \\\\
\\vdots & \\vdots & \\ddots & \\vdots \\\\
\\rho\_\\upsilon \\, \\sigma\_\\upsilon^2 & \\rho\_\\upsilon \\, \\sigma\_\\upsilon^2 & \\ldots & \\sigma\_\\upsilon^2 \\\\ 
\\end{array} \\right\]\_{j \\times j} $$

**Question: I think *ρ*<sub>*υ*</sub> (autocorrelation among subpopulations) is different from *ρ* (temporal autocorrelation), but this is not clear in Holt et al. (2018 CSAS) equations (F7). Check.** Here, *σ*<sub>*υ*</sub> is the standard deviation in residuals without autocorrelation (Rickery 1975, Holt and Bradford 2011) and ***ρ*<sub>*υ*</sub> is the "spatial" autocorrelation among subpopulations.** In the initialization loop, we assumed that *ϕ*<sub>*y* = 1, *j*</sub> = 0.

The calculation of recruits is done in the model using the `rickerModel` function, and the resulting recruits and recruitment deviations are stored in matrices `recruitsBY` and `phi`, respectively.

### Observation sub-model

Spawner abundances are observed with probability, *p*<sub>*y*, *j*</sub>. Depending on the scenario, this probability may be constant over time, or incorporate some change in monitoring effort (e.g., a decline in the probability of being sampled at some point in time or over a period of time). **For now, I have assumed that the sampling probability is the same for indicator and non-indicator streams, but this may need to be changed. Note that the observed spawner abundance for non-indicator streams is not used in the calculation of aggregate spawner abundance to the CU except for the calculation of the second expansion factor. ** In R, this sampling design is applied by the function `samplingDesign`, which takes a base value for *p*, as well as a change in *p* and a start and end year over which that change is applied. Observed spawner abundance incorporates log-normal observation error:

$$ \\hat{S}\_{y,j} = z\_{y,j} \\; \[ S\_{y,j} \\; \\exp (\\delta\_{y,j}) \] $$

where *z*<sub>*y*, *j*</sub> ∼ Bernoulli(prob = *p*<sub>*y*, *j*</sub>), *δ*<sub>*y*, *j*</sub> ∼ *N*(0, *σ*<sub>*δ*</sub><sup>2</sup>), and *σ*<sub>*δ*</sub> is the standard deviation in observation error of spawner abundances.

The observed catch to the entire CU in return year *t* is observed with log-normal error:

$$ \\hat{C}\_{t} = C\_{t} \\; \\exp (\\chi\_{t}) $$

where *χ*<sub>*t*</sub> ∼ *N*(0, *σ*<sub>*χ*</sub><sup>2</sup>), and *σ*<sub>*χ*</sub> is the standard deviation in catch error. **If we wanted to incorporate a bias in catch** to simulate a scenario where 1. fish are caught from other CUs, or 2. fish from the focal CU that were caught in other fisheries, then we could introduce a positive or negative mean in catch error, but for now I have left as zero.

For this project, we are assuming average age-at-return is applied in run reconstruction, but this "observed age-at-return" may have error associated with it. Previous models (e.g., Holt et al. 2018 CSAS) have included error in the "estimated age-at-return" for each return year, where the mean is the true age-at-return, but for the central coast chum, annual age-at-return data are rarely available and so the average is used. We apply observation error to get an average age-at-return that is applied in run reconstruction:

$$ \\hat{pr}\_{g} = \\frac{\\bar{p}\_g \\! \\exp ({\\bar{\\omega}\_{pr} \\, \\varepsilon\_{g}})}{\\sum\_{G = 3}^{5} \\bar{p}\_G \\! \\exp ({\\bar{\\omega}\_{pr} \\, \\varepsilon\_{G}})} $$

where $\\bar{p}\_g$ is the average age-at-return applied in the population submodel. **There are several other possible approaches here that might make more sense:** 1. minor variation: use the average of the true ages `ppnAge` as the means; 2. sub-sample from the true porportion `ppnAge` in some way that reflects the actual number of years that are used in these calculations. 3. sub-sample from the true porportion `ppnAge` AND apply observation error, then calculate mean. **Discuss with Cam.**

### Assessment sub-model

#### Expansion factors for run reconstruction

For stock-recruit benchmarks, time-series of aggregate spawners to the CU and catch are needed in order to calculate recruitment for a given brood year. For chum salmon that can return as 3, 4, or 5 year-olds, these time-series need to be continuous as one year of missing data will result in multiple years of missing stock-recruit pairs. In order to calculate harvest rates from the total onserved catch, $\\hat{C}\_{t}$, an estimate of the total number of spawners is required from indicator streams, non-indicator streams, and streams not monitored. To get this, the observed number of spawners from indicator streams is multiplied by three expansion factors. In the model, each of these expansion factors is calculated and applied separately using functions derived from LGL's North & Central Coast Salmon Database (NCCSDB) package (Challenger 2018).
