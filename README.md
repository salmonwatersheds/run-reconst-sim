# Evaluating the consequences of assumptions in run-reconstructions when assessing the biological status of salmon
#### Stephanie J. Peacock <stephanie.j.peacock@gmail.com>
#### December 21, 2018

This repo contains R code related to the PSF project "Evaluating the consequences of assumptions in run-reconstructions when assessing the biological status of salmon"

## Model overview

This project will develop and apply a stochastic simulation model of salmon population dynamics that allows control over various biological and management factors that may influence the accuracy of status assessments.  This model will be based on previous studies by Carrie Holt and colleagues (e.g., [Holt and Folkes 2015](http://dx.doi.org/10.1016/j.fishres.2015.01.002), Holt et al. 2018) that developed a simulation model comprised of sub-models for salmon population dynamics, observation of spawners, assessment, harvest, and performance (Fig. 1).

![Fig. 1. Schematic of the simulation model, including sub-models for population dynamics, harvest, observation, assessment, and performance.  The entire process will be repeated for different autocorrelation in residuals among sub-populations, inter-annual variability in age-at-return, bias in harvest, and observation errors (see research questions, above). Adapted from [Holt et al. (2016)](https://www.psc.org/fund-project/adapting-benchmarks/).](model.png)

## Model equations

### Population sub-model
The population dynamics of multiple sub-populations, $j$, designated as indicator or non-indicator streams, are simulated within a single hypothetical CU following a Ricker type stock-recruitment relationship, with parameters based on observations from central coast chum CUs ([Connors et al. 2018](https://salmonwatersheds.ca/library/lib_442/)). The true population dynamics are simulated in a loop over `nYears`, with calculations done in vectors across subpopulations.

The first loop is an initialization that calculates the recruits by brood year `recruitsBY`, $R'_{y,j}$,  for years 1 to `gen + 2`, where `gen` is the number of different ages fish can return at. These values are needed in order to calculate the first recruits by return year `recruitsRY`, which includes recruits that return as age 3, 4, or 5 year olds in the case of chum salmon (the model is flexible to incorporate the possible of 2 and 6 year olds returning too). For each year in this initialization, we assumed that the number of spawners was equal to , 20% of equilibrium spawner abundance, $S_{j}^* = a_j/b$ for subpopulation $j$ (Holt et al. 2018 CSAS).

The second loop simulates the true population dynamics from year `gen + 3` to `nYears`. The number of salmon returning for a given year $y$ and subpopulation $j$, $R_{y,j}$, is calculated as:
$$ R_{y,j} = R'_{y-3,j} \dot p_{y-3, 3} +  R'_{y-4,j} \dot p_{y-4, 4} + R'_{y-5,j} \dot p_{ y-5, 5}.$$

where $p_{y,g}$ is the proportion of recruits from brood year $y$ returning as $g$ year olds. Note that we assume that the proportion of recruits returning at a given age is the same among subpopulations, but incorporate interannual variability as in Holt et al. (2018 CSAS):
$$ p_{y,g} = \frac{\bar{p}_g e^{\bar{\omega}*\epsilon_{y,g}}}{\sum_{G = 3}^5 \bar{p}_G e^{\bar{\omega}*\epsilon_{y,G}}} $$

In the model, this is coded as:
```
recruitsRY[y, ] <- ppnAge[cbind(y - ages, 1:simPar$gen)] %*% recruitsBY[y - ages,]
```
where the `ppnAge` matrix incorporating natural interannual variability is calculated prior to the population dynamics loop using the `ppnAgeErr` function.

The number of spawners returning in year $y$ is the number of returning salmon $R_{y,j}$ times $1-h_y$, where $h_y$ is the realized harvest rate for year $y$. We assume a predetermined target harvest rate, $h'_y$, but incorporate normally distirbuted error around that:
$$ h_y \sim N(h'_y, \sigma_h) $$
with the constraint that $0<h_y<1$. If $h_y \leq 0$ or $h_y \geq 1$, we drew another random $h_y$ until this condition was satisdied. As with the `ppnAge`, we assume the same harvest rate among subpopulatins and the vector of harvest rates `harvestRate` is calculated prior to the population dynamics loop and then applied within. 

The number of fish returning to spawn in year $y$ and subpopulation $j$ is calculated as:
$$ S_{y,j} = R_{y,j} \dot (1- h_y) $$
or in R code:
```
spawners[y, ] <- harvestRate[y] * recruitsRY[y, ]
```
For this project, we ignore potential straying of returning adults among subpopulations that was included in some previous analyses (e.g., Holt et al. 2018 and [Peacock & Holt 2012](http://www.nrcresearchpress.com/doi/full/10.1139/f2012-004)).

Finally, to calculate the number of recruits from brood year $y$ and subpopulation $j$, we apply the Ricker model:
$$ R'_{y,j} = S_{y,j} \exp{a_j - b S_{y,j}} \dot \exp{phi_y,j} $$
where $a_j$ is the log recruits per spawner at low spawner abundance (i.e., productivity), which is assumed to be normally distributed among subpopulations with some mean $\bar{a}$ and variance $\sigma_a^2$, density-dependence parameter $b$, which is assumed to be the same among all subpopulations, and $\phi_{y,j}$ are the recruitment deviations for year $y$ and subpopulation $j$. We incorporated temporal autocorrelation in recruitment residuals:
$$ \phi_{y,j} = \rho \dot \phi_{y-1,j} + \upsilon_{y,j} $$
where $\rho$ is the temporal autocorrelation coefficient and $\upsilon_{y,j}$ is drawn from a multivariate normal distribution with mean $-\sigma_\upsilon^2 / 2$, to make the arithmetic mean of the lognormally distributed recruitment variation equal to 1 (as in Holt and Peterman 2008), and variance-covariance matrix
$$ \Sigma_{j \times j} = \left\[ \begin{array} 
\sigma_\upsilon^2 & \rho_\upsilon \sigma_\upsilon^2 & \ldots & \rho_\upsilon \sigma_\upsilon^2 \\
\rho_\upsilon \sigma_\upsilon^2 & \sigma_\upsilon^2 & \ldots & \rho_\upsilon \sigma_\upsilon^2 \\
\vdots & \vdots & \ddots & vdots \\
\rho_\upsilon \sigma_\upsilon^2 & \rho_\upsilon \sigma_\upsilon^2 & \ldots & \sigma_\upsilon^2 & 
\end{array} \right\] $$ 
*Question: I think $\rho_\upsilon$ (autocorrelation among subpopulations) is different from $\rho$ (temporal autocorrelation), but this is not clear in Holt et al. (2018 CSAS) equations (F7). Check.*
Here, $\sigma_\upsilon$ is the standard deviation in residuals without autocorrelation (Rickery 1975, Holt and Bradford 2011) and *$\rho_\upsilon$ is the "spatial" autocorrelation among subpopulations.

The calculation of recruits is done in the model using the `rickerModel` function, and the resulting recruits and recruitment deviations are stored in matrices `recruitsBY` and `phi`, respectively.

### Observation sub-model
