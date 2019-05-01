Evaluating the consequences of commmon assumptions in run-reconstructions when quantifying the biological status of Pacific salmon
================
Stephanie Peacock
2019-05-01

Motivation
==========

The [Pacific Salmon Foundation (PSF)](http://www.psf.ca) has recently led efforts to quantify the biological status of Pacific salmon Conservation Units (CUs) on British Columbia’s (BC) north and central coast, following the approaches outlined under Canada’s Wild Salmon Policy. In this work, the PSF has evaluated current biological status using two metrics and associated benchmarks: 1) historic spawner abundance and 2) the shape of the stock-recruitment relationship ([Connors et al. 2018](https://salmonwatersheds.ca/library/lib_442/)).

In order to evaluate biological status using stock-recruitment benchmarks, a continuous time-series of spawner abundance and recruitment is needed for each CU. However, not every stream within a CU is monitored each year and, in most cases, run reconstruction analyses must be undertaken in order to estimate the total number of spawners within a CU, determine harvest rates, and generate CU-level estimates of recruitment. The PSF has contracted [LGL Ltd.](http://www.lgl.com) for the past 10 years to work with north coast Fisheries and Oceans Canada (DFO) stock assessment biologists to reconstruct estimates of spawner abundance and recruitment for each salmon CU in the region ([English et al. 2012](https://salmonwatersheds.ca/library/lib_1/), [2016](https://salmonwatersheds.ca/library/lib_435/)).

The approach to generating these run reconstructions varies among species, and among CUs within a species. In the most basic case, the total escapement to the CU is obtained by expanding the observed number of spawners in streams that are monitored to account for: (I) indicator streams that are not monitored in a given year, (II) non-indicator streams that are rarely monitored, and (III) observer efficiency. In order to estimate harvest rates, assumptions about the migration routes and run-timing of harvested salmon are made to assign catch – recorded at the level of Statistical Areas (SAs) – to different CUs. Finally, assumptions about the age-at-return of spawners in any given year may be required to calculate recruitment. The combined influence of these assumptions on our ability to accurately assess the status of CUs is unknown. The objective of this project is to understand the consequences of basic assumptions in run reconstructions surrounding expansion factors, harvest rates, and age-at-return on status assessments and the circumstances under which they may bias the kinds of assessments of biological status undertaken by the PSF.

Assumptions in run reconstruction
---------------------------------

The specific assumptions made in reconstructing spawner, harvest, and resulting recruitment time-series may differ among species and even among CUs within a species due to differences in, for example, the spatial resolution of catch data, run-timing variability, or additional data such as from CWTs (see Appendix E of [English et al. (2016)](https://salmonwatersheds.ca/library/lib_435/) for details). At the most basic level, as a starting point for further analysis, we will consider a generic central coast chum CU where harvest rules are relatively simple and a constant age-at-return is assumed. The main assumptions in the run reconstruction that this project will investigate are:

1.  The relative contribution of each indicator stream to aggregate spawner abundance for the CU does not vary through time. Under this assumption, observed spawner abundance is multiplied by Expansion Factor I to account for indicator streams that were not monitored in a given year. Expansion Factor I is calculated based on the current decade or, if an indicator stream is not monitored in that decade, on the closest decade that has at least one escapement estimate from each indicator stream.

2.  The relative contribution of each non-indicator stream to aggregate spawner abundance for the CU does not vary through time. Under this assumption, the estimated spawner abundance to all indicator-streams is multiplied by Expansion Factor II to account for non-indicator streams. Expansion Factor II is calculated based on the decade with the best survey coverage (indicator and non-indicator) for that CU.

3.  Observer efficiency is constant over time for a given CU. Under this assumption, the estimated spawner abundance in both indicator and non-indicator streams is multiplied by Expansion Factor III to account for imperfect observer efficiency. Expansion Factor III is based on expert opinion of regional DFO staff familiar with escapement monitoring techniques and is constant over time within a CU.

4.  Salmon harvested within a SA were destined to spawn in streams within that SA in proportion to the escapement to each stream, and salmon returning to a SA were not caught in other SAs. Under this assumption, harvest rates are calculated from catch within a single SA and estimated total escapement to the CU. Violation of this assumption would result in over- or under-estimation of harvest rates.

5.  Average age composition based on available data is representative of year-specific age composition of returns.

The project will quantify the consequences of these assumptions for the bias and precision of benchmarks and, ultimately, status assessments under a base case that corresponds to the historical conditions for central coast chum (e.g., historical monitoring coverage, trends in covariance among subpopulations, and changes in harvest rates).

Files and folders
=================

model
-----

All functions to run a simulation are in the `model` folder, along with a complete description of the model equations. Functions are documented in `roxygen2` style to facilitate conversion to an R package at some point. The main function that runs a single MC simulation of salmon population dynamics, observation, and assessment, and returns the observed and true status as well as performance metrics is the `reconstSim` function described in the `reconstrSimulator.R` file.

The functions called by `reconstrSim` that comprise the simulation model are organized into files based on the submodel that they correspond to:

-   `populationSubmodFns.R` constains functions
    -   `rickerModel` that simulates the Ricker model
    -   `ppnAgeErr` that simulates the proportion of fish returning at ages 2-6 (by brood year) with error
-   `obsSubmodFns.R` contains functions
    -   `samplingDesign` that simulates the partial monitoring of indicator and non-indicator streams within the CU, with the option of incorporating change in monitoring effort over time
-   `expansionFactors.R` contains functions
    -   `refDecade` that selects a reference decade from which to calculate relevant expansion factors (as per `ExpFactor1RefDecade` from the NCCSDB packaage)
    -   `ExpFactor1` that calculates the value of Expansion Factor 1 for each year and decade
    -   `ExpFactor2` that calculates the value of Expansion Factor 2 for each decade
-   `benchmarkFns.R` contains functions
    -   `calcSmsy` that calculates *S*<sub>*M**S**Y*</sub> (i.e., upper SR benchmark) given Ricker parameters, using the explicit solution from Scheuerell (2016)
    -   `Sgen.optim` that calculates the likelihood of residuals between the projected recruits from an estimated *S*<sub>*G**E**N*1</sub> (`Sgen.hat`) and a known value of *S*<sub>*M**S**Y*</sub>, to be used in the optimization of *S*<sub>*G**E**N*1</sub>.
    -   `calcSgen` that calculates *S*<sub>*G**E**N*1</sub> by optimizing `Sgen.optim`
    -   `assessMetric` that returns the status given current abundance and upper and lower benchmarks on a given metric
    -   `assessPop` that returns a population assessment based on both SR and percentile metrics of spawners abundance, to be applied to either observed or true spawners and recruits
    -   `assessTruePop` that returns a population assessment based on the "true" spawner abundance and the "true" SR benchmarks as calculated from the underlynig SR parameters.
-   `performanceFns.R` contains function
    -   `perfStatus` that takes the true and observed status output from `assessPop` and returns the bias (raw mean error) in benchmarks (observed - true) and the status code 1-9, which tells whether the final status (green (1), amber (2), or red (3)) was the same for true and observed or different, and if they were different, how they differed.

Details of all functions can be found in their respective files.

runSims
-------

This folder contains the code to run simulations specific to the run-reconstruction project and reproduce output described in the associated paper. It may be useful as a reference for future studies, and could be altered to explore sensitivity of assessments to additional model parameters. This folder contains the following files:

-   `runSim.R` produces the three base-case results and loops through the three base cases, performing all sensitivity analyses under each base case. The output of each sensitivity analysis is saved as a RDS file that can be loaded for plotting.
-   `runSensitivity.R` contains three functions that run multiple MC simulations over a range of parameter values for sensitivity analyses. It contains the following functions:
    -   `makeParList` takes a range of values for a single parameter and returns a list of parameter sets with the focal parameter changed in each element. Useful for simple univariate sensitivity analyses.
    -   `runSensitivity` takes a list of parameter sets (can be from `makeParList` above) and runs multiple MC simulations of `reconstrSim()` for each parameter set (i.e., element) in that list.
    -   `delistSensitivity` delists the output from `runSensitivity` so that RB and proportion of simulations with correct or misclassified status can be easily plotted over the different parameter sets.
-   `nSims.R` constains code to run 10 000 MC simulations and determine the number of simulations necessary to ensure that the error in performance metrics is less than 3%.
-   `figures.R` contains code to source the RDS output from `runSim.R` and plot the figures included in the paper.

data
----

This folder contains the data used to parameterize the simulation model and the full suite of base case parameters applied.

-   `baseSimPar.csv` contains the base-case parameters applied in simulations.
-   `ageAtReturn.R` reads in data from the NCCSDBV2 (`nccdbv2_age-composition-data.csv`) and calculates average for central coast chum.
-   `rickerParams.R` reads in river-level stock-recruitment data for north and central coast chum salmon (`NCC_chum_streams_SR_data.csv`) and calculates Ricker parameters for simulating true subpopulation dynamics based on central coast chum.
-   `monitoring.R` looks at the four real monitoring coverage scenarios (A-D) calling on data from `English2016_monitoringCoverage.csv` and `NCC_chum_streams_SR_data.csv` and plots the coverage over time under these different scenarios.
-   `harvestRate.R` takes data on exploitation rates downloaded from the Pacific Salmon Explorer (`PSE_chumExploitationRate.csv`) and plots trends over time AND calculates parameters of the simple harvest control rule using data from the Pacific Salmon Explorer on total run size (`PSE_chumSpawnerAbundanceAndTotalRunSize.csv`) and from the NCCSDBV2 on escapement (`nccdbv2_NCCStreamEscapement.csv`)
