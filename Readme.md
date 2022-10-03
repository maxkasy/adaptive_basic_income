# Adaptive maximization of welfare with negative income taxes

This archive contains code to implement our algorithm for maximization of social welfare in an adaptive negative income tax experiment.
This is work in progress.

## 0_prepare_soep.R
Preparation of SOEP data.
Extraction of wage, participation indicator, and covariates (years of education, gender, relationship status).

## 1_maximum_likelihood_MCMC_labor_exponential.R
This is the core file for estimation of the structural labor supply model.

This file defines the functions *log_likelihood* (corresponding to our model), *log_posterior* (adding log prior terms for regularization), *run_metropolis* (MCMC algorithm to sample from the posterior), *parallel_metropolis* (parallelizing MCMC), and *metropolis_master* (initializing MCMC by finding the maximum a posteriori and the corresponding Hessian), in addition to several helper functions.

A sample from the posterior is saved in the global variable *use_chain*.

## 2_labor_supply_simulate_functions.R
This is the core file for welfare evaluations.
It sets the normative parameters *welfare_exponent* and *break_even_u*, and the set of policies to consider *sim_policies*.

*simulate_data* simulates from the model, taking as arguments model parameters and a matrix of covariate values.
*simulate_sample* simulates data for the different policies.
*counterfactual_welfare* aggregates these data for welfare evaluation of the alternative policies.
*thompson_probabilities* repeats this for multiple draws of the model parameters coming from the mcmc chain, and returns the posterior probability for each policy to be optimal.

## 3_labor_supply_test.R
This file tests simulation of data and estimation of the structural model with fictitious parameter values.

## 4_soep_calibration.R
This file fits the structural model to the SOEP data.
The parameters *beta*, *tau* and *rho* are not identified without policy variation; *log_posterior* is redefined to regularize these parameters.

