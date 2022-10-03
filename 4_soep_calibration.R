library(tidyverse)
library(patchwork)

source("1_maximum_likelihood_MCMC_labor_exponential.R")
source("2_labor_supply_simulate_functions.R")

soep_data = read_csv("Data_processed/SOEP_with_income.csv") 
    # sample_n(4000) # subsample for faster initial calibration
dim_x = 4

# Load estimates from previous calibration round
load("Data_processed/calibrated_estimates.RDat")
start_value = unlist(estimates$Bayes) # starting value for MLE
dim_par = length(unlist(start_value)) - 2 # last two entries are log posterior and log likelihood
start_value = start_value[1:dim_par]
true_par = start_value |> # convert to list format 
    par_vector_to_list()



# Fit model to data ----

# construct informative prior
log_posterior = function(par_vector){
    log_likelihood_vectorarg(par_vector) + # Add log prior for regularization
        + 2 * (log(max(par_vector[3],0)) - 
                   log(max(par_vector[2],0))) + #- # uniform pseudo-prior for log sigma
        dnorm(par_vector[4], sd = .1, log = T) + # informative prior for optimization error
    dnorm(par_vector[1], mean = 0.2, sd = .01, log = T) + # informative prior for beta
    dnorm(par_vector[3], mean = .001, sd = .0001, log = T) # informative prior for tau
}


# Run calibration
if (F){
    metropolis_master(soep_data, dim_x, start_value)
    MCMC_diagnostics()
}

# Simulate data from calibrated values ----

par = estimates$Bayes[1:dim_par] |> 
    par_vector_to_list()

sim_data = simulate_sample(x = x_soep, par = par, n_per_policy = 4000)
sim_data |>
    filter(y < 5000) |>
    counterfactual_plots(filename = "Figures/soep_simulated_plots.png")

# Compare simulated and actual data


baseline_sim_data = sim_data |> 
    filter(y < 5000, w==1, y0==0) 

baseline_soep_data = soep_data |>
    filter(y < 5000) 

p1 = baseline_sim_data |> 
    ggplot(aes(x=y)) + geom_histogram(bins = 20) +
    theme_minimal() + labs(title = "Calibrated")

p2 = baseline_soep_data |> 
    ggplot(aes(x=y)) + geom_histogram(bins = 20) +
    theme_minimal() + labs(title = "SOEP")

ggsave("Figures/data_vs_calibration.png", 
       p1 / p2, width = 5, height = 8)

summary(baseline_sim_data$y)
summary(baseline_soep_data$y)

# Counterfactual welfare, based on estimates ----
welfare =
    estimates |> 
    slice(1:dim_par) |>
    pull(Bayes) |> 
    par_vector_to_list() |> 
    counterfactual_welfare(x = x_soep) |> 
    arrange(regret)

print(welfare)

load("Data_processed/sampled_chains.RDat")
use_chain <<- combine_chains(chains, 2000, 5000)
plan(multisession, workers = 8)
p = thompson_probabilities(x_soep, use_chain, 200)
print(p)
