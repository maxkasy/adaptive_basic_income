library(tidyverse)

source("1_maximum_likelihood_MCMC_labor_exponential.R")
source("2_labor_supply_simulate_functions.R")

# parameters which work nicely
true_par = list(
    beta = 0.2,
    sigma = .5,
    tau = 1/1000,
    rho = .1,
    gamma_alpha = c(7,1,-1),
    gamma_eta =  c(-.5, -1, 1)
)

dim_par = length(unlist(true_par))
dim_x = length(true_par$gamma_alpha)
start_value = unlist(true_par) * runif(dim_par, min= .5, max = 1.5) # starting value for MLE

x = sample_x(1000,dim_x)
sim_random = simulate_sample(x = x, par = true_par, n_per_policy = 200)

sim_random |>
    filter(y <5000) |>
    counterfactual_plots(filename = "Figures/test_plots.png")


metropolis_master(sim_random, dim_x, start_value)

MCMC_diagnostics()


welfare =
    estimates |> 
    select(truth, MAP, Bayes) |> 
    slice(1:dim_par) |> 
    map(~ counterfactual_welfare(par_vector_to_list(.x),
                                x = x))

print(welfare)

regret_comparison = 
    map(welfare, "regret") |> 
    as_tibble() |> 
    mutate(lab = paste0(welfare[[1]]$y0, ", ", welfare[[1]]$w))
     

regret_comparison |> 
    ggplot(aes(x=truth, y=Bayes, label = lab)) +
    geom_abline(slope=1,intercept=0) +
    geom_text() +
    # geom_text(aes(y=MAP), color = "red") +
    theme_minimal()
    

thompson_probabilities(x, use_chain, 200)


