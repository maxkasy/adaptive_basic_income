library(tidyverse)
library(patchwork) # for combining plots

# set.seed(2046) # for consistent simulations

# Simulation parameters ----
# Set of policies for which we simulate data
sim_policies =
    tibble(y0 = c(0, rep(seq(500, 1500, by = 500), 3)),
           w = c(1, rep(seq(.2, .8, by = .3), each = 3)))

# Covariate distribution from SOEP
x_soep = read_csv("Data_processed/covariates_SOEP.csv") |>
    as.matrix()
n_soep = nrow(x_soep)


# Normative parameters for counterfactual welfare
welfare_exponent = .4 # smaller value means more egalitarian preferences
break_even_u = 1800 # value of u where marginal individual euro is worth same as for funds
value_of_funds = welfare_exponent * break_even_u^(welfare_exponent-1) # value of a marginal euro, in units of u^welfare_exponent


# Simulation functions ----
# functions to simulate individual heterogeneity, based on these parameter values
sample_x = function(n, dim_x){
   x = matrix(c(rep(1, n), runif(n * (dim_x - 1))), ncol = dim_x)
   # demean covariates for greater numerical stability
   for (i in seq(2,dim_x,length = max(0,dim_x-1))) {
       x[,i] = x[,i] - mean(x[,i])
   }
   x
}
sample_alpha = function(x, par, n)
    rnorm(n,
          mean = x %*% par$gamma_alpha,
          sd = par$sigma)
sample_eta = function(x, alpha, par, n)
    rnorm(
        n,
        mean = - (x %*% par$gamma_eta) / par$tau,
        sd =  1 / par$tau
    )

# simulating outcomes as a function of policy parameters
simulate_data = function(y0, w, x, par) {
    # parameters to use for simulation
    if (!is.matrix(x)) {x = matrix(x)}
    n = nrow(x)
    dim_x = ncol(x)
    colnames(x) = paste0("x", 1:dim_x)
    
    alpha = sample_alpha(x, par, n)
    eta = sample_eta(x, alpha, par, n)
    
    y1 = exp(alpha + par$beta * w)
    y2 = exp(alpha + par$beta)
    Delta1 =  exp(alpha) * (exp(par$beta * w) - 1) / par$beta 
    Delta2 =  exp(alpha) * (exp(par$beta) - 1) / par$beta - y0
    
    y = rep(0, n)
    u = rep(0, n)
    for (i in 1:n) {
        if (max(Delta1[i], Delta2[i]) <= eta[i]) {
            u[i] = y0
        } else if ((Delta2[i] + y0 - Delta1[i]) * exp(rnorm(1, sd = par$rho)) < y0) { 
            # note the introduction of optimization error in the choice between schedules
            y[i] = y1[i]
            u[i] = y0 + Delta1[i] - eta[i]
        } else {
            y[i] = y2[i]
            u[i] = y0 + Delta2[i] - eta[i]
        }
    }
    
    tibble(
        y0 = y0,
        w = w,
        alpha = alpha,
        eta = eta,
        Delta1 = Delta1,
        Delta2 = Delta2,
        y = y,
        u = u,
        part = (y > 0),
        transf = pmax(y0 + (w - 1) * y, 0),
        v = u ^ welfare_exponent / value_of_funds,
        swf = v - transf
    ) |> 
        bind_cols(x |> as_tibble())
}


# Generate simulated data for the different policies
simulate_sample = function(x = x_soep, par = true_par, n_per_policy = 200) {
    map(1:nrow(sim_policies),
        ~ simulate_data(
            sim_policies[[.x, "y0"]],
            sim_policies[[.x, "w"]],
            x = x[sample(nrow(x), n_per_policy),],
            par = par)) |> 
        bind_rows()
}


# Using simulated data to calculate policy counterfactuals
counterfactual_welfare = function(x = x_soep,
                                  par = true_par) {
    map(1:nrow(sim_policies),
        ~ simulate_data(
            sim_policies[[.x, "y0"]],
            sim_policies[[.x, "w"]],
            x = x,
            par = par) |> 
            summarise_all(mean) |> 
            select(y0, w, part, y, u, v, transf, swf)) |> 
        bind_rows() |>
        mutate(regret = max(swf)-swf)
    
}


# Thompson probabilities, based on MCMC draws
thompson_probabilities = function(x,
                                  mcmc_chain,
                                  no_draws) {
    
    sample_optimal =
        mcmc_chain[sample(nrow(mcmc_chain),size=no_draws,),] |> 
        split(rep(1:no_draws, times = dim_par)) |> 
        future_map(~ counterfactual_welfare(par_vector_to_list(.x),
                                     x = x)|> 
                       pull(swf) |> which.max()) 
    
    sim_policies |> 
        mutate(p_optimal =
                   tabulate(unlist(sample_optimal), nbins = nrow(sim_policies)) / no_draws)
}



# Plot simulated data for different policies
counterfactual_plots = function(sim_dat,
                                filename = "Figures/counterfactual_plots.png") {
    p1 = sim_dat |>
        ggplot(aes(x = u)) +
        geom_histogram(bins = 15) +
        facet_grid(w ~ y0) +
        theme_light() +
        labs(title = "Distribution of utility",
             x = "", y = "")
    
    p2 = sim_dat |>
        ggplot(aes(x = y)) +
        geom_histogram(bins = 15) +
        facet_grid(w ~ y0) +
        theme_light() +
        labs(title = "Distribution of earnings",
             x = "", y = "")
    
    p3 = sim_dat |>
        ggplot(aes(x = transf)) +
        geom_histogram(bins = 15) +
        facet_grid(w ~ y0) +
        theme_light() +
        labs(
            title = "Distribution of transfers",
            x = "",
            y = "",
            caption = "Columns correspond to different values of the basic income size.
                   Rows correspond to different values of the net-of-tax rate."
        )
    
    ggsave(
        filename,
        p1 / p2 / p3,
        width = 8,
        height = 12
    )
    
}




