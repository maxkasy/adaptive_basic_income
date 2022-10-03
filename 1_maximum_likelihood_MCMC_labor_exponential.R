library(furrr) # for parallel chains

# function to prep evaluation of likelihood, memoizing all redundant calculations
# TBD: clearner alternative to using global environment
prep_likelihood = function(dat){
    # for numerical integration over the normal distribution
    gridpoints <<- 20
    epsilon <<- qnorm(seq(.5/gridpoints, by = 1/gridpoints, length.out = gridpoints))
    
    n_obs <<- nrow(dat)
    x <<- dat  |>  
        select(starts_with("x")) |> 
        as.matrix()
    y <<- dat$y
    y0 <<- dat$y0
    w <<- dat$w
    
    part <<- (y > 0) # indicator for positive earnings
    cutoff <<- ifelse(y0 > 0, y0/(1 - w), 0)
    # top <<- (y > cutoff) # indicator for being in no transfer regime
    log_y  <<-  log(y)
}


# this relies on prep_likelihood being executed previously
log_likelihood = function(param) {
    if (param$beta <= 0 | param$sigma <=0 | param$tau <= 0 | param$rho <=0) 
            {LL = -Inf} # imposing constraint for sensible parameter values
        else {
        LL = 0 # log likelihood (we will keep adding)
        
        mu_alpha = x %*% param$gamma_alpha
        mu_eta = x %*% param$gamma_eta
        
        exp_epsilon_sigma = exp(param$sigma * epsilon) # for numerical integration over the log-normal distribution
        
        # cutoff for choosing the top bracket:
        alpha_cutoff = ifelse(y0 > 0,
                    log(param$beta * y0) - log(exp(param$beta)-exp(param$beta * w)), 
                    -Inf)
        for (i in 1:n_obs) { # loop over all observations
            if (part[i]) { # obs with positive income
                alpha_hat1 = log_y[i] - param$beta * w[i]
                Delta_hat1 = y[i] * (1- exp(-param$beta * w[i])) / param$beta
                
                alpha_hat2 = log_y[i] - param$beta
                Delta_hat2= y[i] * (1- exp(-param$beta))/ param$beta - y0[i]
           
                ll = log(
                        pnorm(alpha_cutoff[i] - alpha_hat1, sd = param$rho) * 
                            dnorm(alpha_hat1, mu_alpha[i], param$sigma) +
                        pnorm(- alpha_cutoff[i] + alpha_hat2, sd = param$rho) * 
                            dnorm(alpha_hat2, mu_alpha[i], param$sigma) +
                            .0001 # to avoid NaNs for small sigma
                       ) +
                    pnorm(max(Delta_hat1, Delta_hat2) * param$tau  + mu_eta[i], log = T)

            } else { # obs without positive income
                L = 0 # likelihood
                Delta_hat_mu = exp(mu_alpha[i]) * # Delta hat, up to multiplication by exp(epsilon)
                    max(
                        (exp(param$beta * w[i]) - 1) / param$beta,
                        (exp(param$beta) - 1) / param$beta - y0[i]
                    )
                for (q in 1:gridpoints) { # loop over quantiles of alpha    
                    L = L + pnorm(- exp_epsilon_sigma[q] * Delta_hat_mu * param$tau  - mu_eta[i])
                }
                
                ll = log(L/gridpoints + .25/gridpoints)
            }
            LL = LL + ll
            if (ll == -Inf | is.nan(ll)) browser()
        }
        }
    
    LL
}

# wrapper functions for optimizatiton
par_vector_to_list =  function(par_vector){
    dim_x = (length(par_vector) - 3)/2
    
    list(
        beta = par_vector[1],
        sigma = par_vector[2],
        tau = par_vector[3],
        rho = par_vector[4],
        gamma_alpha = par_vector[5:(5+dim_x-1)],
        gamma_eta =  par_vector[(5+dim_x):(5+2*dim_x-1)]
    )
}


log_likelihood_vectorarg = function(par_vector) {
    ll = par_vector |>
        par_vector_to_list() |>
        log_likelihood()
    cat(par_vector, "\n",
        ll, "\n\n")
    ll
}



######## Metropolis algorithm ################ ----
# based on: https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/

log_posterior = function(par_vector){
    log_likelihood_vectorarg(par_vector) + # Add log prior for regularization
        + 2 * (log(max(par_vector[3],0)) - 
                   log(max(par_vector[2],0))) + #- # uniform pseudo-prior for log sigma
        dnorm(par_vector[4], sd = .1, log = T) # informative prior for optimization error
         # dnorm(par_vector[3], sd = 5000, log = T) +
         # dnorm(par_vector[4], sd = 5000, log = T) +
         # dnorm(par_vector[5 + dim_x], sd = 5000, log = T) +
         # dnorm(par_vector[5 + dim_x + 1], sd = 5000, log = T)
}

proposalfunction = function(par_vector){
    MASS::mvrnorm(n = 1, 
                  mu = par_vector,
                  Sigma = (5/(dim_par)) * Sigma) 
    # Optimal proposal distribution for multivariate normal, using MAP hessian
}

run_metropolis = function(start_value, iterations){
    dim_par = length(start_value)
    chain = array(dim = c(iterations+1,dim_par))

    chain[1,] = start_value
    lp_old = log_posterior(start_value)
    for (i in 1:iterations){
        proposal = proposalfunction(chain[i,])
        # cat(i, "\n", 
            # proposal, "\n")
        lp_new = log_posterior(proposal)
        
        probab = exp(lp_new - lp_old)
        # cat(paste(lp_old, lp_new, probab),
            # "\n\n")
        
        if (runif(1) < probab){
            chain[i+1,] = proposal
            lp_old = lp_new
        }else{
            chain[i+1,] = chain[i,]
        }
    }
    return(chain)
}


# Implementing parallel chains
# list of start values determines number of chains
parallel_metropolis = function(start_values, iterations){
    no_chains = length(start_values)
    plan(multisession, workers = no_chains)
    
    future_map(start_values, 
               ~run_metropolis(.x, iterations),
               seed = T)
}

combine_chains = function(chains, burn_in, iterations){
    do.call(rbind,
         map(chains, ~.x[burn_in:iterations, ]) 
    )
}

# master function for full analysis, including tuning with the MLE
metropolis_master = function(data, dim_x, 
                             start_value,
                             # MCMC parameters
                             iterations = 5000,
                             burn_in = 2000,
                             no_chains = 8){

    # function to prep evaluation of likelihood, memoizing all redundant calculations
    prep_likelihood(data)
    # demean covariates for greater numerical stability
    for (i in seq(2,dim_x,length = max(0,dim_x-1))) {
        x[,i] <<- x[,i] - mean(x[,i])
    }
    
    # Running MAP (penalized MLE) -----
    MAP <<-
        optim(start_value,
              log_posterior,
              method = "Nelder-Mead",# "L-BFGS-B",
              control=list(fnscale=-1)#,
                           # maxit = 100)#, maximum iterations for debugging. drop later.
              # hessian = T,
              # lower = c(c(0.001,0.001,0.00005,0.001), rep(-Inf,  dim_par - 4))
              )
    
    
    starting_value_MCMC = MAP$par
    # hessian = - MAP$hessian
    hessian = - numDeriv::hessian(log_posterior, starting_value_MCMC)
    
    # small / negative eigenvalues are indicative of weak identification    
    cat("Eigenvalues of the Hessian:\n")
    eigen_hessian =  eigen(hessian)
    
    eigen_hessian$values |> 
        print()

    # Sigma <<- solve(hessian)
    # Find information matrix as nearest positive definite matrix to inverse hessian of likelihood
    Sigma <<-
        eigen_hessian$vectors %*%
        diag(1/pmax(eigen_hessian$values, 10^(-2))) %*%
        t(eigen_hessian$vectors)

    # Printing MAP results for intermediate information
    tibble(
        MAP = MAP$par,
        truth = unlist(true_par),
        MAP_se = sqrt(diag(Sigma))
    ) |> 
        print()
    
    # Running MCMC ----
    # define start values as random draws from asymptotic normal distribution of MAP
    start_values = map(1:no_chains,
                       ~ MASS::mvrnorm(n = 1, 
                                       mu = starting_value_MCMC,
                                       Sigma = Sigma) )

    # Catch for starting values with -Inf log prior
    for (i in 1:no_chains) {
        if (log_posterior(start_values[[i]]) == -Inf) browser()
    }
    
    chains <<- parallel_metropolis(start_values, iterations)
    
    # Store MCMC results for later usage
    save(chains,
         file = "Data_processed/sampled_chains.RDat")
    
    # drop burn in period and merge chains
    use_chain <<- combine_chains(chains, burn_in, iterations)
 
    # summarize estimates
    estimates <<- 
        tibble(
            Parameter = names(start_value),
            truth = unlist(true_par),
            MAP = MAP$par,
            MAP_se = sqrt(diag(Sigma)),
            Bayes = colMeans(use_chain),
            Bayes_se = apply(use_chain, 2, sd)
        ) 
    
    estimates <<- estimates |> 
       bind_rows(tibble(
            Parameter = c("log posterior", "log likelihood"),
            truth = c(log_posterior(estimates$truth),log_likelihood_vectorarg(estimates$truth)),
            MAP = c(log_posterior(estimates$MAP),log_likelihood_vectorarg(estimates$MAP)),
            MAP_se = c(NA, NA),
            Bayes = c(log_posterior(estimates$Bayes),log_likelihood_vectorarg(estimates$Bayes)),
            Bayes_se = c(NA, NA)
        ))
    
    print(estimates)
}



# MCMC diagnostics
MCMC_diagnostics = function() {
    # acceptance rate of MCMC proposals
    cat("Acceptance rate of proposals for Metropolis-Hastings:\n")
    print(1-mean(duplicated(use_chain)))
    
    library(coda)
    mcmc_chain = mcmc(use_chain)

    # png("Data_processed/traceplot.png")
    plot(mcmc_chain)
    # traceplot(mcmc_chain)
    # densplot(mcmc_chain)
    # dev.off()
    
    
    
    library(BayesianTools)
    correlationPlot(data.frame(use_chain))
    
    combinedchains = mcmc.list(map(chains, mcmc))
    # plot(combinedchains)
    gelman.diag(combinedchains)
    gelman.plot(combinedchains)
}
