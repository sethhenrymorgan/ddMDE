# Minimum Detectable Effect Calculation by Simulation & Optimization

# see here: https://egap.org/resource/script-power-analysis-simulations-in-r/
# and here for Stata version: https://github.com/sethhenrymorgan/ddpower/blob/main/ddpowersimu.ado
# also see Stata FAQ: https://www.stata.com/support/faqs/statistics/power-by-simulation/
# on measuring ICC in Stata: https://blogs.worldbank.org/en/impactevaluations/tools-of-the-trade-intra-cluster-correlations
# loneway Stata command for measuring ICC: https://www.stata.com/manuals/rloneway.pdf
# and Stata List on generating clustered data: https://www.statalist.org/forums/forum/general-stata-discussion/general/1460617-how-to-simulate-clustered-data-with-a-specific-intra-class-correlation
# On clustering: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1466680/
# and source of some cluster code: https://stats.stackexchange.com/questions/263451/create-synthetic-data-with-a-given-intraclass-correlation-coefficient-icc
# estimating linear models with clustered standard errors (lm_robust): https://www.r-bloggers.com/2021/05/clustered-standard-errors-with-r/
# optimization: https://search.r-project.org/R/refmans/stats/html/optimize.html
# another nice power by simulation in R tutorial: https://cameronraymond.me/blog/power-simulations-in-r/

library(dplyr)
library(tictoc)
library(estimatr)
library(stats)

# Simulation Parameters
alpha <- 0.05 # standard significance level
sims <- seq(1,100, by=1) # number of simulations
power_target <- 0.8 # targeted statistical power (probably stay at 0.8)

# Sample Parameters (set for theory sample, contraception dataset)
N <- 3313 # set sample size
C <- 130 # number of clusters
n_per_c <- N/C
#treatprob <- .48864307 # fraction of sample in treatment
#afterprob <- 0.76755492 # fraction of sample post-treatment time
treatprob_byc <- 0.3 # fraction of clusters in treatment
afterprob_byc <- .82758621 # fraction of clusters in post-treatment time (2016+2022?)

# Variable Parameters (Currently Modern Contraception Use)
mu <- .2659221  # set outcome mean
#ICC <- 0.18598 # intra-cluster correlation (turns out don't need this if we have the below parameters)
sd_alpha <- .1238547 # variation between clusters
sd_epsilon <- .4244322 # variation within clusters
interval <- c(0.1,0.5)   # a reasonable interval over which to search for the MDE
  
# function to generate clustered data
gen_data_c <- function(c, C, N){
  n_per_c <- N/C
  cluster_c <- data.frame(C=c,
                          alpha_c = rnorm(n = 1, mean = 0, sd = sd_alpha),
                          epsilon_ic = rnorm(n = n_per_c, mean = 0, sd = sd_epsilon),
                          Z.sim = rbinom(n=1, size=1, prob = treatprob_byc),
                          Time = rbinom(n=1, size=1, prob = afterprob_byc)) %>%
    mutate(Y0 = mu + alpha_c + epsilon_ic)
  return(cluster_c)
}

# function to generate vector of significant trials
gen_sig <- function(c, C, N, tau){
  df <- lapply(X = 1:C, FUN = gen_data_c, N=N, C=C) %>% 
    bind_rows() %>%
    mutate(C = as.factor(C)) %>%
    mutate(Y1 = Y0 + tau) %>%
    mutate(Y.sim = Y1*(Z.sim*Time) + Y0*(1-(Z.sim*Time)))
  
  fit.sim <- lm_robust(Y.sim ~ Z.sim + Time + Z.sim:Time, clusters = C, data = df) # Do analysis (Simple regression)
  p.value <- summary(fit.sim)$coefficients[4,4]  # Extract p-values
  significant.experiments <- (p.value <= alpha) # Determine significance according to p <= 0.05
  return(significant.experiments)
}

# Function to minimize: calculates power and distance between power and power target
power_diff <- function(tau, sims, N, C) {
  significant.experiments <- lapply(X = sims, FUN = gen_sig, N=N, C=C, tau=tau)
  power <- mean(as.integer(significant.experiments))
  diff = (power - power_target)^2
}

# Use optimize to get MDE
tic("optimize")
MDE <- optimize(f = power_diff, interval = interval, tol = 0.01, sims = sims, N=N, C=C)
toc()

