# Set directory ; load libraries
setwd("d:/data/ms/DDD/")
library(ape) ; library(DDD)

# Load data
load("sample_trees_test.RData")

# decalre dd_ML arguments
initpars <-  c(60, 0.8, 0.4, 40)
methode = 'lsodes'
#optimmethod = 'subplex'
optimmethod = 'simplex'

outerror = data.frame(lambda = -1,mu = -1,K = -1, loglik = -1, df = -1, conv = -1)
results = NULL

for( i in seq_along(trees)){
  mc <- names(trees)[i]
  print(paste("Optimizing on tree",mc))
  
  # fit DD
  results_mc = try(dd_ML(as.numeric(branching.times(trees[[i]])), initparsopt = initpars[2:4]+1E-6, cond = 1, tol = rep(1E-6,3), 
                         methode = methode, optimmethod = optimmethod))
  if(!is.data.frame(results_mc)) { results_mc = outerror }
  
  # Assemble results
  results_mc <- cbind(results_mc, "mc" = mc)
  results <- rbind(results,results_mc)
}
#save(results, file = ".RData")

