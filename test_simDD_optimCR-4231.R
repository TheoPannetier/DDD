library('DDD')
load("data_simDD_optimCR-4231.RData")

new_res <- NULL
cond = 1
tol = rep(1E-6,3)
methode = "ode45"
optimmethod = "subplex"

for(i in seq_along(trees)){
  
  brts = as.numeric(branching.times(trees[[i]]))
  initparsopt = c(res$init_lambda0[i], res$init_mu0[i], Inf)
  new_res[[i]] <- DDD::bd_ML(brts = brts, 
                             initparsopt = initparsopt[1:2] + 1E-6,
                             cond = cond, 
                             tol = tol, 
                             methode = methode,
                             optimmethod = optimmethod
  )
  
}
