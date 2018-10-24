pars1 = c(0.8,0.2,140)
pars2 = c(100,1,1,0,1,2)
brts = 1:30
missnumspec = 0
methode = 'lsoda'
system.time(r2 <- dd_loglik(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,rhs_func_name = 'dd_loglik_rhs',methode = methode))
system.time(r3 <- dd_loglik(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,rhs_func_name = 'dd_loglik_rhs_FORTRAN',methode = methode))
r3 - r2

#probs <- dd_FORTRAN(initprobs = initprobs, t0 = t0, t = t, totmat = totmat, methode = 'lsoda')

dyn.unload(paste("d:/data/ms/DDD/dd_loglik_rhs_FORTRAN", .Platform$dynlib.ext, sep = ""))
system("R CMD SHLIB d:/data/ms/DDD/dd_loglik_rhs_FORTRAN.f")
dyn.load(paste("d:/data/ms/DDD/dd_loglik_rhs_FORTRAN", .Platform$dynlib.ext, sep = ""))
system.time(r3 <- dd_loglik_new(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,rhs_func_name = 'dd_loglik_rhs_FORTRAN',methode = methode))
