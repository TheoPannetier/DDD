dd_AICb <- function(brts,initparsopt,idparsopt,parsfix,idparsfix,idparsnoshift,simmodel = 'DD_SR',simdatafilename = 'simdata.RData',endmc = 10000,seed = 42,ddmodel = 1,res = 100,cond = 1,btorph = 1,soc = 2,startmc = 1,methode = 'analytical')
{
  set.seed(seed)
  brts <- abs(brts)
  age <- max(brts)
  simdata <- list()
  if(simmodel == 'DD' | simmodel == 'CR')
  {
    simfunction <- dd_sim
    MLfunction <- dd_ML
    LLfunction <- dd_loglik
    numpars = 3
  } else
  {
    simfunction <- dd_SR_sim
    MLfunction <- dd_SR_ML
    LLfunction <- dd_SR_loglik
    numpars = 7
  }
  out <- MLfunction(brts = brts,initparsopt = initparsopt,idparsopt = idparsopt,parsfix = parsfix,idparsfix = idparsfix,idparsnoshift = idparsnoshift,cond = cond,btorph = btorph,res = res,soc = soc,methode = methode)
  MLpars <- as.numeric(out[1:numpars])

  if(startmc == 1)
  {
    for(mc in 1:endmc)
    {
       print(mc)
       flush.console()
       simdata[[mc]] <- simfunction(pars = MLpars,age = age,ddmodel = ddmodel)
       save(brts,out,simdata,file = simdatafilename)
    }
  }
  if(startmc > 1)
  {
    load(simdatafilename)
  }
  MLLo <- as.numeric(out$loglik)
  simresults = list()
  MLLb <- rep(0,endmc)
  loglikbo <- rep(0,endmc)
  for(mc in startmc:endmc)
  {
    simbrts <- simdata[[mc]]$brts
    simresults[[mc]] <- MLfunction(brts = simbrts,initparsopt = MLpars[idparsopt],idparsopt = idparsopt,parsfix = parsfix,idparsfix = idparsfix,idparsnoshift = idparsnoshift,res = res,cond = cond,btorph = btorph,soc = soc,methode = methode)
    MLLb[mc] <- as.numeric(simresults[[mc]]$loglik)
    loglikbo[mc] <- LLfunction(brts = brts,pars1 = as.numeric(simresults[[mc]][1:numpars]),pars2 = c(res,ddmodel,cond,btorph,0,soc),missnumspec = 0,methode = methode)
    save(brts,out,simdata,simresults,MLLo,MLLb,loglikbo,file = simdatafilename)
  }
  AICb <- rep(0,2)
  AICb[1] <- AICb1(MLLo,MLLb,loglikbo)
  AICb[2] <- AICb2(MLLo,loglikbo)
  save(brts,out,simdata,simresults,MLLo,MLLb,loglikbo,AICb,file = simdatafilename)
return(list(brts,out,simdata,simresults,MLLo,MLLb,loglikbo,AICb))
}
