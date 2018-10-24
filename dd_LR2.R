dd_LR2 = function(             
   brts,
   initparsoptSR,
   initparsoptCR,
   missnumspec,   
   outputfilename = NULL,
   seed = 42,
   endmc = 1000,
   alpha = 0.05,
   plotit = TRUE,
   res = 10 * (1 + length(brts) + missnumspec),
   ddmodel = 1,
   cond = 1,
   btorph = 1,
   soc = 2,
   tol = c(1E-3,1E-4,1E-6),
   maxiter = 2000,
   changeloglikifnoconv = FALSE,
   optimmethod = 'subplex',
   methode = 'analytical'   
   )
{
  if(!is.null(seed))
  {
     set.seed(roundn(seed))
  }
  if(cond > 1)
  {
    cat("Conditioning on number of tips is not possible.\n")
    return(NULL)
  }
  age = max(brts)
  cat("\nEstimating parameters under the constant-rate model ...\n")  
  outCRO = dd_ML(brts = brts,initparsopt = initparsoptCR,idparsopt = 1:2,idparsfix = 3,parsfix = Inf,res = res,ddmodel = ddmodel,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv,methode = methode,optimmethod = optimmethod)
  cat("\nEstimating parameters under the diversity-dependent model ...\n")  
  outSRO = dd_SR_ML(brts = brts,initparsopt = initparsoptSR,idparsopt = c(1:3,6:7),res = res,ddmodel = ddmodel,missnumspec = missnumspec,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv,methode = methode,optimmethod = optimmethod)
  LRO = outSRO$loglik - outCRO$loglik
  out = cbind(NA,NA,outCRO,outSRO,NA,NA,NA,NA,NA,NA,NA,LRO)
  out = out[,-c(5,7,12,13,17)]
  newnames = c("model","mc","lambda_CR","mu_CR","LL_CR","conv_CR","lambda_SR1","mu_SR1","K1_SR1","K2_SR1","tshift_SR1","LL_SR1","conv_SR1","lambda_SR2","mu_SR2","K1_SR2","K2_SR2","tshift_SR2","LL_SR2","conv_SR2","LR")
  names(out) = newnames
  if(!is.null(outputfilename))
  {
      save(seed,brts,out,file = outputfilename)
  }
  parsCR <<- as.numeric(outCRO[1:2])
  parsSR <<- as.numeric(outSRO[1:7])
  treeCR = list()
  treeSR = list()
  cat('\nSimulating trees under CR and SR models ...\n')
  for(mc in 1:endmc)
  {
     treeCR[[mc]] = dd_sim(pars = c(parsCR,Inf),age = age,ddmodel = ddmodel)
     treeSR[[mc]] = dd_SR_sim(pars = parsSR,age = age,ddmodel = ddmodel)
  }
  if(!is.null(outputfilename))
  {
      save(seed,brts,out,treeCR,treeSR,file = outputfilename)
  }
  cat('\nPerforming bootstrap to determine critical LR ...\n')  
  for(mc in 1:endmc)
  {
     cat('\nAnalyzing simulation:',mc,'\n')
     brtsCR = branching.times(treeCR[[mc]][[1]])
     K1min = max(parsSR[3],10 + length(which(brtsCR > parsSR[7])))
     K2min = max(parsSR[6],10 + length(brtsCR))
     outCR = dd_ML(brtsCR,initparsopt = parsCR,idparsopt = 1:2,idparsfix = 3,parsfix = Inf,res = res,ddmodel = ddmodel,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv,methode = methode,optimmethod = optimmethod)
     outSR1 = dd_SR_ML(brtsCR,initparsopt = c(parsSR[c(1:2)],K1min,K2min,parsSR[7]),idparsopt = c(1:3,6:7),res = res,ddmodel = ddmodel,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv,methode = methode,optimmethod = optimmethod)
     outSR2 = dd_SR_ML(brtsCR,initparsopt = 0.2 * c(rep(1,4),0) + c(parsSR[c(1:2)],K1min,K2min,parsSR[7]),res = res,ddmodel = ddmodel,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv,methode = methode,optimmethod = optimmethod)
     if(outSR1$conv == -1 & outSR2$conv == -1)
     {
        maxLLSR = outCR$loglik
     } else if(outSR1$conv == -1 & outSR2$conv != -1)
     {
        maxLLSR = outSR2$loglik
     } else if(outSR1$conv != -1 & outSR2$conv == -1)
     {
        maxLLSR = outSR1$loglik
     } else {
        maxLLSR = max(outSR1$loglik,outSR2$loglik)
     }   
     LR = pmax(0,maxLLSR - outCR$loglik)
     outff = cbind(1,mc,outCR,outSR1,outSR2,LR)
     outff = outff[,-c(5,7,12,13,17,22,23,27)]     
     names(outff) = newnames
     out = rbind(out,outff)
     if(!is.null(outputfilename))
     {
        save(seed,brts,out,treeCR,treeSR,file = outputfilename)
     }
  }
  opt = rep(0,endmc)
  cat('\nPerforming bootstrap to determine power ...\n')
  for(mc in 1:endmc)
  {
     cat('\nAnalyzing simulation:',mc,'\n')
     brtsSR = sort(branching.times(treeSR[[mc]][[1]]),decreasing = T)
     K1min = max(parsSR[3],10 + length(which(brtsSR > parsSR[7])))
     K2min = max(parsSR[6],10 + length(brtsSR))
     outCR = dd_ML(brtsSR,initparsopt = parsCR,idparsopt = 1:2, idparsfix = 3,parsfix = Inf,res = res,ddmodel = ddmodel,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv,methode = methode,optimmethod = optimmethod)
     outSR1 = dd_SR_ML(brtsSR,initparsopt = c(parsSR[c(1:2)],K1min,K2min,parsSR[7]),idparsopt = c(1:3,6:7),res = res,ddmodel = ddmodel,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv,methode = methode,optimmethod = optimmethod)
     outSR2 = dd_SR_ML(brtsSR,initparsopt = 0.1 * c(rep(1,4),0) + c(parsSR[c(1:2)],K1min,K2min,parsSR[7]),idparsopt = c(1:3,6:7),res = res,ddmodel = ddmodel,missnumspec = 0,cond = cond,btorph = btorph,soc = soc,tol = tol,maxiter = maxiter,changeloglikifnoconv = changeloglikifnoconv,methode = methode,optimmethod = optimmethod)
     if(outSR1$conv == -1 & outSR2$conv == -1)
     {
        maxLLSR = outCR$loglik
        opt[mc] = 1
     } else if(outSR1$conv != -1 & outSR2$conv == -1)
     {
        maxLLSR = outSR1$loglik
        opt[mc] = 2
     } else if(outSR1$conv == -1 & outSR2$conv != -1)
     {
        maxLLSR = outSR2$loglik
        opt[mc] = 3
     } else {
        maxLLSR = max(outSR1$loglik,outSR2$loglik)
        opt[mc] = 1 + min(which(c(outSR1$loglik,outSR2$loglik) == maxLLSR))
     }
     LR = pmax(0,maxLLSR - outCR$loglik)
     outff = cbind(1,mc,outCR,outSR1,outSR2,LR)
     outff = outff[,-c(5,7,12,13,17,22,23,27)]
     names(outff) = newnames
     out = rbind(out,outff)
     if(!is.null(outputfilename))
     {
        save(seed,brts,out,opt,treeCR,treeSR,file = outputfilename)
     }
  }
  inverse_quantile = function(samples,x)
  {
     samplessort = sort(samples)
     pup = which(samplessort > x)
     if(length(pup) > 0)
     {
        if(length(pup) < length(samplessort))
        {
           pup = min(pup)
           invquant = (pup + (x - samplessort[pup])/(samplessort[pup - 1] - samplessort[pup]))/length(samples)
        } else {
           invquant = 0
        }
     } else {
        invquant = 1
     }
     return(invquant)
  }
  funpvalue = function(samples,x)
  {
     samplessort = sort(samples)
     pup = which(samplessort > x)
     pvalue = (length(pup) + 1)/ (length(samples) + 1)
     return(pvalue)
  }
  funpoweroftest = function(samples,x)
  {
     samplessort = sort(samples)
     pup = which(samplessort > x)
     poweroftest = length(pup)/(length(samples) + 1)
     return(poweroftest)     
  }
  #pvalue = 1 - inverse_quantile(out$LR[2:(endmc + 1)],out$LR[1])
  pvalue = funpvalue(out$LR[2:(endmc + 1)],out$LR[1])
  LRalpha = as.numeric(quantile(out$LR[2:(endmc + 1)],1 - alpha,type = 4))
  #poweroftest = 1 - inverse_quantile(out$LR[(endmc + 2):(2 * endmc + 1)],LRalpha)
  poweroftest = funpoweroftest(out$LR[(endmc + 2):(2 * endmc + 1)],LRalpha)
  if(plotit == TRUE)
  {
      try(dev.off())
      try(dev.off())
      pdffilename = paste(getwd(),'/LR.pdf',sep = '')
      pdf(pdffilename,paper = "a4r", width = 29.7, height = 21)
      al = 0.03
      alw = 2
      alw2 = 1.7
      aa = 45
      par(mfrow = c(2,2),cex = 1, mar = c(5, 4, 3, 1) + 0.1)
      hist(out$LR[2:(1 + endmc)],main = 'Distribution of LLR under CR',xlab = 'LLR', ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(0,max(out$LR[1:(endmc + 1)])))
      arrows(out$LR[1],-1E+120, x1 = out$LR[1],y1 = 0, length = al, angle = aa,lwd = alw, col = 'black')
      arrows(LRalpha,-1E+120, x1 = LRalpha,y1 = 0, length = al, angle = aa,lwd = alw, col = 'blue')
      box()
      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", col = "white", col.lab = 'white', col.axis = 'white')
      hist(out$LR[(endmc + 2):(1 + 2 * endmc)], main = 'Distribution of LLR under SR',xlab = 'LLR', ylab = 'Frequency', col = 'red',probability = T,nclass = 30)
      box()
      arrows(out$LR[1],-1E+120, x1 = out$LR[1],y1 = 0, length = al, angle = aa,lwd = alw, col = 'black')
      arrows(LRalpha,-1E+120, x1 = LRalpha,y1 = 0, length = al, angle = aa,lwd = alw, col = 'blue')

      par(mfrow = c(2,3),cex = 1, mar = c(5, 4, 3, 1) + 0.1)
      lambda = out$lambda_CR[(endmc + 2):(2 * endmc + 1)] * (opt == 1) + out$lambda_SR1[(endmc + 2):(2 * endmc + 1)] * (opt == 2) + out$lambda_SR2[(endmc + 2):(2 * endmc + 1)] * (opt == 3)
      mu = out$mu_CR[(endmc + 2):(2 * endmc + 1)] * (opt == 1) + out$mu_SR1[(endmc + 2):(2 * endmc + 1)] * (opt == 2) + out$mu_SR2[(endmc + 2):(2 * endmc + 1)] * (opt == 3)
      K1 = 1E+120 * (opt == 1) + pmin(1E+120,out$K1_SR1[(endmc + 2):(2 * endmc + 1)]) * (opt == 2) + pmin(1E+120,out$K1_SR2[(endmc + 2):(2 * endmc + 1)]) * (opt == 3)
      K2 = 1E+120 * (opt == 1) + pmin(1E+120,out$K2_SR1[(endmc + 2):(2 * endmc + 1)]) * (opt == 2) + pmin(1E+120,out$K2_SR2[(endmc + 2):(2 * endmc + 1)]) * (opt == 3)
      tshift = 0 * (opt == 1) + out$tshift_SR1[(endmc + 2):(2 * endmc + 1)] * (opt == 2) + out$tshift_SR2[(endmc + 2):(2 * endmc + 1)] * (opt == 3)
      hist(lambda,main = NULL, xlab = expression(lambda), ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(0,max(lambda)))
      arrows(out$lambda_SR1[1],-1E+120, x1 = out$lambda_SR1[1],y1 = 0, length = al, angle = aa,lwd = alw2, col = 'black')
      box()
      hist(mu,main = NULL, xlab = expression(mu), ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(0,max(mu)))
      arrows(out$mu_SR1[1],-1E+120, x1 = out$mu_SR1[1],y1 = 0, length = al, angle = aa,lwd = alw2, col = 'black')
      box()
      hist(K1,main = NULL, xlab = 'K1', ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(min(K1),max(K1)))
      arrows(out$K1_SR1[1],-1E+120, x1 = out$K1_SR1[1],y1 = 0, length = al, angle = aa,lwd = alw2, col = 'black')
      box()
      hist(K2,main = NULL, xlab = 'K2', ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(min(K2),max(K2)))
      arrows(out$K2_SR1[1],-1E+120, x1 = out$K2_SR1[1],y1 = 0, length = al, angle = aa,lwd = alw2, col = 'black')
      box()
      hist(tshift,main = NULL, xlab = 't_shift', ylab = 'Frequency', col = 'red',probability = T,nclass = 30, xlim = c(min(tshift),max(tshift)))
      arrows(out$tshift_SR1[1],-1E+120, x1 = out$tshift_SR1[1],y1 = 0, length = al, angle = aa,lwd = alw2, col = 'black')
      box()
      
      try(dev.off())
      try(dev.off())
      os = .Platform$OS.type
      if(os == "windows")
      {
          shell.exec(pdffilename)
      }
      if(os == "unix")
      {
          system(paste("open",pdffilename,sep = " "))
      }
  }
  if(!is.null(outputfilename))
  {
      save(seed,brts,out,opt,treeCR,treeSR,pvalue,LRalpha,poweroftest,file = outputfilename)
  }
  return(list(treeCR = treeCR,treeSR = treeSR,out = out,pvalue = pvalue,LRalpha = LRalpha,poweroftest = poweroftest))
}