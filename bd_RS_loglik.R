expfun <- function(la,mu,t,t0)
{
  res <- exp(-(la - mu) * (t - t0))
  return(res)
}

logH <- function(la,mu,t,tp = 0)
{
  ef <- expfun(la,mu,t = tp,t0 = t)
  if(la > mu)
  {
     res <- 2 * log(la - mu) - (la - mu) * (tp - t) - 2 * log(la - mu * ef)
  } else if(la < mu)
  {
     res <- 2 * log(mu - la) - (la - mu) * (tp - t) - 2 * log(mu * ef - la)
  } else
  {
     res <- -2 * log(1 + mu * (tp - t))
  }
  return(res)
}

logLL <- function(pars1,brts,tshift)
{
  brts <- -abs(brts)
  loglik <- (length(unique(brts)) - 1) * log(pars1[1])
  for(i in 1:length(brts))
  {
     loglik <- loglik + logH(la = pars1[1],mu = pars1[2],t = brts[i])
  }
  if(length(tshift) > 0)
  {
     tshift <- -abs(tshift)
     for(i in 1:length(tshift))
     {
       loglik <- loglik - logH(la = pars1[1],mu = pars1[2],t = tshift[i])
     }
  }
  return(loglik)
}

bd_RS_loglik <- function(pars1,pars2 = NULL,brts,tshift)
{
   loglik <- 0
   for(i in 1:length(brts))
   {
     loglik <- loglik + logLL(pars1 = pars1[[i]],brts = brts[[i]],tshift = tshift[[i]])
   }
   return(loglik)
}

bd_RS_loglik_BAMM <- function(BAMMtable)
{
  BAMMtable <- as.matrix(BAMMtable)
  loglik <- 0
  tp <- max(BAMMtable[,3])
  for(i in 1:dim(BAMMtable)[1])
  {
    loglik <- loglik + 
      logH(la = BAMMtable[i,5],mu = BAMMtable[i,6],t = BAMMtable[i,2],tp = tp) -
      logH(la = BAMMtable[i,5],mu = BAMMtable[i,6],t = BAMMtable[i,3],tp = tp)
  }
  sortBAMMtable <- pracma::sortrows(BAMMtable,2)
  brpts <- which(duplicated(sortBAMMtable[,2]))
  loglik <- loglik + sum(log(sortBAMMtable[brpts[2:length(brpts)],5]))
  return(as.numeric(loglik))
}
