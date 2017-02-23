expfun <- function(la,mu,t,t0)
{
  res <- exp(-(la - mu) * (t - t0))
  return(res)
}

logP0 <- function(la,mu,t0,t,f)
{
  if(la > mu)
  {
    res <- log(mu) + log(1 - f) - log(la - mu * f)
  } else if(la < mu)
  {
    res <- log(mu) + log(f - 1) - log(mu * f - la)
  } else
  {
    res <- log(mu) + log(t - t0) - log(1 + mu * (t - t0))
  }
  return(res)
}

log1minusP0 <- function(la,mu,t,t0,f)
{
  if(la > mu)
  {
    res <- log(la - mu) - log(la - mu * f)
  } else if(la < mu)
  {
    res <- log(mu - la) - log(mu * f - la)
  } else
  {
    res <- -log(1 + mu * (t - t0))
  }
  return(res)
}

logu <- function(la,mu,t,t0,f)
{
  if(la > mu)
  {
    res <- log(la) + log(1 - f) - log(la - mu * f)
  } else if(la < mu)
  {
    res <- log(la) + log(f - 1) - log(mu * f - la)
  } else
  {
    res <- log(mu) + log(t - t0) - log(1 + mu * (t - t0))
  }
  return(res)
}

log1minusu <- function(la,mu,t,t0,f)
{
  if(la > mu)
  {
    res <- log(la - mu) + log(f) - log(la - mu * f)
  } else if(la < mu)
  {
    res <- log(mu - la) + log(f) - log(mu * f - la)    
  } else
  {
    res <- -log(1 + mu * (t - t0))
  }
  return(res)
}

take.matrix.element <- function(mat,rownum = 1,colnum = 1)
{
  return(mat[rownum,colnum])
}

#The data are assumed to be a list where every element
#of the list refers to a single broken fragment after 
#the breaking the tree procedure of Nee et al. 1994
#Each element is a matrix with in the first row the 
#branching times and shift times and in the second
#and third row the corresponding speciation and extinction
#rates. The speciation event that gave rise to the fragment
#is assumed to have the first set of rates.

SR_bd_loglik <- function(dataparslist,cond = 1,soc = 2)
{
  loglik <- 0
  for(i in 1:length(dataparslist))
  {
    times <- c(-abs(dataparslist[[i]][1,]),0)
    la <- dataparslist[[i]][2,]
    mu <- dataparslist[[i]][3,]
    nj <- length(times) - 1
    for(j in 1:nj)
    {
      f <- expfun(la = la[j],mu = mu[j],t = times[j + 1],t0 = times[j])
      loglik <- loglik + log(la[j]) + 
        log1minusP0(la = la[j],mu = mu[j],t = times[j + 1],t0 = times[j],f = f) +
        log1minusu(la = la[j],mu = mu[j],t = times[j + 1],t0 = times[j],f = f)
      if(j < nj)
      {
        f2 <- expfun(la = la[j],mu = mu[j],t = times[j + 2],t0 = times[j + 1])
        loglik <- loglik - 2 * log(1 - exp
          (logu(la = la[j],mu = mu[j],t = times[j + 1],t0 = times[j],f = f) +
          logP0(la = la[j],mu = mu[j],t = times[j + 2],t0 = times[j + 1],f = f2))
        )
      }
    }
  }
  cf <- which.max(unlist(lapply(dataparslist,take.matrix.element, rownum = 1,colnum = 1))) #crown branch
  cla <- dataparslist[[cf]][2,1]
  if(cond == 1)
  {
     loglik <- loglik - soc * log(cla) #remove ancestral lineage contribution              
  }
  return(loglik)
}

SR_bd_loglik2 <- function(t0,t1,ts,t2,la1,la2,mu1,mu2)
{
  f1 <- expfun(la = la1,mu = mu1,t = ts,t0 = t0)
  f2 <- expfun(la = la1,mu = mu1,t = ts,t0 = t1)
  f3 <- expfun(la = la1,mu = mu1,t = t2,t0 = ts)
  f4 <- expfun(la = la2,mu = mu2,t = t2,t0 = ts)
  loglik <- log(la1) + 
            log1minusP0(la = la1,mu = mu1,t = ts,t0 = t0,f = f1) + 
            log1minusu(la = la1,mu = mu1,t = ts,t0 = t0,f = f1) + 
            log1minusP0(la = la1,mu = mu1,t = ts,t0 = t1,f = f2) + 
            log1minusu(la = la1,mu = mu1,t = ts,t0 = t1,f = f2) + 
            log1minusP0(la = la1,mu = mu1,t = t2,t0 = ts,f = f3) + 
            log1minusu(la = la1,mu = mu1,t = t2,t0 = ts,f = f3) + 
            log1minusP0(la = la2,mu = mu2,t = t2,t0 = ts,f = f4) + 
            log1minusu(la = la2,mu = mu2,t = t2,t0 = ts,f = f4) - 
            2 * log(1 - exp
            (logu(la = la1,mu = mu1,t = ts,t0 = t0,f = f1) +
             logP0(la = la1,mu = mu1,t = t2,t0 = ts,f = f3))) -
            2 * log(1 - exp  
            (logu(la = la1,mu = mu1,t = ts,t0 = t1,f = f2) +
             logP0(la = la1,mu = mu1,t = t2,t0 = ts,f = f3))) - 
            log(2 - exp
            (logu(la = la1,mu = mu1,t = ts,t0 = t0,f = f1) + 
             logP0(la = la1,mu = mu1,t = t2,t0 = ts,f = f3)) - exp
            (logu(la = la1,mu = mu1,t = ts,t0 = t1,f = f2) + 
             logP0(la = la1,mu = mu1,t = t2,t0 = ts,f = f3)))
  return(loglik)
}
