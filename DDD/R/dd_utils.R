brts2phylo <- function(times,root=FALSE,tip.label=NULL)
{
# This code is taken from the package TESS by Sebastian Hoehna, where the function is called tess.create.phylo
# It takes a set of branching times and adds a random topology.
  times = sort(times)
  n <- as.integer(length(times))+1
  if ( root ) {
    n <- n-1
  }
  nbr <- 2*n - 2

  # create the data types for edges and edge-lengths
  edge <- matrix(NA, nbr, 2)
  edge.length <- numeric(nbr)

  h <- numeric(2*n - 1) # initialized with 0's
  pool <- 1:n
  # VERY VERY IMPORTANT: the root MUST have index n+1 !!!
  nextnode <- 2L*n - 1L
  if ( n > 1) {
    for (i in 1:(n - 1)) {
      # sample two nodes that have no parent yet
      y <- sample(pool, size = 2)
      # compute the edge indices (we just order the edges from 1 to 2n-2)
      ind <- (i - 1)*2 + 1:2
      # set the source node of the new edges (i.e. the new internal node)
      edge[ind, 1] <- nextnode
      # set the destination of the new edges (i.e. the two sampled nodes)
      edge[ind, 2] <- y
      # compute the edge length from the difference between the node heights (child <-> parent)
      edge.length[ind] <- times[i] - h[y]
      # store the node height of this new internal node
      # we cannot use times because then we would get into trouble with the indices and we would need to check for tip nodes ...
      h[nextnode] <- times[i]
      # reset the pool of available nodes to merge
      pool <- c(pool[! pool %in% y], nextnode)
      # increase the node index counter
      nextnode <- nextnode - 1L
    }
  }

  phy <- list(edge = edge, edge.length = edge.length)
  if (is.null(tip.label))
    tip.label <- paste("t", 1:n, sep = "")
  phy$tip.label <- sample(tip.label)
  phy$Nnode <- n - 1L

  if ( root ) {
    phy$root.edge <- times[n] - times[n-1]
    phy$root <- times[n] - times[n-1]
  }

  class(phy) <- "phylo"

  phy <- reorder(phy)
  ## to avoid crossings when converting with as.hclust:
  phy$edge[phy$edge[, 2] <= n, 2] <- 1:n

  return(phy)
}

conv = function(x,y)
{
   lx = length(x)
   ly = length(y)
   lxy = length(x) + length(y)
   x = c(x,rep(0,lxy - lx))
   y = c(y,rep(0,lxy - ly))
   cvxy = rep(0,lxy)
   for(i in 2:lxy)
   {
      cvxy[i] = crossprod(x[(i-1):1],y[1:(i-1)])
   }
   return(cvxy[2:lxy])
}

flavec = function(ddep,la,mu,K,r,lx,kk,n0)
{
   nn = (0:(lx - 1)) + kk
   if(ddep == 1 | ddep == 5)
   {
       lavec = pmax(rep(0,lx),la - 1/(r + 1) * (la - mu)/K * nn)
   }
   if(ddep == 1.3)
   {
       lavec = pmax(rep(0,lx),la * (1 - nn/K))
   }
   if(ddep == 1.4)
   {
     lavec = pmax(rep(0,lx),la * nn/(nn + K))
   }
   if(ddep == 1.5)
   {
       lavec = pmax(rep(0,lx),la * nn/K * (1 - nn/K))
   }
   if(ddep == 2 | ddep == 2.1 | ddep == 2.2)
   {
       x = -(log(la/mu)/log(K + n0))^(ddep != 2.2)
       lavec = pmax(rep(0,lx),la * (nn + n0)^x)
   }
   if(ddep == 2.3)
   {
       x = K
       lavec = pmax(rep(0,lx),la * (nn + n0)^x)
   }
   if(ddep == 3 | ddep == 4 | ddep == 4.1 | ddep == 4.2)
   {
       lavec = la * rep(1,lx)
   }
   return(lavec)
}

L2phylo = function(L,dropextinct = T)
# makes a phylogeny out of a matrix with branching times, parent and daughter species, and extinction times
{
   L = L[order(abs(L[,3])),1:4]
   age = L[1,1]
   L[,1] = age - L[,1]
   L[1,1] = -1
   notmin1 = which(L[,4] != -1)
   L[notmin1,4] = age - L[notmin1,4]
   if(dropextinct == T)
   {
      sall = which(L[,4] == -1)
      tend = age
   } else {
      sall = which(L[,4] >= -1)
      tend = (L[,4] == -1) * age + (L[,4] > -1) * L[,4]
   }
   L = L[,-4]
   linlist = cbind(data.frame(L[sall,]),paste("t",abs(L[sall,3]),sep = ""),tend)
   linlist[,4] = as.character(linlist[,4])
   names(linlist) = 1:5
   done = 0
   while(done == 0)
   {
      j = which.max(linlist[,1])
      daughter = linlist[j,3]
      parent = linlist[j,2]
      parentj = which(parent == linlist[,3])
      parentinlist = length(parentj)
      if(parentinlist == 1)
      {
         spec1 = paste(linlist[parentj,4],":",linlist[parentj,5] - linlist[j,1],sep = "")
         spec2 = paste(linlist[j,4],":",linlist[j,5] - linlist[j,1],sep = "")
         linlist[parentj,4] = paste("(",spec1,",",spec2,")",sep = "")
         linlist[parentj,5] = linlist[j,1]
         linlist = linlist[-j,]
      } else {
         #linlist[j,1:3] = L[abs(as.numeric(parent)),1:3]
         linlist[j,1:3] = L[which(L[,3] == parent),1:3]
      }
      if(nrow(linlist) == 1) { done = 1 }
   }
   linlist[4] = paste(linlist[4],":",linlist[5],";",sep = "")
   phy = read.tree(text = linlist[1,4])
   tree = as.phylo(phy)
   return(tree)
}

phylo2L = function(phy)
{
  emdata = phy
  # compute the relative branching times 
  brt = branching.times(emdata)
  if(min(brt) < 0)
  {
    brt = brt + abs(min(brt))
  }
  # number of species including extinct species.
  num.species = emdata$Nnode+1
  brt_preL = c(brt[emdata$edge[,1] - length(emdata$tip.label)])
  # check if the relative branching times are equal to the real branching times.
  # if not correct it to the real branching times.
  if(min(brt_preL) == 0)
  {
    correction = max(emdata$edge.length[which(brt_preL==0)])
    brt_preL = brt_preL+correction
  }
  # preliminary L table
  pre.Ltable = cbind(brt_preL,emdata$edge,emdata$edge.length,brt_preL-emdata$edge.length)
  # identify the extant species and the extinct species
  extantspecies.index = pre.Ltable[which(pre.Ltable[,5]<=1e-10),3]
  tipsindex = c(1:num.species)
  extinct.index3 = subset(tipsindex,!(tipsindex %in% extantspecies.index))
  # assigen the extinct species with extinct times; the extant species with -1 and
  # the internal nodes with 0.
  eeindicator = matrix(0,length(emdata$edge.length),1)
  eeindicator[match(extantspecies.index,pre.Ltable[,3])]=-1
  ext.pos = match(extinct.index3,pre.Ltable[,3])
  eeindicator[ext.pos] = pre.Ltable[ext.pos,5]
  pre.Ltable = cbind(pre.Ltable,eeindicator)
  
  sort.L = pre.Ltable[order(pre.Ltable[,1],decreasing = TRUE),]
  nodesindex = unique(emdata$edge[,1])
  L = sort.L
  realL = NULL
  do = 0
  while(do == 0){
    j = which.min(L[,3])
    daughter = L[j,3]
    parent = L[j,2]
    if(parent %in% nodesindex)
    {
      L[which(L[,2]==parent),2] = daughter
      if(length(which(L[,3] == parent)) == 0){
        realL = rbind(realL,L[j,],row.names = NULL)
        L = L[-j,,drop=FALSE]
      } else {
        L[which(L[,3] == parent),6] = L[j,6]
        L[which(L[,3] == parent),3] = daughter
        L = L[-j,,drop=FALSE]
      }
    } else {
      realL = rbind(realL,L[j,],row.names = NULL)
      L = L[-j,,drop=FALSE]
    }
    
    if(nrow(L) == 0){
      do = 1
    }
  }
  realL = realL[order(realL[,1],decreasing = T),]
  L = realL[,c(1,2,3,6)]
  
  daughter.index = L[,3]
  daughter.realindex = c(1:nrow(L))
  parent.index = L[,2]
  parent.realindex = match(parent.index, daughter.index)
  
  L[,2] = parent.realindex
  L[,3] = daughter.realindex
  L[1,2] = 0
  L[1,3] = -1
  L[2,2] = -1
  for(i in c(2:nrow(L)))
  {
    if(L[i-1,3] < 0){
      mrows = which(L[,2] == abs(L[i-1,3]))
      L[mrows,2] = L[i-1,3]
      L[mrows,3] = -1 * L[mrows,3]
    }
  }
  dimnames(L) = NULL
  return(L)
}

L2brts = function(L,dropextinct = T)
# makes a phylogeny out of a matrix with branching times, parent and daughter species, and extinction times
{
   brts = NULL
   L = L[order(abs(L[,3])),1:4]
   age = L[1,1]
   L[,1] = age - L[,1]
   L[1,1] = -1
   notmin1 = which(L[,4] != -1)
   L[notmin1,4] = age - L[notmin1,4]
   if(dropextinct == T)
   {
      sall = which(L[,4] == -1)
      tend = age
   } else {
      sall = which(L[,4] >= -1)
      tend = (L[,4] == -1) * age + (L[,4] > -1) * L[,4]
   }
   L = L[,-4]
   linlist = cbind(data.frame(L[sall,]),paste("t",abs(L[sall,3]),sep = ""),tend)
   linlist[,4] = as.character(linlist[,4])
   names(linlist) = 1:5
   done = 0
   while(done == 0)
   {
      j = which.max(linlist[,1])
      daughter = linlist[j,3]
      parent = linlist[j,2]
      parentj = which(parent == linlist[,3])
      parentinlist = length(parentj)
      if(parentinlist == 1)
      {
         spec1 = paste(linlist[parentj,4],":",linlist[parentj,5] - linlist[j,1],sep = "")
         spec2 = paste(linlist[j,4],":",linlist[j,5] - linlist[j,1],sep = "")
         linlist[parentj,4] = paste("(",spec1,",",spec2,")",sep = "")
         linlist[parentj,5] = linlist[j,1]
         brts = c(brts,linlist[j,1])
         linlist = linlist[-j,]
      } else {
         #linlist[j,1:3] = L[abs(as.numeric(parent)),1:3]
         linlist[j,1:3] = L[which(L[,3] == parent),1:3]
      }
      if(nrow(linlist) == 1) { done = 1 }
   }
   #linlist[4] = paste(linlist[4],":",linlist[5],";",sep = "")
   brts = rev(sort(age - brts))
   return(brts)
}

roundn = function(x, digits = 0)
{
    fac = 10^digits
    n = trunc(fac * x + 0.5)/fac
    return(n)
}

sample2 = function(x,size,replace = FALSE,prob = NULL)
{
    if(length(x) == 1)
    { 
        x = c(x,x)
        prob = c(prob,prob)
    }
    sam = sample(x,size,replace,prob)
    return(sam)
}

simplex = function(fun,trparsopt,optimpars,...)
{
  numpar = length(trparsopt)
  reltolx = optimpars[1]
  reltolf = optimpars[2]
  abstolx = optimpars[3]
  maxiter = optimpars[4]

  ## Setting up initial simplex
  v = t(matrix(rep(trparsopt,each = numpar + 1),nrow = numpar + 1))
  for(i in 1:numpar)
  {
      parsoptff = 1.05 * untransform_pars(trparsopt[i])
      trparsoptff = transform_pars(parsoptff)
      fac = trparsoptff/trparsopt[i]
      if(v[i,i + 1] == 0)
      {
         v[i,i + 1] = 0.00025
      } else {
         v[i,i + 1] = v[i,i + 1] * min(1.05,fac)
      }
  }
  
  fv = rep(0,numpar + 1)
  for(i in 1:(numpar + 1))
  {
     fv[i] = -fun(trparsopt = v[,i], ...)
  }
  
  how = "initial"
  itercount = 1
  string = itercount
  for(i in 1:numpar)
  {
     string = paste(string, untransform_pars(v[i,1]), sep = " ")
  }
  string = paste(string, -fv[1], how, "\n", sep = " ")
  cat(string)
  flush.console()
  
  tmp = order(fv)
  if(numpar == 1)
  {
     v = matrix(v[tmp],nrow = 1,ncol = 2)
  } else {
     v = v[,tmp]
  }
  fv = fv[tmp]
  
  ## Iterate until stopping criterion is reached
  rh = 1
  ch = 2
  ps = 0.5
  si = 0.5
  
  v2 = t(matrix(rep(v[,1],each = numpar + 1),nrow = numpar + 1))
  
  while(itercount <= maxiter & ( ( is.nan(max(abs(fv - fv[1]))) | (max(abs(fv - fv[1])) - reltolf * abs(fv[1]) > 0) ) + ( (max(abs(v - v2) - reltolx * abs(v2)) > 0) | (max(abs(v - v2)) - abstolx > 0) ) ) )
  { 
     ## Calculate reflection point
  
     if(numpar == 1)
     {
         xbar = v[1]
     } else {
         xbar = rowSums(v[,1:numpar])/numpar
     }
     xr = (1 + rh) * xbar - rh * v[,numpar + 1]
     fxr = -fun(trparsopt = xr, ...)
   
     if(fxr < fv[1])
     {
         ## Calculate expansion point
         xe = (1 + rh * ch) * xbar - rh * ch * v[,numpar + 1]
         fxe = -fun(trparsopt = xe, ...)
         if(fxe < fxr)
         {
             v[,numpar + 1] = xe
             fv[numpar + 1] = fxe
             how = "expand"
         } else {
             v[,numpar + 1] = xr
             fv[numpar + 1] = fxr
             how = "reflect"
         }
     } else {
         if(fxr < fv[numpar])
         {      
             v[,numpar + 1] = xr
             fv[numpar + 1] = fxr
             how = "reflect"
         } else {
             if(fxr < fv[numpar + 1])
             {
                ## Calculate outside contraction point
                xco = (1 + ps * rh) * xbar - ps * rh * v[,numpar + 1]
                fxco = -fun(trparsopt = xco, ...)
                if(fxco <= fxr)
                {
                   v[,numpar + 1] = xco
                   fv[numpar + 1] = fxco            
                   how = "contract outside"
                } else {
                   how = "shrink"
                }
             } else {
                ## Calculate inside contraction point
                xci = (1 - ps) * xbar + ps * v[,numpar + 1]
                fxci = -fun(trparsopt = xci, ...)
                if(fxci < fv[numpar + 1])
                {  
                   v[,numpar + 1] = xci
                   fv[numpar + 1] = fxci
                   how = "contract inside"
                } else {
                   how = "shrink"
                }
             }
             if(how == "shrink")
             {
                 for(j in 2:(numpar + 1))
                 {
  
                     v[,j] = v[,1] + si * (v[,j] - v[,1])
                     fv[j] = -fun(trparsopt = v[,j], ...)
                 }
             }
         }
     }
     tmp = order(fv)
     if(numpar == 1)
     {
        v = matrix(v[tmp],nrow = 1,ncol = 2)
     } else {
        v = v[,tmp]
     }
     fv = fv[tmp]
     itercount = itercount + 1
     string = itercount;
     for(i in 1:numpar)
     {
         string = paste(string, untransform_pars(v[i,1]), sep = " ")
     }
     string = paste(string, -fv[1], how, "\n", sep = " ")
     cat(string)
     flush.console()
     v2 = t(matrix(rep(v[,1],each = numpar + 1),nrow = numpar + 1))
  }
  if(itercount < maxiter)
  {
     cat("Optimization has terminated successfully.","\n")
  } else {
     cat("Maximum number of iterations has been exceeded.","\n")
  }
  out = list(par = v[,1], fvalues = -fv[1], conv = as.numeric(itercount > maxiter))
  invisible(out)
}

optimizer = function(optimmethod = 'simplex',optimpars = c(1E-4,1E-4,1E-6,1000),fun,trparsopt, ...)
{
  if(optimmethod == 'simplex')
  {
    out = simplex(fun = fun,trparsopt = trparsopt,optimpars = optimpars,...)
  }
  if(optimmethod == 'subplex')
  {
    minfun = function(fun,trparsopt,...)
    {           
      return(-fun(trparsopt = trparsopt,...))
    }
    
    trySubplex <- function(trparsopt = trparsopt, optimpars = optimpars, ...){
      tryCatch(
        # main function to run
        subplex::subplex(
          par = trparsopt,
          fn = minfun,
          control = list(abstol = optimpars[3],
                         reltol = optimpars[1],
                         maxit = optimpars[4]
          ),
          fun = fun, ...
        ),
        error = function(e){
          print("subplex has encountered an error")
          # output in case of error
          list(
            par = rep(-1, length(trparsopt)),
            value = -Inf,
            convergence = -3 # other values like -1 are already used by subplex for other situations
          )
        }
      )
    }
    out <- trySubplex(trparsopt = trparsopt, optimpars = optimpars, ...)
    
    out = list(par = out$par, fvalues = -out$value, conv = out$convergence)
  }
  return(out)
}

transform_pars <- function(pars)
{
  trpars1 <- sign(pars) * pars/(sign(pars) + pars);
  trpars1[which(pars == 0)] <- 0;
  trpars1[which(pars == -Inf)] <- -1;
  trpars1[which(pars == Inf)] <- 1;
  return(trpars1);
}

untransform_pars <- function(trpars1)
{
  pars <- sign(trpars1) * trpars1/(sign(trpars1) - trpars1);
  pars[which(trpars1 == 0)] <- 0;
  pars[which(trpars1 == 1)] <- Inf;
  pars[which(trpars1 == -1)] <- -Inf;
  return(pars)
}
