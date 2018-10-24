\name{bd_RS_loglik_BAMM}
\alias{bd_RS_loglik_BAMM}
\title{Loglikelihood for a constant-rate birth-death process with shift in the parameters on specific lineages
at times tshift }
\description{
  This function computes loglikelihood of a constant-rate birth-death model for a given set of branching times and parameter values where the parameters are allowed to shift at times tshift.
}
\usage{bd_RS_loglik_BAMM(BAMMtable)}
\arguments{
  \item{BAMMtable}{A table containing a full description of a birth-death tree with rate shifts.
Each row is a branch segment governed by a particular (lambda, mu) parameterization. A branch that has a shift will be governed by two parameter sets and each segment has its own row.
Columns of the table: 
node        = the (tipwards) node index in ape format for the branch segment
start_time  = absolute start time (from root) of segment for particular lambda/mu params
end_time    = absolute end time (from root) for segment
shift_index = reference index of rate regime that links back to original eventdata matrix
lam1        = lambda, the speciation rate
mu1         = mu, the extinction rate
}}
\value{The loglikelihood}
\references{
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi: 10.1086/667574
}
\author{Rampal S. Etienne & Bart Haegeman }
\seealso{
 \code{\link{bd_loglik}},
 \code{\link{bd_RS_loglik}},
 \code{\link{dd_KI_loglik}}
}
\examples{
dd0 <- bd_RS_loglik(
   pars1 = list(c(0.2,0.05),c(1,0.1)),
   pars2 = NULL,
   brts = list(c(20,20,18,16,12,10,6,4),c(9.5,9,8,7,5,2)),
   tshift = list(c(9.5),NULL)
   )
dd1 <- dd_KI_loglik(
    pars1 = c(0.2,0.05,Inf,1,0.1,Inf,9.5),
    pars2 = c(500,1,0,10,1,2),
    brtsM = c(20,18,16,12,10,6,4),
    brtsS = c(9,8,7,5,2),
    missnumspec = 0,
    method = 'ode45'
    )
BAMMtable <- matrix(0,nrow = 25,ncol = 6)
BAMMtable[,2] <- -c(20,18,18,16,16,12,12,10,10,6,6,4,4,20,9.5,9,9,8,8,7,7,5,5,2,2)
BAMMtable[,3] <- -c(18,0,16,0,12,0,10,0,6,0,4,0,0,9.5,9,0,8,0,7,0,5,0,2,0,0)
BAMMtable[,5] <- c(rep(0.2,14),rep(1,11))
BAMMtable[,6] <- c(rep(0.05,14),rep(0.1,11))
dd2 <- bd_RS_loglik_BAMM(BAMMtable)
print(dd0)
print(dd1)
print(dd2)
bd0 <- bd_loglik(
   brts = c(20,18,16,12,10,6,4,9,8,7,5,2),
   pars1 = c(0.2,0.05),
   pars2 = c(0,0,1,0,2,100),
   missnumspec = 0
   )
bd1 <- bd_RS_loglik(
    pars1 = list(c(0.2,0.05),c(0.2,0.05)),
    pars2 = NULL,
    brts = list(c(20,20,18,16,12,10,6,4),c(9.5,9,8,7,5,2)),
    tshift = list(c(9.5),NULL)
    )
BAMMtable[,5] <- rep(0.2,25)
BAMMtable[,6] <- rep(0.05,25)
bd2 <- bd_RS_loglik_BAMM(BAMMtable)
print(bd0)
print(bd1)
print(bd2)
}
\keyword{models}
