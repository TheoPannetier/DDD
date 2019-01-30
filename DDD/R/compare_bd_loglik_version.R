compare_bd_loglik_version <- function(brts, pars1, pars2, missnumspec = 0, methode  ='lsoda', tol = 0.0)
{
  # loglik with dd_ode_FORTRAN()
  loglik_3.8 <- bd_loglik(
    pars1 = pars1,
    pars2 = pars2,
    missnumspec = missnumspec,
    methode = methode
    )
  # loglik with ode()
  loglik_3.2 <- DDD:::bd_loglik_3.2(
    pars1 = pars1,
    pars2 = pars2,
    missnumspec = missnumspec,
    methode = methode
  )
  
  is_loglik_identical = (abs(loglik_3.8 - loglik_3.2) <= tol)
  return(is_loglik_identical)
}