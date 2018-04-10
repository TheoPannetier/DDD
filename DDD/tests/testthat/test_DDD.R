context("test_DDD")

test_that("the Fortran code and the R code give the same answer", {
  pars1 = c(0.8,0.2,140)
  pars2 = c(100,1,1,0,0,2)
  brts = 1:30
  missnumspec = 0
  methode = 'lsoda'
  r1 <- dd_loglik(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,rhs_func_name = 'dd_loglik_rhs',methode = methode)
  r2 <- dd_loglik(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,rhs_func_name = 'dd_loglik_rhs_FORTRAN',methode = methode)
  testthat::expect_equal(r1,r2,tolerance = .00001)
  testthat::expect_equal(-109.177029,r1,tolerance = .000001)
  
  r3 <- dd_SR_loglik(pars1 = c(0.2,0.1,50,0.2,0.1,70,5), pars2 = c(100,1,1,1,0,2), brts = 1:10, missnumspec = 0)
  testthat::expect_equal(-27.37304,r3,tolerance = .000001)
})


