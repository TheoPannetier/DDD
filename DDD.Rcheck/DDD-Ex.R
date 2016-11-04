pkgname <- "DDD"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "DDD-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('DDD')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("L2phylo")
### * L2phylo

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: L2phylo
### Title: Function to convert a table with speciation and extinction
###   events to a phylogeny
### Aliases: L2phylo
### Keywords: models

### ** Examples

sim = dd_sim(c(0.2,0.1,20),10)
phy = L2phylo(sim$L)
plot(phy)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("L2phylo", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bd_ML")
### * bd_ML

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bd_ML
### Title: Maximization of the loglikelihood under the
###   diversity-independent, possibly time-dependent diversification model
### Aliases: bd_ML
### Keywords: models

### ** Examples

cat("Estimating parameters for a set of branching times brts with the default settings:")
brts = 1:20
bd_ML(brts = brts, cond = 1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bd_ML", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bd_loglik")
### * bd_loglik

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bd_loglik
### Title: Loglikelihood for diversity-independent diversification model
### Aliases: bd_loglik
### Keywords: models

### ** Examples
bd_loglik(pars1 = c(0.5,0.1), pars2 = c(0,1,1,0,2), brts = 1:10, missnumspec = 0) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bd_loglik", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dd_KI_ML")
### * dd_KI_ML

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dd_KI_ML
### Title: Maximization of the loglikelihood under a diversity-dependent
###   diversification model with decoupling of a subclade's diversication
###   dynamics from the main clade's dynamics
### Aliases: dd_KI_ML
### Keywords: models

### ** Examples

cat("This will estimate parameters for two sets of branching times brtsM, brtsS\n")
cat("without conditioning.\n")
cat("The tolerance of the optimization is set high so runtime is fast in this example.\n")
cat("In real applications, use the default or more stringent settins for tol.\n")
brtsM = 4:10
brtsS = seq(0.1,3.5,0.7)
tsplit = 5
dd_KI_ML(brtsM = brtsM, brtsS = brtsS, tsplit = tsplit, idparsopt = c(1:3,6,7),
          initparsopt = c(0.885, 2e-14, 6.999, 6.848, 4.001), idparsfix = NULL, parsfix = NULL,
          idparsnoshift = c(4,5), cond = 0, tol = c(3E-1,3E-1,3E-1))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dd_KI_ML", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dd_KI_loglik")
### * dd_KI_loglik

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dd_KI_loglik
### Title: Loglikelihood for diversity-dependent diversification models
###   with decoupling of a subclade from a main clade at time t = t_d
### Aliases: dd_KI_loglik
### Keywords: models

### ** Examples

pars1 = c(0.25,0.12,25.51,1.0,0.16,8.61,9.8)
pars2 = c(200,1,0,18.8,1,2)
missnumspec = 0
brtsM = c(25.2,24.6,24.0,22.5,21.7,20.4,19.9,19.7,18.8,17.1,15.8,11.8,9.7,8.9,5.7,5.2)
brtsS = c(9.6,8.6,7.4,4.9,2.5)
dd_KI_loglik(pars1,pars2,brtsM,brtsS,missnumspec,method = 'ode45')



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dd_KI_loglik", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dd_KI_sim")
### * dd_KI_sim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dd_KI_sim
### Title: Function to simulate a key innovation in macro-evolution with
###   the innovative clade decoupling from the diversity-dependent
###   diversification dynamics of the main clade
### Aliases: dd_KI_sim
### Keywords: models

### ** Examples
 dd_KI_sim(c(0.2,0.1,20,0.1,0.05,30,4),10) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dd_KI_sim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dd_ML")
### * dd_ML

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dd_ML
### Title: Maximization of the loglikelihood under a diversity-dependent
###   diversification model
### Aliases: dd_ML
### Keywords: models

### ** Examples

cat("Estimating the intrinsic speciation rate lambda and the carrying capacity K")
cat("for a fixed extinction rate of 0.1, conditioning on clade survival and two missing species:")
brts = 1:5
dd_ML(brts = brts,initparsopt = c(1.3078,7.4188), idparsopt = c(1,3), parsfix = 0.1,
      cond = 1, missnumspec = 2, tol = c(1E-3,1E-3,1E-4), optimmethod = 'simplex')



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dd_ML", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dd_MS_ML")
### * dd_MS_ML

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dd_MS_ML
### Title: Maximization of the loglikelihood under a diversity-dependent
###   diversification model with decoupling of a subclade's diversication
###   dynamics from the main clade's dynamics
### Aliases: dd_MS_ML
### Keywords: models

### ** Examples

cat("This will estimate parameters for two sets of branching times brtsM, brtsS\n")
cat("without conditioning.\n")
cat("The tolerance of the optimization is set high so runtime is fast in this example.\n")
cat("In real applications, use the default or more stringent settins for tol.\n")
brtsM = 4:10
brtsS = seq(0.1,3.5,0.7)
tsplit = 5
dd_MS_ML(brtsM = brtsM, brtsS = brtsS, tsplit = tsplit, idparsopt = c(1:3,6),
          initparsopt = c(0.885, 2e-14, 10, 4.001), idparsfix = NULL, parsfix = NULL,
          idparsnoshift = c(4,5), cond = 0, tol = c(3E-1,3E-1,3E-1))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dd_MS_ML", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dd_MS_loglik")
### * dd_MS_loglik

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dd_MS_loglik
### Title: Loglikelihood for macro-evolutionary succession under
###   diversity-dependent diversification with the key innovation at time t
###   = t_d
### Aliases: dd_MS_loglik
### Keywords: models

### ** Examples

pars1 = c(0.2,0.1,40,1.0,0.1,9.8)
pars2 = c(200,1,0,18.8,1,2)
missnumspec = 0
brtsM = c(25.2,24.6,24.0,22.5,21.7,20.4,19.9,19.7,18.8,17.1,15.8,11.8,9.7,8.9,5.7,5.2)
brtsS = c(9.6,8.6,7.4,4.9,2.5)
dd_MS_loglik(pars1,pars2,brtsM,brtsS,missnumspec,methode = 'ode45')



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dd_MS_loglik", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dd_MS_sim")
### * dd_MS_sim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dd_MS_sim
### Title: Function to simulate the macro-evolutionary succession process
###   assuming diversity-dependent diversification
### Aliases: dd_MS_sim
### Keywords: models

### ** Examples
 dd_MS_sim(c(0.2,0.1,20,0.1,0.05,4),10) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dd_MS_sim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dd_SR_ML")
### * dd_SR_ML

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dd_SR_ML
### Title: Maximization of the loglikelihood under a diversity-dependent
###   diversification model with a shift in the parameters
### Aliases: dd_SR_ML
### Keywords: models

### ** Examples

cat("This will estimate parameters for a sets of branching times brts without conditioning.\n")
cat("The tolerance of the optimization is set ridiculously high to make runtime fast.\n")
cat("In real applications, use the default or more stringent settings for tol.\n")
brts = 1:10
dd_SR_ML(brts = brts, initparsopt = c(0.4581, 1E-6, 17.69, 11.09, 8.9999), idparsopt = c(1:3,6,7),
         idparsfix = NULL, parsfix = NULL, idparsnoshift = c(4,5), cond = 0,
         tol = c(1E-1,1E-1,1E-1),optimmethod = 'simplex'
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dd_SR_ML", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dd_SR_loglik")
### * dd_SR_loglik

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dd_SR_loglik
### Title: Loglikelihood for diversity-dependent diversification models
###   with a shift in the parameters at time t = tshift
### Aliases: dd_SR_loglik
### Keywords: models

### ** Examples
dd_SR_loglik(pars1 = c(0.2,0.1,50,0.2,0.1,70,5), pars2 = c(100,1,1,1,0,2),
   brts = 1:10, missnumspec = 0) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dd_SR_loglik", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dd_SR_sim")
### * dd_SR_sim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dd_SR_sim
### Title: Function to simulate the diversity-dependent diversification
###   process with a shift in one or more of the parameters
### Aliases: dd_SR_sim
### Keywords: models

### ** Examples
 dd_SR_sim(c(0.2,0.1,20,0.2,0.1,40,5),10) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dd_SR_sim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dd_loglik")
### * dd_loglik

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dd_loglik
### Title: Loglikelihood for diversity-dependent diversification models
### Aliases: dd_loglik
### Keywords: models

### ** Examples
dd_loglik(pars1 = c(0.5,0.1,100), pars2 = c(100,1,1,1,0,2), brts = 1:10, missnumspec = 0) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dd_loglik", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dd_sim")
### * dd_sim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dd_sim
### Title: Function to simulate the diversity-dependent diversification
###   process
### Aliases: dd_sim
### Keywords: models

### ** Examples
 dd_sim(c(0.2,0.1,20),10) 


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dd_sim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("optimizer")
### * optimizer

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: optimizer
### Title: Carries out optimization (finding a minimum)
### Aliases: optimizer
### Keywords: models

### ** Examples

cat("No examples")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("optimizer", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("roundn")
### * roundn

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: roundn
### Title: Rounds up in the usual manner
### Aliases: roundn
### Keywords: models

### ** Examples

round(2.5)
roundn(2.5)
round(3.5)
roundn(3.5)
round(2.65,digits = 1)
roundn(2.65,digits = 1)
round(2.75,digits = 1)
roundn(2.75,digits = 1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("roundn", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sample2")
### * sample2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sample2
### Title: Takes samples in the usual manner
### Aliases: sample2
### Keywords: models

### ** Examples

sample(x = 10,size = 5,replace = TRUE)
sample2(x = 10,size = 5,replace = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sample2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
