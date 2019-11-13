---
title: "Combining single- and repeated-visit occupancy models to make the best of monitoring surveys: Appendix III"
date: "`r Sys.Date()`"
author: "Valentin Lauret, Hélène Labach, Matthieu Authier, Olivier Gimenez"
---

In this Appendix, we present how we built the multi-method single-visit (SV) occupancy models analysing the aerial surveys and the at-sea surveys with as single-visit occupancy models.


# Single visit occupancy

In SV occupancy modelling, _J = 1_ . Then, each detection `y[i]` at site _i_ is taken in a Bernoulli distribution with parameter `mu[i]`.  

`y[i]~dbern(mu[i])`  

The `mu` parameter is defined as `mu[i] = p[i] * z[i]`.  
`p[i]` is the probability of detecting the species at site i.

The latent ecological state remains unchanged from a _repeated-visit_ occupancy model. The occupancy status, of site i is taken in a Bernoulli distribution of parameter `psi`  

`z[i]~dbern(psi[i])`.

Then `psi` represents the occupancy probability of site _i_.   

`psi` is modelled as a logistic regression of the environmental covariates: bathymetry and SST.`p` is modeled as a logistic regression of sampling effort. `alpha.psi`, `alpha.p`, `beta.sst`, `beta.bathy`, and `beta.eff` are unknown parameters:  

 `logit(psi[i]) = alpha.psi + beta.sst * SST[i] + beta.bathy * BATHY[i]`  
 
 `logit(p[i]) = alpha.p + beta.eff * EFF[i]`
 
 `SST`, `BATHY`, and `EFF` are the dataframes that contain the numerical values of respectively SST, bathymetry, and sampling effort. In SV occupancy modelling the sampling effort at site _i_, i.e. `EFF[i]` corresponds to the transect length (in km) prospected by the monitoring method at site _i_ during the whole duration of the monitoring protocol.

# Single-method SV modelling

## Aerial SV modelling  

We considered the detection/non-detection matrix `y` (1 column) that contains 0/1 for every site i of the study area.

There are 2 trials with Bernoulli distributions, with their associated `p`.


  * `y[i]~dbern(mu[i])` , with `mu[i] = z[i] * p[i]`  
  
With  ` logit(p[i]) = alpha.p + beta.eff*eff[i]`

The latent ecological state process remains unchanged from the classical occupancy modelling.
  
Hereafter, you would find the JAGS formulation of this occupancy model. 
 
```{r eval = FALSE}
# Specify model in BUGS language
library(R2jags)

sink("Aerial_SV.jags")
cat("
    model {
    # priors
    alpha.psi ~ dnorm(0,0.444) # occupancy intercept
    alpha.p ~ dnorm(0,0.444) # detection intercept 
    beta.sst ~ dnorm(0,0.444) # slope sst effect
    beta.bathy ~ dnorm(0,0.444)	# slope bathy effect
    beta.eff ~ dnorm(0,0.444)	# slope survey effort effect
    
    for (i in 1:nsite){
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- lpsi[i]
    lpsi[i] <- alpha.psi + beta.sst * SST[i] + beta.bathy * BATHY[i]
    
    mu.p[i] <- z[i] * p[i] 
    logit(p[i]) <- lp[i] 
    lp[i] <- alpha.p + beta.eff * EFF[i]
    y[i] ~ dbern(mu.p[i])
    } #i
    }
    ", fill = TRUE)
sink()
```

# Data & run

Hereafter, we show how we loaded and format the data for occupancy modelling

```{r}
load("Lauret_et_al_Appendix.rdata")


out <- SAMM_sv
library(mcmcplots)
denplot(out, c("alpha.psi","beta.bathy", "beta.sst"))
traplot(out,"beta.sst")
out.m <- as.mcmc(out)
effectiveSize(out.m[[1]][,'beta.bathy'])
str(out.m)
names(out.m)
out.m[[1]][,'beta.bathy']
```
# References 

Article Ref
