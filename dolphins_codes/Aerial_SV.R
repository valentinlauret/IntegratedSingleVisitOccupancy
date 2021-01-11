### SAMM single visit

# Specify model in BUGS language
library(R2jags)

sink("SAMM_SV.jags")
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


#setwd("/home/ogimenez/valentin") # For the cluster

#setwd("/Users/valentinlauret/Documents/Google Drive/These/Work/Occupancy/Single visit")

load("SAMM_ready2run.rdata")

head(samm.effort)
samm.effort$eff.ind[ is.na(samm.effort$eff.ind )== T] <- 0


# Process data
ls()
dim(y)
str(y)
surveyed <- samm.effort$eff.ind # whether a site was surveyed or not
mask <- surveyed # sites that were not surveyed have 0s
y <- y[mask!=0,] # filter out sites that were not surveyed at all throughout the study period
effort <- samm.effort$eff.tot # survey effort
effort <- effort[mask!=0] # get rid of effort for sites that were never surveyed
surveyed <- surveyed[mask!=0] # update the survey indicator
dim(y)
dim(effort)
det <- y$obs # detection and non-detections

# Bundle data
# Build dataset
str(win.data <- list(
  y = det, 
  EFF = effort,
  BATHY = as.numeric(y$bathy.sc),
  SST = as.numeric(y$sst.sc), 
  nsite = nrow(y)))

# Initial values
# Observed occurrence as inits for z, prend le max (ce qui est 1) des quatres occasions et sort du coup la vraie occupation indépendemment de la proba d'observation, si par exemple un site est 0101 au cours des quatres occasions il en sortira 1, et 0000 il en sortira 0. 
# Initial values
inits <- function() {list(
  z = det, 
  alpha.psi = rnorm(1, 0, 1.5),
  alpha.p = rnorm(1, 0, 1.5),
  beta.sst = rnorm(1, 0, 1.5),
  beta.bathy = rnorm(1, 0, 1.5),
  beta.eff = rnorm(1, 0, 1.5))}

# Parameters to be monitored
params <- c("alpha.psi",
            "alpha.p",
            "beta.eff",
            "beta.bathy", 
            "beta.sst","z") 

# MCMC settings
ni <- 20000
nt <- 5
nb <- 1000
nc <- 2

# Call JAGS from R
#ptm <- proc.time()
out <- jags(win.data, params, "SAMM_SV.jags", inits = inits, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
#x <- proc.time() -  ptm

SAMM_sv <- out

save(SAMM_sv, file="RES_SAMM_sv.rdata")
#github
save(win.data,inits,params,file ="AtSea_SV.Rdata")
