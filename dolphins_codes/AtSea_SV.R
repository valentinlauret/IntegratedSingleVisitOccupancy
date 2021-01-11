#############################
## AtSea SV occupancy model 
#############################


############################# 
# Specify model in BUGS language
#############################

library(R2jags)


sink("AtSea_SV.jags")
cat("
    model {

    # priors
   
    alpha.psi ~ dnorm(0,0.444) # occupancy intercept
    alpha.p ~ dnorm(0,0.444) # detection intercept 
    beta.sst ~ dnorm(0,0.444) # slope sst effect
    beta.bathy ~ dnorm(0,0.444)	# slope bathy effect
    beta.eff ~ dnorm(0,0.444)	# slope survey effort effect
    
    for (i in 1:nsite){
    
    # latent ecollogical model

    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- lpsi[i]
    lpsi[i] <- alpha.psi + beta.sst * SST[i] + beta.bathy * BATHY[i]
    
    # observation model 

    mu.p[i] <- z[i] * p[i] 
    logit(p[i]) <- lp[i] 
    lp[i] <- alpha.p + beta.eff * EFF[i]
    y[i] ~ dbern(mu.p[i])
    
    } #i
  }
    ", fill = TRUE)
sink()

#############################
# Data & run
#############################

load("AtSea_SV.rdata")

# see loaded data 

str(win.data)

# see initial values 

inits

# see saved parameters

params

# MCMC settings

ni <- 300000
nt <- 5
nb <- 150000
nc <- 3

# Call JAGS from R

out <- jags(win.data, params, "AtSea_SV.jags", inits = inits, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

AtSea_SV <- out

save(AtSea_SV, file="RES_AtSea_SV.rdata")

#############################
# Outputs
#############################

library(mcmcplots)

load("RES_AtSea_SV.rdata")

head(AtSea_SV$BUGSoutput$summary)

# see posterior distribution of parameters

denplot(AtSea_SV, c("alpha.psi","beta.bathy", "beta.sst"))

# see parameters values for saved iterations

traplot(AtSea_SV,"beta.bathy")
traplot(AtSea_SV,"beta.sst")
traplot(AtSea_SV,"beta.eff")