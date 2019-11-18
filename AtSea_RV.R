#############################
## At-Sea RV occupancy model  
#############################


############################# 
# Specify model in BUGS language
#############################

library(R2jags)

sink("AtSea_RV.jags")
cat("
    model {

    # priors

    alpha.psi ~ dnorm(0,0.444) # occupancy intercept
    alpha.p ~ dnorm(0,0.444) # detection intercept 
    beta.sst ~ dnorm(0,0.444) # slope sst effect
    beta.bathy ~ dnorm(0,0.444)	# slope bathy effect
    beta.eff ~ dnorm(0,0.444)	# slope survey effort effect
    beta.occ2 ~ dnorm(0,0.444) # occasion effect
    beta.occ3 ~ dnorm(0,0.444) # occasion effect
    beta.occ4 ~ dnorm(0,0.444) # occasion effect	
    
    for (i in 1:nsite){

      z[i] ~ dbern(psi[i])
      logit(psi[i]) <- lpsi[i]
      lpsi[i] <- alpha.psi + beta.sst * SST[i] + beta.bathy * BATHY[i]
      
      for (j in 1:nrep){

        mu.p[i,j] <- z[i] * p[i,j] 
        logit(p[i,j]) <- lp[i,j] 
        lp[i,j] <- alpha.p + beta.eff * EFF[i,j] + beta.occ2 * equals(j,2) + beta.occ3 * equals(j,3) + beta.occ4 * equals(j,4) 
        y[i,j] ~ dbern(mu.p[i,j])

      } #j
    } #i
    
}
    ", fill = TRUE)
sink()

#############################
# Data & run
#############################

load("AtSea_RV.rdata")

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

out <- jags(win.data, params, "AtSea_RV.jags", inits = inits, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

AtSea_RV <- out

save(AtSea_RV, file="RES_AtSea_RV.rdata")

#############################
# Outputs
#############################

library(mcmcplots)

load("RES_AtSea_RV.rdata")

head(AtSea_RV$BUGSoutput$summary)

# see posterior distribution of parameters

denplot(AtSea_RV, c("alpha.psi","beta.bathy", "beta.sst"))

# see parameters values for saved iterations

traplot(AtSea_RV,"beta.bathy")
traplot(AtSea_RV,"beta.sst")
traplot(AtSea_RV,"beta.eff")