#############################
## Multi-method RV occupancy model  
#############################


############################# 
# Specify model in BUGS language
#############################

library(R2jags)

sink("MM_RV.jags")
cat("
    model {

# Note that SAMM is the aerial monitoring program, GDEGeM is the at-sea monitoring program
    # priors
    
    alpha.psi ~ dnorm(0,0.444) # occupancy intercept
    alpha.p_samm ~ dnorm(0,0.444) # detection samm intercept 
    alpha.p_gdegem ~ dnorm(0,0.444) # detection gdegem intercept 
    
    beta.sst ~ dnorm(0,0.444) # slope sst effect
    beta.bathy ~ dnorm(0,0.444)	# slope bathy effect
    beta.eff.samm ~ dnorm(0,0.444)	# slope SAMM survey effort effect
    beta.eff.gd ~ dnorm(0,0.444)	# slope GDEGeM survey effort effect
    beta.occ2 ~ dnorm(0,0.444) # occasion effect
    beta.occ3 ~ dnorm(0,0.444) # occasion effect
    beta.occ4 ~ dnorm(0,0.444) # occasion effect	
    
    for (i in 1:nsite){
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- lpsi[i]
    lpsi[i] <- alpha.psi + beta.sst * SST[i] + beta.bathy * BATHY[i]
    
    for (j in 1:nrep){
    
    # p_samm
    logit(p_samm[i,j]) <- lp_samm[i,j] 
    lp_samm[i,j] <- alpha.p_samm + beta.eff.samm * eff.samm[i,j] + beta.occ2 * equals(j,2) + beta.occ3 * equals(j,3) + beta.occ4 * equals(j,4) 
    
    # p_samm
    logit(p_gdegem[i,j]) <- lp_gdegem[i,j] 
    lp_gdegem[i,j] <- alpha.p_gdegem + beta.eff.gd * eff.gd[i,j] + beta.occ2 * equals(j,2) + beta.occ3 * equals(j,3) + beta.occ4 * equals(j,4) 
    
    mu.p[i,j,1] <- 1 - z[i] + z[i] * (1 - ind.samm[i,j] * p_samm[i,j]) * (1 - ind.gd[i,j] * p_gdegem[i,j]) # P(O=0) : no detection
    mu.p[i,j,2] <- z[i] * ind.samm[i,j] * p_samm[i,j] * (1 - ind.gd[i,j] * p_gdegem[i,j])                  # P(O=1) : detection by samm only
    mu.p[i,j,3] <- z[i] * ind.gd[i,j] * p_gdegem[i,j] * (1 - ind.samm[i,j] * p_samm[i,j])                 # P(O=2) : detection by gdegem only
    mu.p[i,j,4] <- z[i] * ind.gd[i,j] * ind.samm[i,j] * p_samm[i,j] * p_gdegem[i,j]                        # P(O=3) : detection by both samm and gdegem
    
    y[i,j] ~ dcat(mu.p[i,j,])
    } #j
    } #i
    }
    ", fill = TRUE)
sink()

#############################
# Data & run
#############################

load("MM_RV.rdata")

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

out <- jags(win.data, params, "MM_RV.jags", inits = inits, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

MM_RV <- out

save(MM_RV, file="RES_MM_RV.rdata")

#############################
# Outputs
#############################

library(mcmcplots)

load("RES_MM_RV.rdata")

head(MM_RV$BUGSoutput$summary)

# see posterior distribution of parameters

denplot(MM_RV, c("alpha.psi","beta.bathy", "beta.sst"))

# see parameters values for saved iterations

traplot(MM_RV,"beta.bathy")
traplot(MM_RV,"beta.sst")
traplot(MM_RV,"beta.eff")