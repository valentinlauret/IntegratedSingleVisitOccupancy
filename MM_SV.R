#############################
## Multi-method SV occupancy model  
#############################


############################# 
# Specify model in BUGS language
#############################

library(R2jags)

sink("MM_SV.jags")
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
    
    for (i in 1:nsite){

      z[i] ~ dbern(psi[i])
      logit(psi[i]) <- lpsi[i]
      lpsi[i] <- alpha.psi + beta.sst * SST[i] + beta.bathy * BATHY[i]
    
      # p_samm
      logit(p_samm[i]) <- lp_samm[i] 
      lp_samm[i] <- alpha.p_samm + beta.eff.samm * eff.samm[i] 
      
      # p_samm
      logit(p_gdegem[i]) <- lp_gdegem[i] 
      lp_gdegem[i] <- alpha.p_gdegem + beta.eff.gd * eff.gd[i] 
      
      
      mu.p[i,1] <- 1 - z[i] + z[i] * (1 - ind.samm[i] * p_samm[i]) * (1 - ind.gd[i] * p_gdegem[i]) # P(O=0) : no detection
      mu.p[i,2] <- z[i] * ind.samm[i] * p_samm[i] * (1 - ind.gd[i] * p_gdegem[i])                  # P(O=1) : detection by samm only
      mu.p[i,3] <- z[i] * ind.gd[i] * p_gdegem[i] * (1 - ind.samm[i] * p_samm[i])                 # P(O=2) : detection by gdegem only
      mu.p[i,4] <- z[i] * ind.gd[i] * ind.samm[i] * p_samm[i] * p_gdegem[i]                        # P(O=3) : detection by both samm and gdegem
    
    y[i] ~ dcat(mu.p[i,])
    
} #i
    
}
    ", fill = TRUE)
sink()

#############################
# Data & run
#############################

load("MM_SV.rdata")

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

out <- jags(win.data, params, "MM_SV.jags", inits = inits, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

MM_SV <- out

save(MM_SV, file="RES_MM_SV.rdata")

#############################
# Outputs
#############################

library(mcmcplots)

load("RES_MM_SV.rdata")

head(MM_SV$BUGSoutput$summary)

# see posterior distribution of parameters

denplot(MM_SV, c("alpha.psi","beta.bathy", "beta.sst"))

# see parameters values for saved iterations

traplot(MM_SV,"beta.bathy")
traplot(MM_SV,"beta.sst")
traplot(MM_SV,"beta.eff")