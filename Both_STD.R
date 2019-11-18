# Modele d'occupancy GDEGeM ! 

# Specify model in BUGS language
library(R2jags)
sink("Both_std.jags")
cat("
    model {
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

#Process data
setwd("~/Documents/Google Drive/These/Work/Occupancy/Occupancy Standard/")
setwd("/home/ogimenez/valentin") # For the cluster

# load gdegem data and effort
load("GDEGeM_ready2run.rdata")
y.gd <- y
# load samm data and effort
load("SAMM_ready2run.rdata")
y.sa <- y

## effort for SAMM
eff.samm <- samm.effort[,2:5]
ind.samm <- samm.effort[,11:14]
## effort for GDEGEM
eff.gd <- cbind(as.numeric(gd.effort$occ1),as.numeric(gd.effort$occ2),as.numeric(gd.effort$occ3),as.numeric(gd.effort$occ4))
ind.gd <- cbind(as.numeric(gd.effort$occ1.ind),as.numeric(gd.effort$occ2.ind),as.numeric(gd.effort$occ3.ind),as.numeric(gd.effort$occ4.ind))
colnames(ind.gd) <- c("occ1.ind","occ2.ind","occ3.ind","occ4.ind")
# Create the tibble for combined detections 
y <-  data.frame(obs = y.sa$obs, occ1 = y.sa$occ1, occ2 = y.sa$occ2, occ3 = y.sa$occ3, occ4 =y.sa$occ4, bathy.sc = as.numeric(y.sa$bathy.sc), sst.sc = as.numeric(y.sa$sst.sc))
# we now have to add gdegem detections

y$obs[which(y.gd$obs==1)] <- 1

y$occ1[which(ind.gd[,1]==1)] <- ifelse(is.na(y$occ1[which(ind.gd[,1]==1)])==T,0,y$occ1[which(ind.gd[,1]==1)])
y$occ1[which(y.gd$occ1==1)] <- ifelse(y$occ1[which(y.gd$occ1==1)]==1,3,2)

y$occ2[which(ind.gd[,2]==1)] <- ifelse(is.na(y$occ2[which(ind.gd[,2]==1)])==T,0,y$occ2[which(ind.gd[,2]==1)])
y$occ2[which(y.gd$occ2==1)] <- ifelse(y$occ2[which(y.gd$occ2==1)]==1,3,2)

y$occ3[which(ind.gd[,3]==1)] <- ifelse(is.na(y$occ3[which(ind.gd[,3]==1)])==T,0,y$occ3[which(ind.gd[,3]==1)])
y$occ3[which(y.gd$occ3==1)] <- ifelse(y$occ3[which(y.gd$occ3==1)]==1,3,2)

y$occ4[which(ind.gd[,4]==1)] <- ifelse(is.na(y$occ4[which(ind.gd[,4]==1)])==T,0,y$occ4[which(ind.gd[,4]==1)])
y$occ4[which(y.gd$occ4==1)] <- ifelse(y$occ4[which(y.gd$occ4==1)]==1,3,2)

# Process data, part 2, by OG
dim(y)
str(y)
surveyed <- ind.gd + ind.samm # whether a site was surveyed or not
mask <- apply(surveyed,1,sum) # sites that were not surveyed have 0s
y <- y[mask!=0,] # filter out sites that were not surveyed at all throughout the study period
eff.samm <- eff.samm[mask!=0,] # get rid of effort for sites that were never surveyed
eff.gd <- eff.gd[mask!=0,] # get rid of effort for sites that were never surveyed
ind.samm <- ind.samm[mask!=0,]
ind.gd <- ind.gd[mask!=0,]
surveyed <- surveyed[mask!=0,] # update the survey indicator
dim(y)
det <- cbind(y$occ1, y$occ2, y$occ3, y$occ4) # detection and non-detections
det[surveyed==0] <- NA # for occasions with no effort, use NA

# Build dataset
win.data <- list(
  y = det+1, 
  eff.gd = eff.gd,
  eff.samm = eff.samm,
  ind.gd = ind.gd,
  ind.samm = ind.samm,
  BATHY = as.numeric(y$bathy.sc),
  SST = as.numeric(y$sst.sc), 
  nsite = nrow(y),
  nrep = 4)

# Initial values
# Initial values
zst <- ifelse(y$obs!=0,1,0)
inits <- function() {list(
  z = zst, 
  alpha.psi = rnorm(1, 0, 1.5),
  alpha.p_samm = rnorm(1, 0, 1.5),
  alpha.p_gdegem = rnorm(1, 0, 1.5),
  beta.sst = rnorm(1, 0, 1.5),
  beta.bathy = rnorm(1, 0, 1.5),
  beta.eff.samm = rnorm(1, 0, 1.5),
  beta.eff.gd = rnorm(1, 0, 1.5),
  beta.occ2 = rnorm(1, 0, 1.5),
  beta.occ3 = rnorm(1, 0, 1.5),
  beta.occ4 = rnorm(1, 0, 1.5))}

# Parameters monitored
params <- c("alpha.psi",
            "alpha.p_samm",
            "alpha.p_gdegem",
            "beta.eff.gd",
            "beta.eff.samm",
            "beta.occ2",
            "beta.occ3",
            "beta.occ4",
            "beta.bathy",
            "beta.sst","z") 

# MCMC settings
ni <- 20000
nt <- 5
nb <- 1000
nc <- 2

# Call JAGS from R
#ptm <- proc.time()
out <- jags(win.data, params, "Both_std.jags", inits = inits, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
#x <- proc.time() -  ptm

Both_std <-  out

save(Both_std, file ="RES_Both_std.rdata")

#github
save(win.data,inits,params, file = "MM_RV.RData")
