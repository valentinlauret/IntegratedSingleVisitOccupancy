# all simulations 
# CONF1

# pqckages 

library(R2jags)
# load data

# load("~/Google Drive/These/Work/Occupancy/Simulation/function_sim.Rdata")
# load("~/Google Drive/These/Work/Occupancy/Simulation/loop.Rdata")
setwd("/home/vlauret/")
getwd()
# load functionS
sim_state_cov <- function(nb_sites=nsites,init_occ = covparam, covenv = covenv){
  
  # define various quantities
  R = nb_sites # number of sites
  
  l.psi <- init_occ[1] + init_occ[2] * covenv + init_occ[3] * covenv * covenv
  psi1 <-  1 / ( 1 + exp( -l.psi ))# quadratic relation of the covariate with initial occupancy
  
  # pre-allocate memory
  site <- 1:R # Sites
  
  z <- array(dim = c(R)) # Expected and realized occurrence
  
  # define state process
  # first year/season
  z <- rbinom(R, 1, psi1) # Initial occupancy state
  
  
  return(z)
}
sim_obs_cov <- function(z=z, nb_occ = nocc,Seff, obsparam = obsparam){
  
  # define various quantities
  R = length(z) # number of sites
  K = nb_occ # number of occasions
  
  y <- array(NA, dim = c(R,K))
  p <- array(NA, dim = c(R,K))
  
  # define observation process
  
  for(i in 1:R){
    for(k in 1:K){
      
      l.p <- obsparam[1] + obsparam[2] * Seff[i,k] + obsparam[3] * Seff[i,k] * Seff[i,k]# quadratic relation of the covariate with initial occupancy
      
      p[i,k] <- ( 1 / ( 1 + exp(- l.p )))
      
      prob <- z[i] * p[i,k]
      y[i,k] <- rbinom(1,1, prob= prob)
    } #k
  } #i
  
  
  # format data
  return(y)
}
loop <- function(nsites, nocc, nsim, covparam, obsparam1,obsparam2, Seff, Seff2, Covenv){
  
  
  #### allocate memories for posterior estimates ####
  
  p_i <- rep(NA, nsim)
  psi_i <- rep(NA, nsim)

  bcov_i <- rep(NA, nsim)
  bseff_i <- rep(NA, nsim)

  p_ii <- rep(NA, nsim)
  psi_ii <- rep(NA, nsim)

  bcov_ii <- rep(NA, nsim)
  bseff_ii<- rep(NA, nsim)

  p1_iii <- rep(NA, nsim)
  p2_iii <- rep(NA, nsim)
  psi_iii <-rep(NA, nsim)

  bcov_iii <- rep(NA, nsim)
  bseff1_iii<- rep(NA, nsim)
  bseff2_iii<- rep(NA, nsim)

  p1_iv <- rep(NA, nsim)
  p2_iv <- rep(NA, nsim)
  psi_iv <-rep(NA, nsim)

  bcov_iv <- rep(NA, nsim)
  bseff1_iv<- rep(NA, nsim)
  bseff2_iv<- rep(NA, nsim)

  # bis high p 
  p_ibis <- rep(NA, nsim)
  psi_ibis <- rep(NA, nsim)

  bcov_ibis <- rep(NA, nsim)
  bseff_ibis <- rep(NA, nsim)

  p_iibis <- rep(NA, nsim)
  psi_iibis <- rep(NA, nsim)

  bcov_iibis <- rep(NA, nsim)
  bseff_iibis<- rep(NA, nsim)

  # allocate memory to saved parameters
  
  rbpsi_i    <- rep(NA, nsim)
  rmsepsi_i   <- rep(NA, nsim)
  rbcov_i     <- rep(NA, nsim)
  rmsecov_i <- rep(NA, nsim)
  rbp_i     <- rep(NA, nsim)
  rmsep_i    <- rep(NA, nsim)

  rbseff_i    <- rep(NA, nsim)
  rmseseff_i  <- rep(NA, nsim)
  # ii<- rep(NA, nsim)
  rbpsi_ii    <- rep(NA, nsim)
  rmsepsi_ii  <- rep(NA, nsim)
  
  rbcov_ii    <- rep(NA, nsim)
  rmsecov_ii  <- rep(NA, nsim)
 
  rbp_ii      <- rep(NA, nsim)
  rmsep_ii    <- rep(NA, nsim)

  rbseff_ii   <- rep(NA, nsim)
  rmseseff_ii <- rep(NA, nsim)

  # iii   <- rep(NA, nsim)
  rbpsi_iii   <- rep(NA, nsim)
  rmsepsi_iii <- rep(NA, nsim)
  
  rbcov_iii   <- rep(NA, nsim)
  rmsecov_iii <- rep(NA, nsim)
  
  rbp1_iii    <- rep(NA, nsim)
  rmsep1_iii  <- rep(NA, nsim)
  rbseff1_iii <- rep(NA, nsim)
  rmseseff1_iii<- rep(NA, nsim)
  
  rbp2_iii    <- rep(NA, nsim)
  rmsep2_iii  <- rep(NA, nsim)
  rbseff2_iii <- rep(NA, nsim)
  rmseseff2_iii<- rep(NA, nsim)
                         # iv
  rbpsi_iv    <- rep(NA, nsim)
  rmsepsi_iv  <- rep(NA, nsim)
  rbcov_iv    <- rep(NA, nsim)
  rmsecov_iv  <- rep(NA, nsim)
  rbp1_iv     <- rep(NA, nsim)
  rmsep1_iv   <- rep(NA, nsim)
  rbseff1_iv  <- rep(NA, nsim)
  rmseseff1_iv<- rep(NA, nsim)
  rbp2_iv     <- rep(NA, nsim)
  rmsep2_iv   <- rep(NA, nsim)
  rbseff2_iv  <- rep(NA, nsim)
  rmseseff2_iv<- rep(NA, nsim)
  
  # bis highp    
  p_ibis <- rep(NA, nsim)
  psi_ibis <- rep(NA, nsim)

  bcov_ibis <- rep(NA, nsim)
  bseff_ibis <- rep(NA, nsim)

  p_iibis <- rep(NA, nsim)
  psi_iibis <- rep(NA, nsim)

  bcov_iibis <- rep(NA, nsim)
  bseff_iibis<- rep(NA, nsim)
  
  # parameters
  rbpsi_ibis    <- rep(NA, nsim)
  rmsepsi_ibis   <- rep(NA, nsim)
  rbcov_ibis     <- rep(NA, nsim)
  rmsecov_ibis <- rep(NA, nsim)
  rbp_ibis     <- rep(NA, nsim)
  rmsep_ibis    <- rep(NA, nsim)
  
  rbseff_ibis    <- rep(NA, nsim)
  rmseseff_ibis  <- rep(NA, nsim)
  # ii<- rep(NA, nsim)
  rbpsi_iibis    <- rep(NA, nsim)
  rmsepsi_iibis  <- rep(NA, nsim)
  
  rbcov_iibis    <- rep(NA, nsim)
  rmsecov_iibis  <- rep(NA, nsim)
  
  rbp_iibis      <- rep(NA, nsim)
  rmsep_iibis    <- rep(NA, nsim)
  
  rbseff_iibis   <- rep(NA, nsim)
  rmseseff_iibis <- rep(NA, nsim)

  # JAGS models 
  
  sink("Sim_occu_cov.jags")
  cat("
    model {
    
    # Specify priors
    
    alpha.psi ~ dnorm(0,0.444) # intercept occupancy
    alpha.p ~ dnorm(0,0.444) # intercept detection

    beta.cov ~ dnorm(0,0.444) # slope cov effect
    beta.qcov ~ dnorm(0,0.444)	# slope quadratic cov effect
    
    beta.seff ~ dnorm(0,0.444)	# slope Sampling effort effect
    beta.qseff ~ dnorm(0,0.444)	# slope quadratic sampling effort
    
    beta.occ2 ~ dnorm(0,0.444) # occasion effect
    beta.occ3 ~ dnorm(0,0.444) # occasion effect
    beta.occ4 ~ dnorm(0,0.444) # occasion effect
    
    
    ### Ecological states
    
    for (i in 1:nsite){
    
      z[i] ~ dbern(muZ[i])
      logit(muZ[i]) <- l.muZ[i]
      
      l.muZ[i] <- alpha.psi + beta.cov * covenv[i] + beta.qcov * covenv[i] * covenv[i]
    
    } # i
  
    # Observation model
   
    for (i in 1:nsite){
    
      for (k in 1:nocc){
    
          logit(p[i,k]) <- -l.p[i,k]
          l.p[i,k] <- alpha.p + beta.seff * Seff[i,k] + beta.qseff * Seff[i,k] * Seff[i,k] + beta.occ2 * equals(k,2) + beta.occ3 * equals(k,3) + beta.occ4 * equals(k,4)
          

    y[i,k] ~ dbern(z[i]*p[i,k])    
    
      } #k
    } #i
    
    
    }
    ",fill = TRUE)
  sink()
  
  sink("sim_IOM_cov.jags")
  cat("
    model {
    # priors

    alpha.psi ~ dnorm(0,0.444) # occupancy intercept
    alpha.p1 ~ dnorm(0,0.444) # detection samm intercept 
    alpha.p2 ~ dnorm(0,0.444) # detection gdegem intercept 
    
    # slope priors 
    
    beta.cov ~ dnorm(0,0.444) # slope cov effect
    beta.qcov ~ dnorm(0,0.444)	# slope quadratic cov effect
    
    beta.seff1 ~ dnorm(0,0.444)	# slope Sampling effort effect
    beta.qseff1 ~ dnorm(0,0.444)	# slope quadratic sampling effort
    
    beta.seff2 ~ dnorm(0,0.444)	# slope Sampling effort effect
    beta.qseff2 ~ dnorm(0,0.444)	# slope quadratic sampling effort
    
    beta.occ2 ~ dnorm(0,0.444) # occasion effect
    beta.occ3 ~ dnorm(0,0.444) # occasion effect
    beta.occ4 ~ dnorm(0,0.444) # occasion effect

    # NO COV
    
    # Ecological state 
    
    for (i in 1:nsite){
    
       z[i] ~ dbern(muZ[i])
      logit(muZ[i]) <- l.muZ[i]
      
        l.muZ[i] <- alpha.psi + beta.cov * covenv[i] + beta.qcov * covenv[i] * covenv[i]

      for (j in 1:nocc){

        # p1
        logit(p1[i,j]) <- lp1[i,j] 
        lp1[i,j] <- alpha.p1  + beta.seff1 * Seff1[i,j] + beta.seff1 * Seff1[i,j] * Seff1[i,j] + beta.occ2 * equals(j,2) + beta.occ3 * equals(j,3) + beta.occ4 * equals(j,4)
        
         # p2
        logit(p2[i,j]) <- lp2[i,j] 
        lp2[i,j] <- alpha.p2 + beta.seff2 * Seff2[i,j] + beta.qseff2 * Seff2[i,j] * Seff2[i,j] + beta.occ2 * equals(j,2) + beta.occ3 * equals(j,3) + beta.occ4 * equals(j,4)
          
        
        mu.p[i,j,1] <- 1 - z[i] + z[i] * (1 - p1[i,j]) * (1 - p2[i,j]) # P(O=0) : no detection
        mu.p[i,j,2] <- z[i] * p1[i,j] * (1 - p2[i,j])                  # P(O=1) : detection by samm only
        mu.p[i,j,3] <- z[i] * p2[i,j] * (1 - p1[i,j])                 # P(O=2) : detection by gdegem only
        mu.p[i,j,4] <- z[i] * p1[i,j] * p2[i,j]                        # P(O=3) : detection by both samm and gdegem

        y[i,j] ~ dcat(mu.p[i,j,])
    } #j
    } #i
    }
    ", fill = TRUE)
  sink()
  
  #### loop over nsim ####
  for( s in 1:nsim){
    
    # simulate datasets s
    
    z <- sim_state_cov(nb_sites = nsites, covenv= covenv, init_occ= covparam)
    y <- sim_obs_cov(z,nocc,Seff= Seff, obsparam= obsparam1)
    
    y1 <- y # the first detection dataset
    p1 <- obsparam1
    p2 <- obsparam2
    y2bis <- y2 <- sim_obs_cov(z,nocc,obsparam = obsparam2, Seff =Seff2) # another dataset with different p
    y2bis[y2>0] <- 2 
    y_int <- y1 + y2bis # 0 -> no detection, 1 -> detection by monitoring 1, 2-> detection by 2, 3 -> detection by both methods
    
    
    
    # perform analysis
    # Specify MCMC settings one for all
    ni <-  4000
    nb <- 1000
    nc <- 1
    
    #### i)  One monitoring RV low p ####
    
    
    # Bundle data
    
    win.data <- list(y = y1,nsite = dim(y)[1], nocc = dim(y)[2], Seff = Seff, covenv = covenv)
    
    # Initial values
    zst <- apply(y, 1, max)	# Observed occurrence as inits for z
    psi = covparam[1]    # sum(z)/nsites
    p = obsparam1[1]
    bcov <- covparam[2]
    bseff <- obsparam1[2]
    
    inits <- function(){ list(z = zst, alpha.psi = psi, alpha.p= p, beta.cov = bcov, beta.seff= bseff)}
    
    # Parameters monitored
    params <- c( "alpha.p","alpha.psi", "beta.cov", "beta.qcov","beta.seff" ,"beta.qseff") 
    
    
    # Call JAGS from R (BRT 3 min)
    out <- jags(win.data, inits, params, "Sim_occu_cov.jags", n.chains = nc, n.iter = ni, n.burnin = nb, working.directory = getwd())
    
    # denplot(out)
    # head(out$BUGSoutput$summary)
    
    # save results 
    
    # intercept
    p_i[s] <- out$BUGSoutput$mean$alpha.p
    psi_i[s] <- out$BUGSoutput$mean$alpha.psi

    # slope linear estimate (Q: do we have to save the quadratic estimates-)
    
    bcov_i[s] <-    out$BUGSoutput$mean$beta.cov
    bseff_i[s] <-   out$BUGSoutput$mean$beta.seff

    #### i bis )  One monitoring RV high p ####
    
    
    # Bundle data
    
    win.data <- list(y = y2,nsite = dim(y2)[1], nocc = dim(y2)[2], Seff = Seff2, covenv = covenv)
    
   
    # Initial values
    zst <- apply(y2, 1, max)	# Observed occurrence as inits for z
    psi = covparam[1]   # sum(z)/nsites
    p = obsparam2[1]
    bcov <- covparam[2]
    bseff <- obsparam2[2]
    
    inits <- function(){ list(z = zst, alpha.psi = psi, alpha.p= p, beta.cov = bcov, beta.seff= bseff)}
    
    # Parameters monitored
    params <- c( "alpha.p","alpha.psi", "beta.cov", "beta.qcov","beta.seff" ,"beta.qseff") 
    
    
    # Call JAGS from R (BRT 3 min)
    out <- jags(win.data, inits, params, "Sim_occu_cov.jags", n.chains = nc, n.iter = ni, n.burnin = nb, working.directory = getwd())
    
   # head(out$BUGSoutput$summary)
    
    # save results 
    
    # intercept
    p_ibis[s] <- out$BUGSoutput$mean$alpha.p
    psi_ibis[s] <- out$BUGSoutput$mean$alpha.psi

    # slope linear estimate (Q: do we have to save the quadratic estimates-)
    
    bcov_ibis[s] <- out$BUGSoutput$mean$beta.cov
    bseff_ibis[s] <- out$BUGSoutput$mean$beta.seff

    
     #### ii)  One monitoring SV ####
    
    # same model jags than i)
    
    # Bundle data
    
    ysv <-  matrix(apply(y,1,max))
    Seffsv <- matrix(apply(Seff,1,sum))
    
    win.data <- list(y = ysv,nsite = dim(ysv)[1], nocc = dim(ysv)[2], Seff = Seffsv, covenv = covenv)
    
    
    # Initial values
    zst <- apply(ysv, 1, max)	# Observed occurrence as inits for z
    psi = covparam[1]   # sum(z)/nsites
    p = obsparam1[1]
    bcov <- covparam[2]
    bseff <- obsparam1[2]
    
    inits <- function(){ list(z = zst, alpha.psi = psi, alpha.p= p, beta.cov = bcov, beta.seff= bseff)}
    
    
    # Parameters monitored
    params <- c("z", "alpha.p","alpha.psi", "beta.cov", "beta.seff","beta.qcov", "beta.qseff") 
    
    
    # Call JAGS from R (BRT 3 min)
   
    out <- jags(win.data, inits, params, "Sim_occu_cov.jags", n.chains = nc,  n.iter = ni,  n.burnin = nb, working.directory = getwd())
    
    # head(out$BUGSoutput$summary)
    
    p_ii[s] <- out$BUGSoutput$mean$alpha.p
    psi_ii[s] <- out$BUGSoutput$mean$alpha.psi

    # slope linear estimate (Q: do we have to save the quadratic estimates-)
    
    bcov_ii[s] <- out$BUGSoutput$mean$beta.cov
    bseff_ii[s] <- out$BUGSoutput$mean$beta.seff

    #### ii bis )  One monitoring SV high p ####
    
    # same model jags than i)
    
    # Bundle data
    
    y2sv <-  matrix(apply(y2,1,max))
    
    Seff2sv <- matrix(apply(Seff2,1,max))
    win.data <- list(y = y2sv,nsite = dim(y2sv)[1], nocc = dim(y2sv)[2], Seff = Seff2sv, covenv = covenv)
    
    
    # Initial values
    zst <- apply(y2sv, 1, max)	# Observed occurrence as inits for z
    psi = covparam[1]   # sum(z)/nsites
    p = obsparam2[1]
    bcov <- covparam[2]
    bseff <- obsparam2[2]
    
    inits <- function(){ list(z = zst, alpha.psi = psi, alpha.p= p, beta.cov = bcov, beta.seff= bseff)}
    
    
    
    # Parameters monitored
    params <- c("z", "alpha.p","alpha.psi", "beta.cov", "beta.seff") 
    
    
    # Call JAGS from R (BRT 3 min)
    out <- jags(win.data, inits, params, "Sim_occu_cov.jags", n.chains = nc,  n.iter = ni,  n.burnin = nb, working.directory = getwd())
    

    p_iibis[s] <- out$BUGSoutput$mean$alpha.p
    psi_iibis[s] <- out$BUGSoutput$mean$alpha.psi

    # slope linear estimate (Q: do we have to save the quadratic estimates-)
    
    bcov_iibis[s] <- out$BUGSoutput$mean$beta.cov
    bseff_iibis[s] <- out$BUGSoutput$mean$beta.seff

    
    #### iii)  Two monitoring RV ####
    
    # we work with y_int
    win.data <- list(y = y_int+1,nsite = dim(y_int)[1], nocc = dim(y_int)[2], covenv = covenv, Seff1 = Seff1, Seff2 = Seff2)
    

    # Initial values
    zst <- apply(y_int, 1, max)
    zst[zst>0] <- 1 # Observed occurrence as inits for z
    psi = covparam[1]   # sum(z)/nsites
    p1 = obsparam1[1]
    p2 = obsparam2[1]
    bcov <- covparam[2]
    bseff1 <- obsparam1[2]
    bseff2 <- obsparam2[2]
    
    inits <- function(){ list(z = zst, alpha.psi = psi, alpha.p1=p1, alpha.p2=p2,beta.cov = bcov, beta.seff1= bseff1, beta.seff2= bseff2)}
    
    # Parameters monitored
    params <- c("z", "alpha.p1","alpha.psi","alpha.p2","beta.cov","beta.qcov","beta.seff1","beta.seff2","beta.qseff1","beta.qseff2") 
    
    
    # Call JAGS from R (BRT 3 min)
    out <- jags(win.data, inits, params, "sim_IOM_cov.jags", n.chains = nc,  n.iter = ni,n.burnin = nb, working.directory = getwd())
    
    p1_iii[s] <-out$BUGSoutput$mean$alpha.p1
    p2_iii[s] <- out$BUGSoutput$mean$alpha.p2

    psi_iii[s] <- out$BUGSoutput$mean$alpha.psi

    # save results slope 
    
    bseff1_iii[s] <-out$BUGSoutput$mean$beta.seff1

    bseff2_iii[s] <- out$BUGSoutput$mean$beta.seff2

    bcov_iii[s] <- out$BUGSoutput$mean$beta.cov

    
    #### iv)  Two monitoring SV ####
    
    # same model jags than iii)
    
    # Bundle data
    
    y1sv <-  matrix(apply(y1,1,max))
    y2sv <-  matrix(apply(y2,1,max))
    y2sv[y2sv>0] <- 2
    
    ysv_int <- y1sv + y2sv
    
    Seff1sv <-  matrix(apply(Seff1,1,sum))
    Seff2sv <-  matrix(apply(Seff2,1,sum))
    
    win.data <- list(y = ysv_int+1,nsite = dim(ysv_int)[1], nocc = dim(ysv_int)[2], Seff1= Seff1sv, Seff2 = Seff2sv, covenv = covenv)
    
    # Initial values
    zst <- apply(ysv_int, 1, max)
    zst[zst>0] <- 1 # Observed occurrence as inits for z
    psi = covparam[1]   # sum(z)/nsites
    p1 = obsparam1[1]
    p2 = obsparam2[1]
    bcov <- covparam[2]
    bseff1 <- obsparam1[2]
    bseff2 <- obsparam2[2]
    
    inits <- function(){ list(z = zst, alpha.psi = psi, alpha.p1=p1, alpha.p2=p2,beta.cov = bcov, beta.seff1= bseff1, beta.seff2= bseff2)}
    
    # Parameters monitored
    params <- c("z", "alpha.p1","alpha.psi","alpha.p2","beta.cov","beta.qcov","beta.seff1","beta.seff2","beta.qseff1","beta.qseff2") 
    
    
    # Call JAGS from R (BRT 3 min)
    out <- jags(win.data, inits, params, "sim_IOM_cov.jags", n.chains = nc,  n.iter = ni, n.burnin = nb, working.directory = getwd())
    
    p1_iv[s] <-out$BUGSoutput$mean$alpha.p1
    p2_iv[s] <- out$BUGSoutput$mean$alpha.p2

    psi_iv[s] <-out$BUGSoutput$mean$alpha.psi

    # save results slope 
    
    bseff1_iv[s] <-out$BUGSoutput$mean$beta.seff1
    bseff2_iv[s] <- out$BUGSoutput$mean$beta.seff2
    bcov_iv[s] <- out$BUGSoutput$mean$beta.cov

    # calcul and save RMSE and RB
    
    # i 
    rbpsi_i[s]        <- (psi_i[s]-covparam[1])/covparam[1]   # RB alpha.psi
    rmsepsi_i[s]      <- (psi_i[s] - covparam[1])^2 # RMSE alpha.psi
    
    rbcov_i[s]        <- (bcov_i[s]-covparam[2])/covparam[2]   # RB cov
    rmsecov_i[s]      <- (bcov_i[s] - covparam[2])^2 # RMSE cov
    
    rbp_i[s]          <- (p_i[s]-obsparam1[1])/obsparam1[1]       # RB alpha p 
    rmsep_i[s]        <- (p_i[s] - obsparam1[1])^2 # RMSE alpha p
    
    rbseff_i[s]       <- (bseff_i[s]-obsparam1[2])/obsparam1[2]   # RB obs
    rmseseff_i[s]     <- (bseff_i[s] - obsparam1[2])^2 # RMSE obs
    
    
    # ii
    rbpsi_ii[s]       <- (psi_ii[s]-covparam[1])/covparam[1]   # RB alpha.psi
    rmsepsi_ii[s]     <- (psi_ii[s] - covparam[1])^2 # RMSE alpha.psi
    
    rbcov_ii[s]       <- (bcov_ii[s]-covparam[2])/covparam[2]  # RB cov
    rmsecov_ii[s]     <- (bcov_ii[s] - covparam[2])^2 # RMSE cov
    
    rbp_ii[s]         <- (p_ii[s]-obsparam1[1])/obsparam1[1]      # RB alpha p 
    rmsep_ii[s]       <- (p_ii[s] - obsparam1[1])^2 # RMSE alpha p
    
    rbseff_ii[s]      <- (bseff_ii[s]-obsparam1[2])/obsparam1[2]   # RB obs
    rmseseff_ii[s]    <- (bseff_ii[s] - obsparam1[2])^2 # RMSE obs
    
    
    # iii   
    rbpsi_iii[s]      <- (psi_iii[s]-covparam[1])/covparam[1]  # RB alpha.psi
    rmsepsi_iii[s]    <- (psi_iii[s] - covparam[1])^2 # RMSE alpha.psi
    
    rbcov_iii[s]      <- (bcov_iii[s]-covparam[2])/covparam[2]   # RB cov
    rmsecov_iii[s]    <- (bcov_iii[s] - covparam[2])^2 # RMSE cov
    
    rbp1_iii[s]       <- (p1_iii[s]-obsparam1[1])/obsparam1[1]        # RB alpha p 
    rmsep1_iii[s]     <- (p1_iii[s] - obsparam1[1])^2 # RMSE alpha p
    
    rbseff1_iii[s]     <- (bseff1_iii[s]-obsparam1[2])/obsparam1[2]  # RB obs
    rmseseff1_iii[s]   <- (bseff1_iii[s] - obsparam1[2])^2 # RMSE obs
    
    rbp2_iii[s]        <- (p2_iii[s]-obsparam2[1])/obsparam2[1]       # RB alpha p 
    rmsep2_iii[s]      <- (p2_iii[s] - obsparam2[1])^2 # RMSE alpha p
    
    rbseff2_iii[s]     <- (bseff2_iii[s]-obsparam2[2])/obsparam2[2]  # RB obs
    rmseseff2_iii[s]   <- (bseff2_iii[s] - obsparam2[2])^2 # RMSE obs
    
    
    # iv
    rbpsi_iv[s]        <- (psi_iv[s]-covparam[1])/covparam[1]  # RB alpha.psi
    rmsepsi_iv[s]      <- (psi_iv[s] - covparam[1])^2 # RMSE alpha.psi
    
    rbcov_iv[s]        <- (bcov_iv[s]-covparam[2])/covparam[2]   # RB cov
    rmsecov_iv[s]      <- (bcov_iv[s] - covparam[2])^2 # RMSE cov
    
    rbp1_iv[s]         <- (p1_iv[s]-obsparam1[1])/obsparam1[1]        # RB alpha p 
    rmsep1_iv[s]       <- (p1_iv[s] - obsparam1[1])^2 # RMSE alpha p
    
    rbseff1_iv[s]      <- (bseff1_iv[s]-obsparam1[2])/obsparam1[2]   # RB obs
    rmseseff1_iv[s]    <- (bseff1_iv[s] - obsparam1[2])^2 # RMSE obs
    
    rbp2_iv[s]         <- (p2_iv[s]-obsparam2[1])/obsparam2[1]       # RB alpha p 
    rmsep2_iv[s]       <- (p2_iv[s] - obsparam2[2])^2 # RMSE alpha p
    
    rbseff2_iv[s]      <- (bseff2_iv[s]-obsparam2[2])/obsparam2[2]   # RB obs
    rmseseff2_iv[s]    <- (bseff2_iv[s] - obsparam2[2])^2 # RMSE obs
    
    # ibis  
    rbpsi_ibis[s]        <- (psi_ibis[s]-covparam[1])/covparam[1]   # RB alpha.psi
    rmsepsi_ibis[s]      <- (psi_ibis[s] - covparam[1])^2 # RMSE alpha.psi
    
    rbcov_ibis[s]        <- (bcov_ibis[s]-covparam[2])/covparam[2]   # RB cov
    rmsecov_ibis[s]      <- (bcov_ibis[s] - covparam[2])^2 # RMSE cov
    
    rbp_ibis[s]          <- (p_ibis[s]-obsparam2[1])/obsparam2[1]       # RB alpha p 
    rmsep_ibis[s]        <- (p_ibis[s] - obsparam2[1])^2 # RMSE alpha p
    
    rbseff_ibis[s]       <- (bseff_ibis[s]-obsparam2[2])/obsparam2[2]   # RB obs
    rmseseff_ibis[s]     <- (bseff_ibis[s] - obsparam2[2])^2 # RMSE obs
    
    
    # ii bis
    rbpsi_iibis[s]       <- (psi_iibis[s]-covparam[1])/covparam[1]   # RB alpha.psi
    rmsepsi_iibis[s]     <- (psi_iibis[s] - covparam[1])^2 # RMSE alpha.psi
    
    rbcov_iibis[s]       <- (bcov_iibis[s]-covparam[2])/covparam[2]  # RB cov
    rmsecov_iibis[s]     <- (bcov_iibis[s] - covparam[2])^2 # RMSE cov
    
    rbp_iibis[s]         <- (p_iibis[s]-obsparam2[1])/obsparam2[1]      # RB alpha p 
    rmsep_iibis[s]       <- (p_iibis[s] - obsparam2[1])^2 # RMSE alpha p
    
    rbseff_iibis[s]      <- (bseff_iibis[s]-obsparam2[2])/obsparam2[2]   # RB obs
    rmseseff_iibis[s]    <- (bseff_iibis[s] - obsparam2[2])^2 # RMSE obs
    
    
  } # end for 
  
  

  para <- data.frame("nsites" = nsites,"nocc" =  nocc, "nsim" = nsim, "covparam"= covparam, "obsparam1"= obsparam1,"obsparam2" = obsparam2)
  
  return(list(para, # i 
              rbpsi_i,rmsepsi_i , rbcov_i ,
              rmsecov_i,
              rbp_i ,
              rmsep_i ,
              rbseff_i , 
              rmseseff_i ,
              # ii
              rbpsi_ii,
              rmsepsi_ii,
              rbcov_ii,
              rmsecov_ii,
              rbp_ii,
              rmsep_ii,
              rbseff_ii,
              rmseseff_ii ,
              # iii
              rbpsi_iii,
              rmsepsi_iii ,
              rbcov_iii,
              rmsecov_iii ,
              rbp1_iii,
              rmsep1_iii,
              rbseff1_iii ,
              rmseseff1_iii,
              rbp2_iii,
              rmsep2_iii,
              rbseff2_iii ,
              rmseseff2_iii,
              # iv
              rbpsi_iv ,
              rmsepsi_iv,
              rbcov_iv,
              rmsecov_iv,
              rbp1_iv,
              rmsep1_iv,
              rbseff1_iv,
              rmseseff1_iv,
              rbp2_iv ,
              rmsep2_iv,
              rbseff2_iv,
              rmseseff2_iv, 
              # ibis
              rbpsi_ibis,rmsepsi_ibis , rbcov_ibis ,
              rmsecov_ibis,
              rbp_ibis ,
              rmsep_ibis ,
              rbseff_ibis , 
              rmseseff_ibis ,
              # iibis
              rbpsi_iibis,
              rmsepsi_iibis,
              rbcov_iibis,
              rmsecov_iibis,
              rbp_iibis,
              rmsep_iibis,
              rbseff_iibis,
              rmseseff_iibis))
  
} # endfunction

#### ####
#load("/home/vlauret/SimulationIO/loop.Rdata")

nsim = 100

nsites = 100
nocc = 4

covparam.3 <- c(0.5,0.2,0) # psi = 0.3
covparam.1 <- c(-1.9,0.2,0) # psi = 0.1

covparam <- list(covparam.1 = covparam.1,
                 covparam.3 = covparam.3)


obsparam1 <- c(-1.5,0.25,0) # p = 0.15

obsparam2 <- c(0.5,0.26,0) # p = 0.50

Seff <- Seff1 <- data.frame(cbind(rnorm(nsites),rnorm(nsites),rnorm(nsites),rnorm(nsites)))
Seff2 <- data.frame(cbind(rnorm(nsites),rnorm(nsites),rnorm(nsites),rnorm(nsites)))

Covenv <- covenv <- runif(nsites,0,1)

# all simul for loop 
res <- list()
i <- 0
t <- proc.time()
for(c in 1: length(covparam)){  #i <- i+1
  

    i <- i+1 
    
    print(paste("cov",c))
    res[[i]] <- loop(nsites, nocc, nsim, covparam[[c]], obsparam1,obsparam2, Seff, Seff2, Covenv)
    
    
} # end c

run.time <- proc.time() - t

  names(res[[1]]) <-  names(res[[2]]) <-  c("para", "rbpsi_i","rmsepsi_i" , "rbcov_i" ,
  "rmsecov_i",
  "rbp_i" ,
  "rmsep_i" ,
  "rbseff_i" , 
  "rmseseff_i" ,
  "rbpsi_ii",
  "rmsepsi_ii",
  "rbcov_ii",
  "rmsecov_ii",
  "rbp_ii",
  "rmsep_ii",
  "rbseff_ii",
  "rmseseff_ii" ,
  "rbpsi_iii",
  "rmsepsi_iii" ,
  "rbcov_iii",
  "rmsecov_iii" ,
  "rbp1_iii",
  "rmsep1_iii",
  "rbseff1_iii" ,
  "rmseseff1_iii",
  "rbp2_iii",
  "rmsep2_iii",
  "rbseff2_iii" ,
  "rmseseff2_iii",
  "rbpsi_iv",
  "rmsepsi_iv",
  "rbcov_iv",
  "rmsecov_iv",
  "rbp1_iv",
  "rmsep1_iv",
  "rbseff1_iv",
  "rmseseff1_iv",
  "rbp2_iv" ,
  "rmsep2_iv",
  "rbseff2_iv",
  "rmseseff2_iv", 
  # bis
  "rbpsi_ibis","rmsepsi_ibis" , "rbcov_ibis" ,
  "rmsecov_ibis",
  "rbp_ibis" ,
  "rmsep_ibis" ,
  "rbseff_ibis" , 
  "rmseseff_ibis" ,
  "rbpsi_iibis",
  "rmsepsi_iibis",
  "rbcov_iibis",
  "rmsecov_iibis",
  "rbp_iibis",
  "rmsep_iibis",
  "rbseff_iibis",
  "rmseseff_iibis")
 
  names(res) <- c("occ1", "occ3")
  
save(res,run.time,file = "simIOM.rdata")

# ----- plot ---- 
#  library(tidyverse)
#  library(cowplot)
#  
#  nsim = 100
#  
#  # make a tibble to store values 
#  extract.RMSE <- function(res, para){
#    m <- c()
#    for(i in 1:length(res)){
#      m[i] <- sqrt(mean(res[[i]][[para]]))
#    } # i
#    return(m)
#  }
#  
#  extract.RB <- function(res, para){
#    m <- c()
#    for(i in 1:length(res)){
#      m[i] <- mean(res[[i]][[para]])
#    } # i
#    return(m)
#  }
#  
#
#  res1 <- res$occ1
#  res2 <- res$occ3
#  
#  rocc <- tibble(model = rep(c("cov_i","cov_ii","cov_ibis","cov_iibis","cov_iii","cov_iv"), 2),
#               rb.occ = c(mean(res1$rbcov_i),
#                          mean(res1$rbcov_ii),
#                          mean(res1$rbcov_ibis),
#                          mean(res1$rbcov_iibis),
#                          mean(res1$rbcov_iii),
#                          mean(res1$rbcov_iv),
#                          mean(res2$rbcov_i),
#                          mean(res2$rbcov_ii),
#                          mean(res2$rbcov_ibis),
#                          mean(res2$rbcov_iibis),
#                          mean(res2$rbcov_iii),
#                          mean(res2$rbcov_iv)), 
#               rmse.occ = c(sqrt(mean(res1$rmsecov_i)),
#                          sqrt(mean(res1$rmsecov_ii)),
#                          sqrt(mean(res1$rmsecov_ibis)),
#                          sqrt(mean(res1$rmsecov_iibis)),
#                          sqrt(mean(res1$rmsecov_iii)),
#                          sqrt(mean(res1$rmsecov_iv)),
#                          sqrt(mean(res2$rmsecov_i)),
#                          sqrt(mean(res2$rmsecov_ii)),
#                          sqrt(mean(res2$rmsecov_ibis)),
#                          sqrt(mean(res2$rmsecov_iibis)),
#                          sqrt(mean(res2$rmsecov_iii)),
#                          sqrt(mean(res2$rmsecov_iv))))
#    
#    rocc <- rocc %>% mutate(type = rep(c("RV","SV"),length(model)/2),
#           occ = c(rep("Ψ = 0.1",length(model)/2),rep("Ψ = 0.3", length(model)/2)),
#           name = rep(c("SOM p1","SOM p1","SOM p2","SOM p2","IOM","IOM"),2))
#  
#  
#  pocc1 <- rocc %>%
#    ggplot() +
#    aes(x = name, y = rmse.occ, fill = type) + 
#    geom_col(position = "dodge", width=.6) + 
#    labs(x = 'Model',
#         y = 'RMSE in occupancy covariate effect size') + 
#    scale_fill_manual(name = NULL, values=c("#442b48","#32a287"))+ 
#    facet_wrap(~ occ,  strip.position = "top", labeller = label_wrap_gen(multi_line = FALSE)) + 
#    theme_bw(base_size = 14)
#  pocc1
#  
#  pocc2 <- rocc %>%
#    ggplot() +
#    aes(x = name, y = rb.occ, fill = type) + 
#    geom_col(position = "dodge", width=.6) + 
#    labs(x = 'Model',
#         y = 'RB in occupancy covariate effect size') + 
#    scale_fill_manual(name = NULL, values=c("#442b48","#32a287"))+ 
#    facet_wrap(~ occ,  strip.position = "top", labeller = label_wrap_gen(multi_line = FALSE)) + 
#    theme_bw(base_size = 14)
#  
#  # seff
#  rseff <- tibble(model = rep(c("seff_i","seff_ii","seff_ibis","rseff_iibis","seff1_iii","seff12_iii","seff1_iv","seff2_iv"), 2),
#                 rb.seff = c(mean(res1$rbseff_i),
#                            mean(res1$rbseff_ii),
#                            mean(res1$rbseff_ibis),
#                            mean(res1$rbseff_iibis),
#                            mean(res1$rbseff1_iii),
#                            mean(res1$rbseff2_iii),
#                            mean(res1$rbseff1_iv),
#                            mean(res1$rbseff2_iv),
#                            mean(res2$rbseff_i),
#                            mean(res2$rbseff_ii),
#                            mean(res2$rbseff_ibis),
#                            mean(res2$rbseff_iibis),
#                            mean(res2$rbseff1_iii),
#                            mean(res2$rbseff2_iii),
#                            mean(res2$rbseff1_iv),
#                            mean(res2$rbseff2_iv)), 
#                 rmse.seff = c(sqrt(mean(res1$rmseseff_i)),
#                              sqrt(mean(res1$rmseseff_ii)),
#                              sqrt(mean(res1$rmseseff_ibis)),
#                              sqrt(mean(res1$rmseseff_iibis)),
#                              sqrt(mean(res1$rmseseff1_iii)),
#                              sqrt(mean(res1$rmseseff2_iii)),
#                              sqrt(mean(res1$rmseseff1_iv)),
#                              sqrt(mean(res1$rmseseff2_iv)),
#                              sqrt(mean(res2$rmseseff_i)),
#                              sqrt(mean(res2$rmseseff_ii)),
#                              sqrt(mean(res2$rmseseff_ibis)),
#                              sqrt(mean(res2$rmseseff_iibis)),
#                              sqrt(mean(res2$rmseseff1_iii)),
#                              sqrt(mean(res2$rmseseff2_iii)),
#                              sqrt(mean(res2$rmseseff1_iv)),
#                              sqrt(mean(res2$rmseseff2_iv))))
#  
#  rseff <- rseff %>% mutate(type = rep(c("RV","SV","RV","SV","RV","RV","SV","SV"),2),
#                          occ = c(rep("Ψ = 0.1",length(model)/2),rep("Ψ = 0.3", length(model)/2)),
#                          name = rep(c("SOM p1","SOM p1","SOM p2","SOM p2","IOM p1","IOM p2","IOM p1","IOM p2"),2))
#  
#  pp1 <- rseff %>%
#    ggplot() +
#    aes(x = name, y = rmse.seff, fill = type) + 
#    geom_col(position = "dodge", width=.6) + 
#    labs(x = 'Model',
#         y = 'RMSE in sampling effort effect size') + 
#    scale_fill_manual(name = NULL, values=c("#442b48","#32a287"))+ 
#    facet_wrap(~ occ,  strip.position = "top", labeller = label_wrap_gen(multi_line = FALSE)) + 
#    theme_bw(base_size = 14)
#  
#  pp2 <- rseff %>%
#    ggplot() +
#    aes(x = name, y = rb.seff, fill = type) + 
#    geom_col(position = "dodge", width=.6) + 
#    labs(x = 'Model',
#         y = 'RB in sampling effort effect size') + 
#    scale_fill_manual(name = NULL, values=c("#442b48","#32a287"))+ 
#    facet_wrap(~ occ,  strip.position = "top", labeller = label_wrap_gen(multi_line = FALSE)) + 
#    theme_bw(base_size = 14)
#  
#  plot_grid(pocc1,pocc2,pp1,pp2, ncol = 2, nrow = 2)
#  
#  