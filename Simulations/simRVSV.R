# simulations only for RV and SV occupancy
# no mention of integrated occupancy

# pqckages 

library(R2jags)

#library(nimble)
#library(unmarked)
# load data
setwd("/home/vlauret/")

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

# ---- JAGS model ----
sink("Occu_cov.jags")
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

#---- loop for 100 simulations with JAGS ----
loopRVSV<- function(nsites, nocc, nsim, covparam, obsparam, Seff, Covenv){
  
  #### allocate memories for estimated parameters ####
  
  p_i <- rep(NA, nsim)
  psi_i <- rep(NA, nsim)

  bcov_i <- rep(NA, nsim)
  bseff_i <- rep(NA, nsim)

  p_ii <- rep(NA, nsim)
  psi_ii <- rep(NA, nsim)

  bcov_ii <- rep(NA, nsim)
  bseff_ii<- rep(NA, nsim)
  
  # allocate RMSE RB
  rbpsi_i    <- rep(NA, nsim)
  rmsepsi_i   <- rep(NA, nsim)
  rbcov_i     <- rep(NA, nsim)
  rmsecov_i <- rep(NA, nsim)
  rbp_i     <- rep(NA, nsim)
  rmsep_i    <- rep(NA, nsim)
  rbseff_i    <- rep(NA, nsim)
  rmseseff_i  <- rep(NA, nsim)

  rbpsi_ii    <- rep(NA, nsim)
  rmsepsi_ii  <- rep(NA, nsim)
  rbcov_ii    <- rep(NA, nsim)
  rmsecov_ii  <- rep(NA, nsim)
  rbp_ii      <- rep(NA, nsim)
  rmsep_ii    <- rep(NA, nsim)
  rbseff_ii   <- rep(NA, nsim)
  rmseseff_ii <- rep(NA, nsim)
  

  #### loop over nsim ####
  for( s in 1:nsim){
    print(s)
    
    # simulate datasets s
    
    z <- sim_state_cov(nb_sites = nsites, covenv= covenv, init_occ= covparam)
    y <- sim_obs_cov(z,nocc,Seff= Seff, obsparam= obsparam)
    
    p <- obsparam
    
    
# WITH Jags
   
    #### i = RV ####
                    
    # Bundle data
    win.data <- list(y = y,nsite = dim(y)[1], nocc = dim(y)[2], Seff = Seff, covenv = covenv)
    
    # Initial values
    zst <- apply(y, 1, max)	# Observed occurrence as inits for z
    
    # inits used for simulation
    psi = covparam[1]    # sum(z)/nsites
    p = obsparam[1]
    bcov <- covparam[2]
    bseff <- obsparam[2]
    
    inits <- function(){ list(z = zst, alpha.psi = psi, alpha.p= p, beta.cov = bcov, beta.seff= bseff)}
    
    # Parameters monitored
    params <- c( "alpha.p","alpha.psi", "beta.cov", "beta.qcov","beta.seff" ,"beta.qseff") 
    
    ni <- 4000
    nc <-   1
    nb <- 1000
      
    # Call JAGS from R (BRT 3 min)
    out <- jags(win.data, inits, params, "Occu_cov.jags", n.chains = nc, n.iter = ni, n.burnin = nb, working.directory = getwd())
    
     #head(out$BUGSoutput$summary)
     #denplot(out)
  #effectiveSize(as.mcmc(out))
    # intercept
    p_i[s] <- out$BUGSoutput$mean$alpha.p      
    psi_i[s] <- out$BUGSoutput$mean$alpha.psi

    # slope linear estimate (Q: do we have to save the quadratic estimates-)
    
    bcov_i[s] <- out$BUGSoutput$mean$beta.cov   
    bseff_i[s] <- out$BUGSoutput$mean$beta.seff

    #### ii = SV ####
 
    ysv <-  matrix(apply(y,1,max))
    
    Seffsv <- matrix(apply(Seff,1,sum))
    
    # Bundle data
    win.data <- list(y = ysv,nsite = dim(ysv)[1], nocc = dim(ysv)[2], Seff = Seff, covenv = covenv)
    
    # Initial values
    zst <- apply(ysv, 1, max)	# Observed occurrence as inits for z
    
    # inits used for simulation
    psi = covparam[1]    # sum(z)/nsites
    p = obsparam[1]
    bcov <- covparam[2]
    bseff <- obsparam[2]
    
    inits <- function(){ list(z = zst, alpha.psi = psi, alpha.p= p, beta.cov = bcov, beta.seff= bseff)}
    
    # Parameters monitored
    params <- c( "alpha.p","alpha.psi", "beta.cov", "beta.qcov","beta.seff" ,"beta.qseff") 
    

    # Call JAGS from R (BRT 3 min)
    out <- jags(win.data, inits, params, "Occu_cov.jags", n.chains = nc, n.iter = ni, n.burnin = nb, working.directory = getwd())
    
    # head(out$BUGSoutput$summary)
    # denplot(out)
    
    # intercept
    p_ii[s] <- out$BUGSoutput$mean$alpha.p       
    psi_ii[s] <- out$BUGSoutput$mean$alpha.psi 
    
    # slope linear estimate (Q: do we have to save the quadratic estimates-)
    
    bcov_ii[s] <- out$BUGSoutput$mean$beta.cov
    bseff_ii[s] <- out$BUGSoutput$mean$beta.seff

      # save --- RB & SE   
    # i 
    rbpsi_i[s]        <- (psi_i[s]-covparam[1])/covparam[1]   # RB alpha.psi
    rmsepsi_i[s]      <- (psi_i[s] - covparam[1])^2 # RMSE alpha.psi
    
    rbcov_i[s]        <- (bcov_i[s]-covparam[2])/covparam[2]   # RB cov
    rmsecov_i[s]      <- (bcov_i[s] - covparam[2])^2 # RMSE cov
    
    rbp_i[s]          <- (p_i[s]-obsparam[1])/obsparam[1]      # RB alpha p 
    rmsep_i[s]        <- (p_i[s] - obsparam[1])^2 # RMSE alpha p
    
    rbseff_i[s]       <- (bseff_i[s]-obsparam[2])/obsparam[2]   # RB obs
    rmseseff_i[s]     <- (bseff_i[s] - obsparam[2])^2 # RMSE obs
    
    
    # ii
    rbpsi_ii[s]       <- (psi_ii[s]-covparam[1])/covparam[1]   # RB alpha.psi
    rmsepsi_ii[s]     <- (psi_ii[s] - covparam[1])^2 # RMSE alpha.psi
    
    rbcov_ii[s]       <- (bcov_ii[s]-covparam[2])/covparam[2]   # RB cov
    rmsecov_ii[s]     <- (bcov_ii[s] - covparam[2])^2 # RMSE cov
    
    rbp_ii[s]         <- (p_ii[s]-obsparam[1])/obsparam[1]       # RB alpha p 
    rmsep_ii[s]       <- (p_ii[s] - obsparam[1])^2 # RMSE alpha p
    
    rbseff_ii[s]      <- (bseff_ii[s]-obsparam[2])/obsparam[2]   # RB obs
    rmseseff_ii[s]    <- (bseff_ii[s] - obsparam[2])^2 # RMSE obs
    
  
  } # end for 
  
  

  para <- data.frame("nsites" = nsites,"nocc" =  nocc, "nsim" = nsim, "covparam"= covparam, "obsparam"= obsparam)
  
  return(list(para, # i 
              rbpsi_i,
              rmsepsi_i,
              rbcov_i,
              rmsecov_i,
              rbp_i,
              rmsep_i,
              rbseff_i,
              rmseseff_i,
              #,
              rbpsi_ii,
              rmsepsi_ii,
              rbcov_ii,
              rmsecov_ii,
              rbp_ii,
              rmsep_ii,
              rbseff_ii,
              rmseseff_ii))
  
} # endfunction


# ---- initialisation ----

nsim = 100 # number of simulations

nsites = 100 # number of sites
nocc = 4 # number of occasions

# parameters of the logit relation between Psi and fictive covariate logit(Psi) <- a0 + a1 * cov + a2 * cov * cov
covparam.3 <- c(-0.5,0.2,0) # psi = 0.3
covparam.1 <- c(-1.9,0.2,0) # psi = 0.1
covparam.5 <- c(0.5,0.23,0) # psi = 0.5
covparam.9 <- c(2.5,0.15,0) # psi = 0.9

covparam <- list(covparam.1 =  covparam.1,
                 covparam.3 =  covparam.3, 
                 covparam.5 =  covparam.5,
                 covparam.9 <- covparam.9)

# ---- plot relations psi = f(cov) ----
#       covenv <- runif(100,0,1)
#       cov_rel.1 <- as.numeric(covparam.1[1]+covparam.1[2]*scale(covenv)+covparam.1[3]*scale(covenv)*scale(covenv))
#       cov_rel.3 <- as.numeric(covparam.3[1]+covparam.3[2]*scale(covenv)+covparam.3[3]*scale(covenv)*scale(covenv))
#       cov_rel.5 <- as.numeric(covparam.5[1]+covparam.5[2]*scale(covenv)+covparam.5[3]*scale(covenv)*scale(covenv))
#       cov_rel.9 <- as.numeric(covparam.9[1]+covparam.9[2]*scale(covenv)+covparam.9[3]*scale(covenv)*scale(covenv))
#       
#       ggplot() + geom_line(aes(x=covenv,y=1/(1+exp(-cov_rel.3)), color="0.3"), lwd= 2) +
#         geom_line(aes(x=covenv,y=1/(1+exp(-cov_rel.1)),color="0.1"), lwd= 2) +
#         geom_line(aes(x=covenv,y=1/(1+exp(-cov_rel.5)),color="0.5"), lwd= 2) +
#         geom_line(aes(x=covenv,y=1/(1+exp(-cov_rel.9)),color="0.9"), lwd= 2) +
#         ylim(0,1) + ylab("Occupancy probability : Ψ") +xlab("Environmental covariate") + theme_minimal() + scale_color_viridis_d()
# ----

# parameters of the logit relation between p and fictive sampling effort logit(p) <- b0 + ab * seff + a2 * seff * seff

obsparam.5 <- c(0.3,0.26,0) # p = 0.5
obsparam.8 <- c(1.8,0.3,0) # p = 0.8
obsparam.1 <- c(-1.5,0.26,0) # p = 0.5
obsparam.3 <- c(-0.6,0.25,0)

#  plot 
#     obs_rel.1 <- as.numeric(obsparam.1[1]+obsparam.1[2]*Seff.8+obsparam.1[3]*Seff.8*Seff.8)
#     obs_rel.3 <- as.numeric(obsparam.3[1]+obsparam.3[2]*Seff.8+obsparam.3[3]*Seff.8*Seff.8)
#     obs_rel.5 <- as.numeric(obsparam.5[1]+obsparam.5[2]*Seff.8+obsparam.5[3]*Seff.8*Seff.8)
#     obs_rel.8 <- as.numeric(obsparam.8[1]+obsparam.8[2]*Seff.8+obsparam.8[3]*Seff.8*Seff.8)
#     
#     ggplot() + geom_line(aes(x=Seff.8,y=1/(1+exp(-obs_rel.8)), color="0.8"), lwd= 2) +
#       geom_line(aes(x=Seff.8,y=1/(1+exp(-obs_rel.5)),color="0.5"), lwd= 2) +
#       geom_line(aes(x=Seff.8,y=1/(1+exp(-obs_rel.1)),color="0.1"), lwd= 2) +
#       geom_line(aes(x=Seff.8,y=1/(1+exp(-obs_rel.3)),color="0.3"), lwd= 2) +
#       ylim(0,1) + ylab("Detection probability : p") +xlab("Sampling effort covariate") + theme_minimal() + scale_color_viridis_d(begin = 0.2, end = 0.9)

obsparam <- list(obsparam.1 <- obsparam.1,# p = 0.5
                 obsparam.3 <- obsparam.3,
                 obsparam.5 <- obsparam.5, # p = 0.5
                 obsparam.8 <- obsparam.8# p = 0.8
                 )

# generate covariate values 

  Seff <- data.frame(cbind(rnorm(nsites),rnorm(nsites),rnorm(nsites),rnorm(nsites)))
  
  Covenv <- covenv <- runif(nsites,0,1)
  
  res <- list()
  i <- 0
  t <- proc.time()
 for(c in 1: length(covparam)){  #i <- i+1
 
   for(o in 1:length(obsparam)){
    i <- i+1 
    
   print(paste("obs",o))
   print(paste("cov",c))
   res[[i]] <- loopRVSV(nsites = nsites, nocc = nocc, nsim = nsim, covparam = covparam[[c]], obsparam[[o]], Seff, Covenv)

   
   }# end o
   } # end c

  run.time <- proc.time() - t
#
# ------ save results if you want
save(res,run.time, file ="simRVSV.Rdata")

  
# names res element as the 16 combinations -->  "Psi:p"
# names(res) <- c("0.1:0.1","0.1:0.3","0.1:0.5","0.1:0.8",
#                 "0.3:0.1","0.3:0.3","0.3:0.5","0.3:0.8",
#                 "0.5:0.1","0.5:0.3","0.5:0.5","0.5:0.8",
#                 "0.9:0.1","0.9:0.3","0.9:0.5","0.9:0.8")
# 
## ---- verif parametres aux bornes ----
 
 verif <- list(a = res$`0.9:0.1`, b = res$`0.9:0.3`)
 names(verif$a) <-  names(verif$b) <-  c("para", "rbpsi_i","rmsepsi_i" , "rbcov_i" ,
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
                     "rmseseff_ii")
 para.rv <- verif$b$rbcov_i*0.15 + 0.15
 para.rv
 
 para.sv <- verif$b$rbcov_ii*0.15 + 0.15
 para.sv


## ---- plot the results as in the html file ----
# library(tidyverse)
# library(cowplot)
# 
#
#
# nsim = 100
#
## make a tibble to store values 
#extract.RMSE <- function(res, para){
#  m <- c()
#  for(i in 1:length(res)){
#    m[i] <- sqrt(mean(res[[i]][[para]]))
#  } # i
#  return(m)
#}
#
#extract.RB <- function(res, para){
#  m <- c()
#  for(i in 1:length(res)){
#    m[i] <- mean(res[[i]][[para]])
#  } # i
#  return(m)
#}
# 
#RB.psi <- tibble(model =c("0.1:0.1","0.1:0.3","0.1:0.5","0.1:0.8",
#                          "0.3:0.1","0.3:0.3","0.3:0.5","0.3:0.8",
#                          "0.5:0.1","0.5:0.3","0.5:0.5","0.5:0.8",
#                          "0.9:0.1","0.9:0.3","0.9:0.5","0.9:0.8"),
#                 psi.interceptRV = extract.RB(res,2),
#                 psi.interceptSV = extract.RB(res,10),
#                 psi.covRV = extract.RB(res,4),
#                 psi.covSV = extract.RB(res,12)) 
# 
#RMSE.psi <- tibble(model =c("0.1:0.1","0.1:0.3","0.1:0.5","0.1:0.8",
#                            "0.3:0.1","0.3:0.3","0.3:0.5","0.3:0.8",
#                            "0.5:0.1","0.5:0.3","0.5:0.5","0.5:0.8",
#                            "0.9:0.1","0.9:0.3","0.9:0.5","0.9:0.8"),
#                   psi.interceptRV = sqrt(extract.RMSE(res,3)),
#                   psi.interceptSV = sqrt(extract.RMSE(res,11)),
#                   psi.covRV = sqrt(extract.RMSE(res,5)),
#                   psi.covSV = sqrt(extract.RMSE(res,13)))
#
#
#RB.p <- tibble(model =c("0.1:0.1","0.1:0.3","0.1:0.5","0.1:0.8",
#                        "0.3:0.1","0.3:0.3","0.3:0.5","0.3:0.8",
#                        "0.5:0.1","0.5:0.3","0.5:0.5","0.5:0.8",
#                        "0.9:0.1","0.9:0.3","0.9:0.5","0.9:0.8"),
#               p.interceptRV = extract.RB(res,6),
#               p.interceptSV = extract.RB(res,14),
#               p.covRV = extract.RB(res,8),
#               p.covSV = extract.RB(res,16))
#
#RMSE.p <- tibble(model =c("0.1:0.1","0.1:0.3","0.1:0.5","0.1:0.8",
#                          "0.3:0.1","0.3:0.3","0.3:0.5","0.3:0.8",
#                          "0.5:0.1","0.5:0.3","0.5:0.5","0.5:0.8",
#                          "0.9:0.1","0.9:0.3","0.9:0.5","0.9:0.8"),
#                 p.interceptRV = sqrt(extract.RMSE(res,7)),
#                 p.interceptSV = sqrt(extract.RMSE(res,15)),
#                 p.covRV = sqrt(extract.RMSE(res,9)),
#                 p.covSV = sqrt(extract.RMSE(res,17)))
#
#
#  
#r <- tibble(model = rep(RMSE.psi$model, 2),rmse.occ = c(RMSE.psi$psi.covRV,RMSE.psi$psi.covSV), rmse.p = c(RMSE.p$p.covRV,RMSE.p$p.covSV), 
#            rb.occ = c(RB.psi$psi.covRV,RB.psi$psi.covSV), rb.p = c(RB.p$p.covRV,RB.p$p.covSV)) %>% 
#  mutate(type = c(rep("RV",length(model)/2),rep("SV",length(model)/2)),
#         occ = c(rep("Ψ = 0.1", 4),rep("Ψ = 0.3", 4),rep("Ψ = 0.5", 4),rep("Ψ = 0.9", 4),
#                 rep("Ψ = 0.1", 4),rep("Ψ = 0.3", 4),rep("Ψ = 0.5", 4),rep("Ψ = 0.9", 4)),
#         det = rep(c("0.1", "0.3","0.5","0.8"), 8))
#             
#
#pocc1 <- r %>%
#  ggplot() +
#  aes(x = det, y = rmse.occ, fill = type) + 
#  geom_col(position = "dodge", width=.6) + 
#  labs(x = 'Detection probability',
#       y = 'RMSE in occupancy covariate effect size') + 
#  scale_fill_manual(name = NULL, values=c("#442b48","#32a287"))+ 
#    facet_wrap(~ occ,  strip.position = "top", labeller = label_wrap_gen(multi_line = FALSE)) + 
#  theme_bw(base_size = 14)
#
#pocc2 <- r %>%
#  ggplot() +
#  aes(x = det, y = rb.occ, fill = type) + 
#  geom_col(position = "dodge", width=.6) + 
#  labs(x = 'Detection probability',
#       y = 'RB in occupancy covariate effect size') + 
#  scale_fill_manual(name = NULL, values=c("#442b48","#32a287"))+ 
#  facet_wrap(~ occ,  strip.position = "top", labeller = label_wrap_gen(multi_line = FALSE)) + 
#  theme_bw(base_size = 14)
#
#pp1 <- r %>%
#  ggplot() +
#  aes(x = det, y = rmse.p, fill = type) + 
#  geom_col(position = "dodge", width=.6) + 
#  labs(x = 'Detection probability',
#       y = 'RMSE in sampling effort effect size') + 
#  scale_fill_manual(name = NULL, values=c("#442b48","#32a287"))+ 
#  facet_wrap(~ occ,  strip.position = "top", labeller = label_wrap_gen(multi_line = FALSE)) + 
#  theme_bw(base_size = 14)
#
#pp2 <- r %>%
#  ggplot() +
#  aes(x = det, y = rb.p, fill = type) + 
#  geom_col(position = "dodge", width=.6) + 
#  labs(x = 'Detection probability',
#       y = 'RB in sampling effort effect size') + 
#  scale_fill_manual(name = NULL, values=c("#442b48","#32a287"))+ 
#  facet_wrap(~ occ,  strip.position = "top", labeller = label_wrap_gen(multi_line = FALSE)) + 
#  theme_bw(base_size = 14)
#
#plot_grid(pocc1,pocc2,pp1,pp2, ncol = 2, nrow = 2)
# 
#
#
##--------- Intercept -----
#
#r <- tibble(model = rep(RMSE.psi$model, 2),rmse.occ = c(RMSE.psi$psi.interceptRV,RMSE.psi$psi.interceptSV), rmse.p = c(RMSE.p$p.interceptRV,RMSE.p$p.interceptSV), 
#            rb.occ = c(RB.psi$psi.interceptRV,RB.psi$psi.interceptSV), rb.p = c(RB.p$p.interceptRV,RB.p$p.interceptSV)) %>% 
#  mutate(type = c(rep("RV",length(model)/2),rep("SV",length(model)/2)),
#         occ = c(rep("Ψ = 0.1", 4),rep("Ψ = 0.3", 4),rep("Ψ = 0.5", 4),rep("Ψ = 0.9", 4),
#                 rep("Ψ = 0.1", 4),rep("Ψ = 0.3", 4),rep("Ψ = 0.5", 4),rep("Ψ = 0.9", 4)),
#         det = rep(c("0.1", "0.3","0.5","0.8"), 8))
#
#
#pocc1 <- r %>%
#  ggplot() +
#  aes(x = det, y = rmse.occ, fill = type) + 
#  geom_col(position = "dodge", width=.6) + 
#  labs(x = 'Detection probability',
#       y = 'RMSE in occupancy intercept') + 
#  scale_fill_manual(name = NULL, values=c("#442b48","#32a287"))+ 
#  facet_wrap(~ occ,  strip.position = "top", labeller = label_wrap_gen(multi_line = FALSE)) + 
#  theme_bw(base_size = 14)
#
#pocc2 <- r %>%
#  ggplot() +
#  aes(x = det, y = rb.occ, fill = type) + 
#  geom_col(position = "dodge", width=.6) + 
#  labs(x = 'Detection probability',
#       y = 'RB in occupancy intercept') + 
#  scale_fill_manual(name = NULL, values=c("#442b48","#32a287"))+ 
#  facet_wrap(~ occ,  strip.position = "top", labeller = label_wrap_gen(multi_line = FALSE)) + 
#  theme_bw(base_size = 14)
#
#pp1 <- r %>%
#  ggplot() +
#  aes(x = det, y = rmse.p, fill = type) + 
#  geom_col(position = "dodge", width=.6) + 
#  labs(x = 'Detection probability',
#       y = 'RMSE in detection prob. intercept') + 
#  scale_fill_manual(name = NULL, values=c("#442b48","#32a287"))+ 
#  facet_wrap(~ occ,  strip.position = "top", labeller = label_wrap_gen(multi_line = FALSE)) + 
#  theme_bw(base_size = 14)
#
#pp2 <- r %>%
#  ggplot() +
#  aes(x = det, y = rb.p, fill = type) + 
#  geom_col(position = "dodge", width=.6) + 
#  labs(x = 'Detection probability',
#       y = 'RB in detection prob. intercept') + 
#  scale_fill_manual(name = NULL, values=c("#442b48","#32a287"))+ 
#  facet_wrap(~ occ,  strip.position = "top", labeller = label_wrap_gen(multi_line = FALSE)) + 
#  theme_bw(base_size = 14)
#
#
#plot_grid(pocc1,pocc2,pp1,pp2, ncol = 2, nrow = 2)
#
#
## res1 low occupancy 100 sites
#Res1_cov<- tibble(model = c(rep("RV low p",nsim),rep("SV low p",nsim),rep("RV high p",nsim),rep("SV high p",nsim)),
#                  RMSE = c(res1$rmsecov_i^2, res1$rmsecov_ii^2,res1$rmsecov_ibis^2, res1$rmsecov_iibis^2), 
#                  RB = c(res1$rbcov_i, res1$rbcov_ii, res1$rbcov_ibis, res1$rbcov_iibis))
#
#Res1_seff <- tibble(model = c(rep("RV low p",nsim),rep("SV low p",nsim),rep("RV high p",nsim),rep("SV high p",nsim)),
#                   RMSE = c(res1$rmseseff_i^2, res1$rmseseff_ii^2,res1$rmseseff_ibis^2, res1$rmseseff_iibis^2), 
#                   RB = c(res1$rbseff_i, res1$rbseff_ii, res1$rbseff_ibis, res1$rbseff_iibis))
#
## res2 high occupancy 100 sites
#Res2_cov<- tibble(model = c(rep("RV low p",nsim),rep("SV low p",nsim),rep("RV high p",nsim),rep("SV high p",nsim)),
#                  RMSE = c(res2$rmsecov_i^2, res2$rmsecov_ii^2,res2$rmsecov_ibis^2, res2$rmsecov_iibis^2), 
#                  RB = c(res2$rbcov_i, res2$rbcov_ii, res2$rbcov_ibis, res2$rbcov_iibis))
#
#Res2_seff <- tibble(model = c(rep("RV low p",nsim),rep("SV low p",nsim),rep("RV high p",nsim),rep("SV high p",nsim)),
#                    RMSE = c(res2$rmseseff_i^2, res2$rmseseff_ii^2,res2$rmseseff_ibis^2, res2$rmseseff_iibis^2), 
#                    RB = c(res2$rbseff_i, res2$rbseff_ii, res2$rbseff_ibis, res2$rbseff_iibis))
#
## res3 low occupancy 400 sites
#Res3_cov<- tibble(model = c(rep("RV low p",nsim),rep("SV low p",nsim),rep("RV high p",nsim),rep("SV high p",nsim)),
#                  RMSE = c(res3$rmsecov_i^2, res3$rmsecov_ii^2,res3$rmsecov_ibis^2, res3$rmsecov_iibis^2), 
#                  RB = c(res3$rbcov_i, res3$rbcov_ii, res3$rbcov_ibis, res3$rbcov_iibis))
#
#Res3_seff <- tibble(model = c(rep("RV low p",nsim),rep("SV low p",nsim),rep("RV high p",nsim),rep("SV high p",nsim)),
#                    RMSE = c(res3$rmseseff_i^2, res3$rmseseff_ii^2,res3$rmseseff_ibis^2, res3$rmseseff_iibis^2), 
#                    RB = c(res3$rbseff_i, res3$rbseff_ii, res3$rbseff_ibis, res3$rbseff_iibis))
## res4 low occupancy 400 sites
#Res4_cov<- tibble(model = c(rep("RV low p",nsim),rep("SV low p",nsim),rep("RV high p",nsim),rep("SV high p",nsim)),
#                  RMSE = c(res4$rmsecov_i^2, res4$rmsecov_ii^2,res4$rmsecov_ibis^2, res4$rmsecov_iibis^2), 
#                  RB = c(res4$rbcov_i, res4$rbcov_ii, res4$rbcov_ibis, res4$rbcov_iibis))
#
#Res4_seff <- tibble(model = c(rep("RV low p",nsim),rep("SV low p",nsim),rep("RV high p",nsim),rep("SV high p",nsim)),
#                    RMSE = c(res4$rmseseff_i^2, res4$rmseseff_ii^2,res4$rmseseff_ibis^2, res4$rmseseff_iibis^2), 
#                    RB = c(res4$rbseff_i, res4$rbseff_ii, res4$rbseff_ibis, res4$rbseff_iibis))
## tableau of RB and RMSE
#
#
##
##
#RMSE_cov <- tibble( r1 = c(mean(Res1_cov$RMSE[Res1_cov$model=="RV low p"]),
#                        mean(Res1_cov$RMSE[Res1_cov$model=="SV low p"]),
#                        mean(Res1_cov$RMSE[Res1_cov$model=="IOM RV"]),
#                        mean(Res1_cov$RMSE[Res1_cov$model=="IOM SV"])), 
#                 r2 = c(mean(Res2_cov$RMSE[Res2_cov$model=="RV low p"]),
#                        mean(Res2_cov$RMSE[Res2_cov$model=="SV low p"]),
#                        mean(Res2_cov$RMSE[Res2_cov$model=="IOM RV"]),
#                        mean(Res2_cov$RMSE[Res2_cov$model=="IOM SV"])),
#                 r3 = c(mean(Res3_cov$RMSE[Res3_cov$model=="RV low p"]),
#                        mean(Res3_cov$RMSE[Res3_cov$model=="SV low p"]),
#                        mean(Res3_cov$RMSE[Res3_cov$model=="IOM RV"]),
#                        mean(Res3_cov$RMSE[Res3_cov$model=="IOM SV"])),
#                 r4 = c(mean(Res4_cov$RMSE[Res4_cov$model=="RV low p"]),
#                        mean(Res4_cov$RMSE[Res4_cov$model=="SV low p"]),
#                        mean(Res4_cov$RMSE[Res4_cov$model=="IOM RV"]),
#                        mean(Res4_cov$RMSE[Res4_cov$model=="IOM SV"])),
#                 model = c("RV","SV","IOM RV","IOM SV"))
#
#RMSE_seff <- tibble( r1 =c(mean(Res1_seff$RMSE[Res1_seff$model=="RV low p"]),
#                        mean(Res1_seff$RMSE[Res1_seff$model=="SV low p"]),
#                        mean(Res1_seff$RMSE[Res1_seff$model=="IOM RV low p"]),
#                        mean(Res1_seff$RMSE[Res1_seff$model=="IOM RV high p"]),
#                        mean(Res1_seff$RMSE[Res1_seff$model=="IOM SV low p"]),
#                        mean(Res1_seff$RMSE[Res1_seff$model=="IOM SV high p"])), 
#                 r2 = c(mean(Res2_seff$RMSE[Res2_seff$model=="RV low p"]),
#                        mean(Res2_seff$RMSE[Res2_seff$model=="SV low p"]),
#                        mean(Res2_seff$RMSE[Res2_seff$model=="IOM RV low p"]),
#                        mean(Res2_seff$RMSE[Res2_seff$model=="IOM RV high p"]),
#                        mean(Res2_seff$RMSE[Res2_seff$model=="IOM SV low p"]),
#                        mean(Res2_seff$RMSE[Res2_seff$model=="IOM SV high p"])),
#                 r3 = c(mean(Res3_seff$RMSE[Res3_seff$model=="RV low p"]),
#                        mean(Res3_seff$RMSE[Res3_seff$model=="SV low p"]),
#                        mean(Res3_seff$RMSE[Res3_seff$model=="IOM RV low p"]),
#                        mean(Res3_seff$RMSE[Res3_seff$model=="IOM RV high p"]),
#                        mean(Res3_seff$RMSE[Res3_seff$model=="IOM SV low p"]),
#                        mean(Res3_seff$RMSE[Res3_seff$model=="IOM SV high p"])),
#                 r4 = c(mean(Res4_seff$RMSE[Res4_seff$model=="RV low p"]),
#                        mean(Res4_seff$RMSE[Res4_seff$model=="SV low p"]),
#                        mean(Res4_seff$RMSE[Res4_seff$model=="IOM RV low p"]),
#                        mean(Res4_seff$RMSE[Res4_seff$model=="IOM RV high p"]),
#                        mean(Res4_seff$RMSE[Res4_seff$model=="IOM SV low p"]),
#                        mean(Res4_seff$RMSE[Res4_seff$model=="IOM SV high p"])),
#                 model = c("RV","SV","IOM RV low p","IOM RV p 2","IOM SV low p", "IOM SV high p")) 
##
#save(Res1_cov,Res1_seff,
#     Res2_cov,Res2_seff,
#     Res3_cov,Res3_seff,
#     Res4_cov,Res4_seff,file = "simulunmarked.rdata")
# plot 
##100 sites violin
#RMSEcov1 <- bind_rows(Res1_cov,Res2_cov) %>% mutate(occ = c(rep("ψ = 0.1",nrow(Res1_cov)),rep("ψ = 0.3",nrow(Res2_cov))))
##
#p_RMSEcov1 <- ggplot() + geom_violin(data = RMSEcov1, aes(x = model, y = RMSE, fill = model),alpha= 0.7,draw_quantiles = 0.5) + 
# ylim(0,0.5) + scale_fill_viridis_d()  + theme_cowplot() +facet_wrap(~ occ, strip.position = "bottom") +
# scale_x_discrete(labels = rep("",4)) +
# labs(title= "Occupancy covariate effect size",subtitle= paste("100 sites,", nsim,"simulations")) + xlab("")
#p_RMSEcov1
##
#RMSEseff1 <- bind_rows(Res1_seff,Res2_seff) %>% mutate(occ = c(rep("ψ = 0.1",nrow(Res1_seff)),rep("ψ = 0.3",nrow(Res2_seff))))
##
#p_RMSEseff1 <- ggplot() + geom_violin(data = RMSEseff1, aes(x = model, y = RMSE, fill = model),alpha= 0.7,draw_quantiles = 0.5, scale ="width") + 
# ylim(0,0.5) + scale_fill_viridis_d()  + theme_cowplot() +facet_wrap(~ occ,  strip.position = "bottom") +
# scale_x_discrete(labels = rep("",6)) +
# labs(title= "Observation covariate effect size",subtitle= paste("100 sites,", nsim,"simulations"))  + xlab("")
#p_RMSEseff1
## 400 sites violin
#RMSEcov4 <- bind_rows(Res3_cov,Res4_cov) %>% mutate(occ = c(rep("ψ = 0.1",nrow(Res3_cov)),rep("ψ = 0.3",nrow(Res4_cov))))
##
#p_RMSEcov4 <- ggplot() + geom_violin(data = RMSEcov4, aes(x = model, y = RMSE, fill = model),alpha= 0.7,draw_quantiles = 0.5) + 
# ylim(0,0.5) + scale_fill_viridis_d()  + theme_cowplot() +facet_wrap(~ occ, strip.position = "bottom") +
# scale_x_discrete(labels = rep("",4)) +
# labs(title= "",subtitle= paste("400 sites,", nsim,"simulations")) + theme(legend.position = "none") + xlab("")
#p_RMSEcov4
##
#RMSEseff4 <- bind_rows(Res3_seff,Res4_seff) %>% mutate(occ = c(rep("ψ = 0.1",nrow(Res3_seff)),rep("ψ = 0.3",nrow(Res4_seff))))
##
#p_RMSEseff4 <- ggplot() + geom_violin(data = RMSEseff4, aes(x = model, y = RMSE, fill = model),alpha= 0.7,draw_quantiles = 0.5,scale = "width") + 
# ylim(0,0.5) + scale_fill_viridis_d()  + theme_cowplot() +facet_wrap(~ occ , strip.position = "bottom") +
# scale_x_discrete(labels = rep("",6)) +
# labs(title= "",subtitle= paste("400 sites,", nsim,"simulations")) + theme(legend.position = "none") + xlab("")
#p_RMSEseff4
##
#
#library(cowplot)
#plot_grid(p_RMSEcov1,p_RMSEcov4, rel_widths = c(1.25,1))
#plot_grid(p_RMSEseff1,p_RMSEseff4,rel_widths = c(1.25,1))
#
## RB
##100 sites violin
#RBcov1 <- bind_rows(Res1_cov,Res2_cov) %>% mutate(occ = c(rep("ψ = 0.1",nrow(Res1_cov)),rep("ψ = 0.3",nrow(Res2_cov))))
##
#p_RBcov1 <- ggplot() + geom_violin(data = RBcov1, aes(x = model, y = RB, fill = model),alpha= 0.7,draw_quantiles = 0.5)  +
#  scale_fill_viridis_d()  + theme_cowplot() +facet_wrap(~ occ, strip.position = "bottom") +
#  scale_x_discrete(labels = rep("",4)) +
#  labs(title= "Occupancy covariate effect size",subtitle= paste("100 sites,", nsim,"simulations")) + xlab("")
#p_RBcov1
##
#RBseff1 <- bind_rows(Res1_seff,Res2_seff) %>% mutate(occ = c(rep("ψ = 0.1",nrow(Res1_seff)),rep("ψ = 0.3",nrow(Res2_seff))))
##
#p_RBseff1 <- ggplot() + geom_violin(data = RBseff1, aes(x = model, y = RB, fill = model),alpha= 0.7,draw_quantiles = 0.5, scale ="width") + 
#   scale_fill_viridis_d()  + theme_cowplot() +facet_wrap(~ occ,  strip.position = "bottom") +
#  scale_x_discrete(labels = rep("",6)) +
#  labs(title= "Observation covariate effect size",subtitle= paste("100 sites,", nsim,"simulations"))  + xlab("")
#p_RBseff1
## 400 sites violin
#RBcov4 <- bind_rows(Res3_cov,Res4_cov) %>% mutate(occ = c(rep("ψ = 0.1",nrow(Res3_cov)),rep("ψ = 0.3",nrow(Res4_cov))))
##
#p_RBcov4 <- ggplot() + geom_violin(data = RBcov4, aes(x = model, y = RB, fill = model),alpha= 0.7,draw_quantiles = 0.5) + 
#  scale_fill_viridis_d()  + theme_cowplot() +facet_wrap(~ occ, strip.position = "bottom") +
#  scale_x_discrete(labels = rep("",4)) +
#  labs(title= "",subtitle= paste("400 sites,", nsim,"simulations")) + theme(legend.position = "none") + xlab("")
#p_RBcov4
##
#RBseff4 <- bind_rows(Res3_seff,Res4_seff) %>% mutate(occ = c(rep("ψ = 0.1",nrow(Res3_seff)),rep("ψ = 0.3",nrow(Res4_seff))))
##
#p_RBseff4 <- ggplot() + geom_violin(data = RBseff4, aes(x = model, y = RB, fill = model),alpha= 0.7,draw_quantiles = 0.5,scale = "width") + 
# scale_fill_viridis_d()  + theme_cowplot() +facet_wrap(~ occ , strip.position = "bottom") +
#  scale_x_discrete(labels = rep("",6)) +
#  labs(title= "",subtitle= paste("400 sites,", nsim,"simulations")) + theme(legend.position = "none") + xlab("")
#p_RBseff4
##
#
#plot_grid(p_RMSEcov1,p_RMSEcov4, p_RBcov1,p_RBcov4, ncol = 2, nrow = 2,rel_widths = c(1.25,1))
#plot_grid(p_RMSEseff1,p_RMSEseff4, p_RBseff1,p_RBseff4, ncol = 2, nrow = 2,rel_widths = c(1.25,1))
#
## plot the following line when multiple simulations performed
#
#p_RBcov <- ggplot() + geom_boxplot(data = Res1_cov, aes(x = model, y = RB, fill = model),alpha= 0.5) + 
#  ylim(0,max(Res1_cov$RB)) + scale_fill_viridis_d()  + theme_minimal() + labs(title= paste("Relative bias of state cov slope :", nsim,"simulations"))
#p_RBcov
#
#p_RBcov <- ggplot() + geom_boxplot(data = Res1_cov, aes(x = model, y = RB, fill = model),alpha= 0.5) + 
#  ylim(0,max(Res1_cov$RMSE)) + scale_fill_viridis_d()  + theme_minimal() + labs(title= paste("RMSE of state cov slope :", nsim,"simulations"))
#p_RMSEcov
#
## p_CVp <- ggplot() + geom_point(data = Res_p, aes(x = model, y = CV, color = model)) + 
##   ylim(0,2.2) + scale_color_viridis_d() + theme_minimal()
## p_CVp
#
#p_RMSEseff <- ggplot() + geom_boxplot(data = Res1_seff, aes(x = model, y = RMSE, fill = model),alpha= 0.5) + 
#  ylim(0,max(Res1_seff$RMSE)) + scale_fill_viridis_d()  + theme_minimal() + labs(title= paste("RMSE of sampling effort :", nsim,"simulations"))
#p_RMSEseff
#
#p_RBseff <- ggplot() + geom_boxplot(data = Res1_seff, aes(x = model, y = RB, fill = model),alpha= 0.5) + 
#  ylim(0,max(Res1_seff$RB)) + scale_fill_viridis_d()  + theme_minimal() + labs(title= paste("RB of sampling effort :", nsim,"simulations"))
#p_RBseff
#
## make grid 
#
#p_CVRB <- plot_grid(p_RMSEcov, p_RBcov, p_RMSEseff, p_RBseff, ncol =2, nrow = 2)
#p_CVRB
#