
library(tidyverse)
library(PlackettLuce)
library(ClimMobTools)
library(mvtnorm)
library(parallel)
library(gtools)
library(stringi)
library(asreml)
library(emmeans)
library(HelpersMG)
library(complexlm)
library(lme4)
library(gosset)

setwd("Paper_1_Estimation_of_genetic_parameters/")

setwd('/home/joost/Hugo_Projects/Sim_est_genotypic_values/') 

source('Scripts/0. Ranking_pars_FUN.R')

# Add reference variety

data_30_exmple <- read.csv('Results/Plot_inputs/dataset_summary_out_30_ref_case.csv')

data_30_exmple <- data_30_exmple %>% filter(iter==3 & model=='MM') %>% dplyr::select("genotype","value.genotype")


repetition <- expand.grid(iter=1:100,n.rep.genotype.env=c(30),sigma.plot=0.9)


#repetition <- data.frame(seed=1:nrow(experiment),experiment)

mu.genotype <- 3.63

#-------------------------------------------------------------------------------

simulate_multi_env_tricot_V2 <- function(n.genotypes=20,mu.genotype =NULL,n.envs =1,
                                         sigma.genotype =NULL,sigma.env.farm=3.352,sigma.env=0,
                                         sigma.gen.env=0,sigma.plot=0.9,n.rep.genotype.env=30,genotypic.means = NULL ){
  
  #  
  
  if(!is.null(genotypic.means)){
    n.genotypes <- length(genotypic.means)
    mu.genotype <- mean(genotypic.means)
    sigma.genotype <- sd(genotypic.means)
  }
  
  require(PlackettLuce)
  require(ClimMobTools)
  
  n.packages <- round((n.genotypes*n.rep.genotype.env/3))
  
  names.genotype <- paste0('G',1:n.genotypes)
  
  # first random process
  
  geno.comb.matrix <- lapply(seq(n.envs),function(env){
    randomise(npackages=n.packages, itemnames= names.genotype)
  })
  
  geno.comb.matrix <- do.call(rbind,geno.comb.matrix)
  
  n.farms  <- nrow(geno.comb.matrix)
  n.plots  <- ncol(geno.comb.matrix)
  farm     <- paste0('F',1: n.farms ) # farms will be generated across the enviroments with not repetions
  plot     <- paste0('P',1:n.plots)
  env      <- paste0('Env',1:n.envs) # Note the optimal, enviroment must be stratified
  gen.env  <- paste(names.genotype,env,sep='.')
  
  ##set up data frame with assignment of varieties to farms, plots
  farm.vec <- rep(farm,each= n.plots)
  plot.vec <- rep(plot,length= length(farm.vec))
  env.vec  <- rep(env,array(n.packages*3,n.envs))
  
  
  # set up (small) deviations due to farm and plot (residual in this case)
  
  long.design <- data.frame(env.vec=env.vec,farm.vec = farm.vec, plot.vec = plot.vec , env.farm.vec = paste0(env.vec,'.',farm.vec), farm.plot.vec= paste0(env.vec,'.',farm.vec,'.',plot.vec), genotype = c(t(geno.comb.matrix)))
  
  long.design$gen.env.vec <- factor(paste0(long.design$genotype,'.',long.design$env.vec))
  
  long.design$gen.env.vec <- factor(long.design$gen.env.vec,levels = unique(long.design$gen.env.vec))
  
  long.design$env.vec <- factor(long.design$env.vec,levels = unique(env.vec))
  
  long.design$farm.vec <- factor(long.design$farm.vec,levels = unique(farm.vec)) 
  
  long.design$plot.vec <- factor(long.design$plot.vec,levels = unique(plot.vec))
  
  long.design$env.farm.vec <- factor(long.design$env.farm.vec,levels = unique(long.design$env.farm.vec))
  
  long.design$farm.plot.vec <- factor(long.design$farm.plot.vec,levels = unique(long.design$farm.plot.vec))
  
  long.design$genotype <- factor(long.design$genotype,levels = names.genotype)
  
  # Design matrix
  env.m.design <- if(n.envs!=1){model.matrix(~-1+env.vec,data=long.design)}else{matrix(array(0,nrow(long.design)),ncol=1)}
  env.farm.m.design      <- model.matrix(~-1+env.farm.vec,data=long.design)
  farm.plot.vec.m.design <- model.matrix(~-1+farm.plot.vec,data=long.design)
  genotype.m.design      <- model.matrix(~-1+genotype,data=long.design)
  gen.env.m.design       <- model.matrix(~-1+gen.env.vec,data=long.design)
  
  # Coefficients matrix (second random process)
  
  coef.genotypes <- if(is.null(genotypic.means)){matrix(rnorm(n.genotypes,mean=mu.genotype,sd=sigma.genotype),ncol = 1 , nrow = n.genotypes)}else{matrix(genotypic.means,ncol=1)}
  coef.envs      <- matrix(rnorm(n.envs,mean=0,sd=sigma.env) ,ncol = 1 , nrow =  n.envs)
  coef.gen.env   <- matrix(rnorm(n.genotypes*n.envs,mean=0,sd=sigma.gen.env),ncol = 1,nrow = n.genotypes*n.envs)
  coef.env.farms <- matrix(rnorm(n.farms,mean=0,sd=sigma.env.farm) ,ncol = 1 , nrow =  n.farms) # farms only can be dependend to enviroments
  coef.plots     <- matrix(rnorm(n.plots*n.farms,mean=0,sd=sigma.plot),ncol = 1, nrow = n.plots*n.farms)
  
  
  # Sum of effects
  
  trait.value.obs <- genotype.m.design %*% coef.genotypes +  env.m.design %*% coef.envs+ gen.env.m.design %*% coef.gen.env + env.farm.m.design %*% coef.env.farms + farm.plot.vec.m.design %*% coef.plots
  
  # complete effect
  
  data_simulated <-
    data.frame(long.design,value.genotype = genotype.m.design %*% coef.genotypes,
               value.env = env.m.design %*% coef.envs,
               value.gen.env = gen.env.m.design %*% coef.gen.env,
               value.env.farm = env.farm.m.design %*% coef.env.farms,
               value.env.farm.plot = farm.plot.vec.m.design %*% coef.plots,trait.value.obs = trait.value.obs)
  
  
  
  data_simulated$trait.value.rank <- do.call(c,with(data_simulated,
                                                    tapply(-1*trait.value.obs, farm.vec,rank)))
  
  # Output
  
  pars = data.frame(n.rep.genotype.env=n.rep.genotype.env,n.packages=n.packages,n.genotypes=n.genotypes,mu.genotype=mu.genotype,n.envs =n.envs,
                    sigma.genotype=sigma.genotype,sigma.env.farm=sigma.env.farm,sigma.env=sigma.env,sigma.plot=sigma.plot,sigma.gen.env =sigma.gen.env)
  
  list(dataset=data_simulated,rankings= convert_rankings(data_simulated = data_simulated,names.genotype=names.genotype),pars = pars)
  
}



rank_experiment <- function(w){
  print('------------------------------------------------------')
  print(repetition[w,])
  rept <- repetition[w,]
  
  iter    <- rept$iter
  #sigma_g <- rept$sigma.genotype
  sigma.plot <- rept$sigma.plot
  n.rep.genotype.env <- rept$n.rep.genotype.env
  seed <- rept$iter
 # rep <- rept$rep
  #set.seed(rept$iter)
  
  #iter <- 4
  #sigma_g <- 1.2
  #sigma.plot <- 1.1
  #n.rep.genotype.env <- 200
  
  set.seed(seed)
  
  mu.gs <- data_30_exmple$value.genotype
  
  set.seed(seed)
  sim_dataset <- simulate_multi_env_tricot_V2(n.genotypes=20,mu.genotype =NULL,n.envs =1,
                                              sigma.genotype =NULL,sigma.env.farm=3.352,
                                              sigma.env=0,genotypic.means=mu.gs,
                                              sigma.gen.env=0,sigma.plot=sigma.plot, # big mistake
                                              n.rep.genotype.env=n.rep.genotype.env)
  
  
  #   simulate_multi_env_tricot(n.genotypes=20, #20
  #                                          mu.genotype =3.63, 
  #                                          n.envs =1,
  #                                          sigma.genotype =sigma_g,
  #                                          sigma.env.farm=3.352, # Median maize observed data
  #                                          sigma.env=0,
  #                                          sigma.gen.env=0,
  #                                          sigma.plot=sigma.plot,
  #                                          n.rep.genotype.env=n.rep.genotype.env
  # )
  # 
  #sd(unique(sim_dataset$dataset$value.genotype))
  
  rankings <- sim_dataset$rankings
  
  # set.seed(123)
  
  cp= list(maxit = 10000, temp = 50, trace = TRUE,REPORT = 5000)
  
  matMu <- runif(ncol(rankings),-2,2)
  
  mu.names <- dimnames(rankings)[[2]]
  
  names(matMu) <- paste0('mu',mu.names)
  
  
  start_time <- Sys.time()
  th_mus <- optim(par = matMu[-1], tricotNLLMU0,model="t",data= rankings,
                  method =c("Nelder-Mead","SANN","BFGS")[3],control=cp,hessian=T)
  
  th_time <- Sys.time()-start_time
  
  start_time <- Sys.time()
  log_worths <- optim(par = matMu[-1], tricotNLLMU0,model="pl",data= rankings,
                      method =c("Nelder-Mead","SANN","BFGS")[3],control=cp,hessian=T)
  pl_time <- Sys.time()-start_time
  
  outputs <- sim_dataset
  
  outputs[['th_mus']] <- th_mus
  
  outputs[['log_worths']] <- log_worths
  
  outputs[['iter']] <- iter
  
  #outputs[['rep']] <- rep
  
  
  outputs[['th_time']] <- th_time
  
  outputs[['pl_time']] <- pl_time
  
  outputs[['seed']] <- seed
  
  saveRDS(outputs,paste0('standard_error_val_base_case/output_sd_error_sdg0.8_sd0.9_rep30_',iter,'','.rds'))
  
  outputs
}

outputs <- mclapply(seq(nrow(repetition)), rank_experiment,mc.cores = 10)