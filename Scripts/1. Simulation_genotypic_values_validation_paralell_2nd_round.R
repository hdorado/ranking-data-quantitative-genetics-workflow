
# Genotypic values estimation - simulation
# Hugo Dorado
# created 06/10/2024

library(tidyverse)
library(PlackettLuce)
library(ClimMobTools)
library(mvtnorm)
library(parallel)
library(gtools)

# /home/joost/Hugo_Projects/Sim_est_genotypic_values/

# new_pars_sim
# low_sample_sim

setwd('Paper_1_Estimation_of_genetic_parameters/')

setwd('/home/joost/Hugo_Projects/Sim_est_genotypic_values/')

source('Scripts/0. Ranking_pars_FUN.R')

set.seed(123)

# Low sample size

repetition_lown <- expand.grid(iter=1:30,n.rep.genotype.env=c(3,5),sigma.genotype = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2),sigma.plot = c(0.6, 0.8, 0.9, 1.1, 1.4, 1.6))

repetition_ini_sim <- expand.grid(iter=1:30,n.rep.genotype.env=c(10,30,50,100,200,500),sigma.genotype = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2),sigma.plot = c(0.6, 0.8, 0.9, 1.1, 1.4, 1.6))

# Including new parameters

repetition <- expand.grid(iter=1:30,n.rep.genotype.env=c(3,5,10,30,100),
                          sigma.genotype = c(0.05,0.1,0.2, 0.4, 0.6, 0.8, 1.0, 1.2,1.5,1.8,2),
                          sigma.plot = c(0.4,0.6, 0.8, 0.9, 1.1, 1.4, 1.6,2,2.5,3,3.5,4,4.5,5))

# 30*5*11*14 | 30*5*11*14*20*3
 
repetition <- rbind(repetition_lown,repetition_ini_sim,repetition)

repetition <- repetition %>% dplyr::filter(n.rep.genotype.env %in% c(3,5,10,30,100)) # I did a mistake by having repetition <- repetition %>% dplyr::filter(n.rep.genotype.env == c(3,5,10,30,100)) 

repetition <- repetition[!duplicated(repetition),]

repetition <- repetition %>% mutate(file = paste0('output_',iter,'-',sigma.genotype,'-',sigma.plot,'-',n.rep.genotype.env,'.rds'))

gvp <- list.files('Genotypic_values_performance3')

table(repetition$file %in% gvp) #23100 files in total

repetition <- repetition[!(repetition$file %in% gvp),]

nrow(repetition)

#repetition <- repetition[1:8376,]



repetition <- repetition[,-5]

# 5*11*14*30

# sort(unique(c(0.8,0.2,1.2,0.4, 0.6, 1.0 )))
# sort(unique(c(0.9,0.6,1.6,0.8,1.1,1.4)))

# 22,26,28  # 0.4-0.6
# 2,3,5,8 #0.6-0.6
# 25 # 0.6-1.6


#repetition$seed <- 8001:(8000+nrow(repetition)) # seeds updated
repetition$seed <- 17001:(17000+nrow(repetition)) # seeds updated

# w <- 2514

rank_experiment <- function(w){
  print('------------------------------------------------------')
  print(repetition[w,])
  rept <- repetition[w,]
  
  iter    <- rept$iter
  sigma_g <- rept$sigma.genotype
  sigma.plot <- rept$sigma.plot
  n.rep.genotype.env <- rept$n.rep.genotype.env
  seed <- rept$seed
  #set.seed(rept$iter)
  
  #iter <- 4
  #sigma_g <- 1.2
  #sigma.plot <- 1.1
  #n.rep.genotype.env <- 200
  
  set.seed(seed)
  sim_dataset <- simulate_multi_env_tricot(n.genotypes=20, #20
                                           mu.genotype =3.63, 
                                           n.envs =1,
                                           sigma.genotype =sigma_g,
                                           sigma.env.farm=3.352, # Median maize observed data
                                           sigma.env=0,
                                           sigma.gen.env=0,
                                           sigma.plot=sigma.plot,
                                           n.rep.genotype.env=n.rep.genotype.env
  )
  
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
  
  outputs[['rep']] <- rep
  
  
  outputs[['th_time']] <- th_time
  
  outputs[['pl_time']] <- pl_time
  
  outputs[['seed']] <- seed
  
  saveRDS(outputs,paste0('Genotypic_values_performance3/output_',iter,'-',sigma_g,'-',sigma.plot,'-',n.rep.genotype.env,'.rds'))
  
  outputs
}

outputs <- mclapply(seq(nrow(repetition)), rank_experiment,mc.cores = 35)




