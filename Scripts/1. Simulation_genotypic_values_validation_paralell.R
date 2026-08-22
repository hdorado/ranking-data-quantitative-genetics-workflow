
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

setwd('Paper_1_Estimation_of_genetic_parameters/')

source('Scripts/0. Ranking_pars_FUN.R')

set.seed(123)



repetition <- expand.grid(iter=1:30,n.rep.genotype.env=c(10,30,50,100,200,500),sigma.genotype = c(0.8,0.2,1.2),sigma.plot = c(0.9,0.6,1.6) )

repetition <- expand.grid(iter=1:30,n.rep.genotype.env=c(10,30,50,100,200,500),sigma.genotype = c(0.4, 0.6, 1.0 ),sigma.plot = c(0.8,1.1,1.4) )

repetition <- expand.grid(iter=1:30,n.rep.genotype.env=c(10,30,50,100,200,500),sigma.genotype = c(0.8,0.2,1.2),sigma.plot = c(0.8,1.1,1.4) )

repetition <- expand.grid(iter=1:30,n.rep.genotype.env=c(10,30,50,100,200,500),sigma.genotype = c(0.4, 0.6, 1.0 ),sigma.plot = c(0.9,0.6,1.6) )

repetition <- repetition[-c(1:1080),]

repetition <- expand.grid(iter=1:30,n.rep.genotype.env=c(10,30,50,100,200,500),sigma.genotype = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2),sigma.plot = c(0.6, 0.8, 0.9, 1.1, 1.4, 1.6))

# sort(unique(c(0.8,0.2,1.2,0.4, 0.6, 1.0 )))
# sort(unique(c(0.9,0.6,1.6,0.8,1.1,1.4)))

# 22,26,28  # 0.4-0.6
# 2,3,5,8 #0.6-0.6
# 25 # 0.6-1.6

repetition5 <- rbind(expand.grid(iter=c(22,26,28),n.rep.genotype.env=c(500),sigma.genotype = c(0.4),sigma.plot = c(0.6) ),
                    expand.grid(iter=c(2,3,5,8),n.rep.genotype.env=c(500),sigma.genotype = c(0.6),sigma.plot = c(0.6) ),
                    expand.grid(iter=c(25),n.rep.genotype.env=c(500),sigma.genotype = c(0.6),sigma.plot = c(1.6) ))

dim(all_possible_iter)

# Experiment sigma.genotype = c(0.4, 0.6, 1.0 ),sigma.plot = c(0.9,0.6,1.6)

# sever 1

repetition1 <- repetition[repetition$n.rep.genotype.env %in% c(10,30,50,100),] # Completed at the beast (exp 4 missing 100)

# server 2
  
repetition2 <- repetition[repetition$n.rep.genotype.env ==  200,] # Sent to the beast (exp4) #probably duplicate ouput4 200 exp 3
 
# server 3

repetition3 <- repetition[repetition$n.rep.genotype.env == 500 ,][129:149,] # (129>150 incompleto) Sent to linux 92 (exp 4) , before(150:270) (missing repetitions sent to the beast)

repetition3 <- repetition3[21:90,]


repetition$seed <- 1:nrow(repetition)

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
  
  saveRDS(outputs,paste0('Genotypic_values_performance2/output_',iter,'-',sigma_g,'-',sigma.plot,'-',n.rep.genotype.env,'.rds'))
  
  outputs
}

outputs <- mclapply(seq(nrow(repetition)), rank_experiment,mc.cores = 40)




