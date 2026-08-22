
library(PlackettLuce)
library(tidyverse)
library(MASS)

source("Scripts/0. Ranking_pars_FUN.R", echo=TRUE)

set.seed(123)

tricot_simulation <- simulate_multi_env_tricot(n.genotypes=10,mu.genotype =5,n.envs =1,
                          sigma.genotype =0.91,sigma.env.farm=4.39,sigma.env=0,
                          sigma.gen.env=0,sigma.plot=0.8,n.rep.genotype.env=10)


rankings <- tricot_simulation$rankings


ranking_model <- function(rankings,model=c("t",'pl')[1]){

  outputs <- list()
  cp= list(maxit = 10000, temp = 50, trace = TRUE,REPORT = 5000)
  
  matMu <- runif(ncol(rankings),-2,2)
  
  mu.names <- dimnames(rankings)[[2]]
  
  names(matMu) <- paste0('mu',mu.names)
  
  
  start_time <- Sys.time()
  mus <- optim(par = matMu[-1], tricotNLLMU0,model=model,data= rankings,
                  method =c("Nelder-Mead","SANN","BFGS")[3],control=cp,hessian=T)
  
  time <- Sys.time()-start_time
  
  outputs[['rank_means']] <- mus
  
  outputs[['time']] <- time
  
  outputs$par <- c(muG1=0,mus$par)
  
  outputs
  
}

QuasiVarEstimation <- function(ranks_outpts){
  
  out <- list()  
  
  mu_hess <- ranks_outpts$rank_means$hessian
  
  mu_inv_hess <- try(solve(mu_hess),TRUE)
  
  sing_th <- "try-error" %in% class(mu_inv_hess)
  
  if(sing_th){
    th_inv_hess <- ginv(mu_hess)
    colnames(mu_inv_hess) <- colnames(mu_hess) 
    rownames(mu_inv_hess) <- colnames(mu_hess)
    mu_inv_hess
  }
  
  
  vcov_mu <- cbind(muG1=0,rbind(muG1=0,mu_inv_hess))
  
  quasiVar_mu <- try(qvcalc(vcov_mu,estimates = c(muG0=0,ranks_outpts$rank_means$par)),TRUE)
  
  if("try-error" %in% class(quasiVar_mu)){
    quasiVar_mu <- data.frame(estimates = c(muG0=0,ranks_outpts$rank_means$par),SE =NA,    quasiSE =NA,    quasiVar=NA)
  }else{
    quasiVar_mu <- quasiVar_mu$qvframe
  }
  
  rank_analysis <- rbind(data.frame(genotype = row.names(quasiVar_mu),
                                    est_means=quasiVar_mu$estimate,
                                    SE=quasiVar_mu$SE,
                                    quasiSE = quasiVar_mu$quasiSE,
                                    quasiVar=quasiVar_mu$quasiVar))
  
  rank_analysis$genotype <- gsub('^mu','',rank_analysis$genotype)
  
  rank_analysis
  
}

ranks_outpts <- ranking_model(rankings,model='th')

QuasiVarEstimation(ranks_outpts)

