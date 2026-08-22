
# Data aggregation to process the resuls of the simulation
# Hugo Dorado
# 22-07-2024

library(tidyverse)
library(PlackettLuce)
library(asreml)
library(lme4)
library(nlme)
library(MASS)

setwd('/home/joost/Hugo_Projects/Sim_est_genotypic_values/')

# Read final file cosolidated---------------------------------------------------

ls_estimations <- readRDS('Processed_data/output_ls_estimations2.rds')

#------------------------------------------------------------------------------
# 
# Incorporate the results of metric data

summary.asrmlErr <- function(object, ...) {
  df <- data.frame(component=c(NA,NA))
  row.names(df) <- c('genotype','R')  
  return(list(varcomp=df))
}

# Calulate predictor error variances


PEV <- function(dataset,genotype='genotype',means='est_values',invw='quasiSE'){
  require(asreml)
  ds <- dataset[c(genotype,means,invw)]
  names(ds) <- c('genotype','est_means','quasiSE')
  ds$genotype <- factor(ds$genotype )
  
  
  if( !(sum(is.na(ds['quasiSE']))>0 | sum(ds['quasiSE']<=0)>0) ){ # this can be wrong in the rest of code check!
    ds$w2 <- 1/ds[,invw]^2 # Quasi-standard errors blues
    
    modVal <- try(asreml(fixed = est_means ~ 1, random = ~ genotype , weights = w2,
                         family = (asr_gaussian(dispersion = 1)),
                         data = ds,trace=F),TRUE)
    
    if("try-error" %in% class(modVal)){
      #a <- 'err'
      #class(a) <- 'asrmlErr'
      #a
      return(array(NA,nrow(ds)))
    }else{
      mod <- asreml(fixed = est_means ~ 1, random = ~ genotype , weights = w2,
             family = (asr_gaussian(dispersion = 1)),
             data = ds,trace=F)
      BLUPl <- summary(mod, coef=TRUE)$coef.random
      BLUPl <-as.data.frame(BLUPl)
      BLUPl <- BLUPl[grep('genotype_',rownames(BLUPl)),]
      
      return(BLUPl$std.error^2)
      
    }

  }else{array(NA,nrow(ds))}
  
}






#



genotypic_sd_2stages <- function(dataset,weights_qv='quasiSE',weights_cont='SE_cont',weights_smZ='SE_smZ',weights_SE='SE',ref="G1"){
  
  # Heritability function
  
  qH2 <- function(model){
    if(class(model)=='asrmlErr'){return(NA)}
    
    varcomp <- summary(model)$varcomp
    
    BLUPl <- summary(model, coef=TRUE)$coef.random
    BLUPl <-as.data.frame(BLUPl)
    BLUPl <- BLUPl[grep('genotype_',rownames(BLUPl)),]
    
    BLUPl$PEV <- BLUPl$std.error^2
    
    PEV1<-BLUPl$PEV
    
    1-mean(PEV1)/varcomp['genotype','component']
  }
  
  # Estimation of genetic parameters
  
  require(asreml)
  
  dataset$genotype <- factor(dataset$genotype)

  # Valid standard errors/quasi-standard errors
 
  results2stg <- data.frame(model=unique(dataset$model),
             sd_g_w.nW = NA,
             sd_g_w.qV= NA,
             sd_g_w.Cotr = NA,
             sd_g_w.GM = NA,
             sd_g_w.SumZr = NA,
             sd_g_sc.scl = NA,
             
             qh2.2st.nW = NA, 
             qh2.2st.w.qV = NA, 
             qh2.2st.Cotr= NA,
             qh2.2st.w.GM =NA,
             qh2.2st.w.SumZr =NA,
             qh2.2st.scl=NA)
  

  # No weighting
  
  modVal <- try(asreml(fixed = est_means ~ 1, random = ~ genotype , weights = 1,
                       family = (asr_gaussian(dispersion = 1)),
                       data = dataset,trace=F),TRUE)
  
  
  stage2.nW <- if("try-error" %in% class(modVal)){
    a <- 'err'
    class(a) <- 'asrmlErr'
    a
  }else{
    asreml(fixed = est_means ~ 1, random = ~ genotype , weights = 1,
           family = (asr_gaussian(dispersion = 1)),
           data = dataset,trace=F)
  }


  results2stg$sd_g_w.nW = sqrt(summary(stage2.nW)$varcomp['genotype','component'])
  
  results2stg$qh2.2st.nW = qH2(stage2.nW)
  
  
  # Quasi-standard errors
  
  if( !(sum(is.na(dataset$quasiSE))>0 | sum(dataset$quasiSE<=0)>0) ){
    dataset$w2 <- 1/dataset[,weights_qv]^2 # Quasi-standard errors blues
    
    modVal <- try(asreml(fixed = est_means ~ 1, random = ~ genotype , weights = w2,
                         family = (asr_gaussian(dispersion = 1)),
                         data = dataset,trace=F),TRUE)
    
    stage2.qV <- if("try-error" %in% class(modVal)){
      a <- 'err'
      class(a) <- 'asrmlErr'
      a
    }else{
      asreml(fixed = est_means ~ 1, random = ~ genotype , weights = w2,
             family = (asr_gaussian(dispersion = 1)),
             data = dataset,trace=F)
    }
    
    results2stg$sd_g_w.qV= sqrt(summary(stage2.qV)$varcomp['genotype','component'])
    
    results2stg$qh2.2st.w.qV = qH2(stage2.qV)
  }  
  

  
  # Genotypic means
  
  if( !( sum(is.na(dataset$SE))>0 |  sum(dataset$SE<=0)>0 ) ){
  
  dataset$w <- 1/dataset[,weights_SE]^2 # Genotypic means
  
  stage2.GM <- if(unique(dataset$model)!='MM'){
    a <- 'err'
    class(a) <- 'asrmlErr'
    a
    
  }else{
     modVal <-  try(asreml(fixed = est_means ~ 1, random = ~ genotype , weights =w,
                               family = (asr_gaussian(dispersion = 1)),
                               data = dataset,trace=F),TRUE) # 
    
    
     if("try-error" %in% class(modVal)){
      a <- 'err'
      class(a) <- 'asrmlErr'
      a
    }else{    asreml(fixed = est_means ~ 1, random = ~ genotype , weights =w,
                     family = (asr_gaussian(dispersion = 1)),
                     data = dataset,trace=F) # 
          }
  }
  
  results2stg$sd_g_w.GM = sqrt(summary(stage2.GM)$varcomp['genotype','component'])
  
  results2stg$qh2.2st.w.GM =qH2(stage2.GM)
  
  }
  
  # Standard error
  
  if( !( sum(is.na(dataset$SE_cont))>0 | sum(dataset$SE_cont==0)>1 | sum(dataset$SE_cont<0)>0 )){  
  
  dataset$w4 <- 1/dataset[,weights_cont]^2 # Contrasts
  
  if(unique(dataset$model)=='MM'){
    dataset[,'est_cont'][1] <- 0
  }
  
  dataset[,'w4'][1] <- 1 # weight equal to 1, probably doesn't make sense
  
  modVal <- try(asreml(fixed = est_cont ~ 1, random = ~ genotype , weights =w4,
                                family = (asr_gaussian(dispersion = 1)),
                                data = dataset,trace=F),TRUE) # 
  
  stage2.Cotrs <- if("try-error" %in% class(modVal)){
    a <- 'err'
    class(a) <- 'asrmlErr'
    a
    }else{    asreml(fixed = est_cont ~ 1, random = ~ genotype , weights =w4,
                 family = (asr_gaussian(dispersion = 1)),
                 data = dataset,trace=F) # 
    }
  
  results2stg$sd_g_w.Cotr = sqrt(summary(stage2.Cotrs)$varcomp['genotype','component'])
  
  results2stg$qh2.2st.Cotr= qH2(stage2.Cotrs)
  }
  
  # SumZero constraint
  
  if( !( sum(is.na(dataset$SE_smZ))>0 |  sum(dataset$SE_smZ<=0)>1 ) ){
  
  dataset$w3 <- 1/dataset[,weights_smZ]^2 # SumZero only rankings
  
  stage2.SumZr <- if(unique(dataset$model)=='MM'){
    a <- 'err'
    class(a) <- 'asrmlErr'
    a

  }else{
    
  modVal <- try(asreml(fixed = est_smZ ~ 1, random = ~ genotype , weights =w3,
                                family = (asr_gaussian(dispersion = 1)),
                                data = dataset,trace=F),TRUE) # 
  
  
    if("try-error" %in% class(modVal)){
      a <- 'err'
      class(a) <- 'asrmlErr'
      a
    }else{    asreml(fixed = est_smZ ~ 1, random = ~ genotype , weights =w3,
                   family = (asr_gaussian(dispersion = 1)),
                   data = dataset,trace=F) # 
    }
  }
  
  results2stg$sd_g_w.SumZr = sqrt(summary(stage2.SumZr)$varcomp['genotype','component'])
  
  results2stg$qh2.2st.w.SumZr =qH2(stage2.SumZr)
  }
  
  
  # Scaled

  
  modVal <- try(asreml(fixed = scl_means ~ 1, random = ~ genotype , weights = 1,
                       family = (asr_gaussian(dispersion = 1)),
                       data = dataset,trace=F),TRUE)
  
  
  stage2.scl <- if("try-error" %in% class(modVal)){
    a <- 'err'
    class(a) <- 'asrmlErr'
    a
  }else{
    asreml(fixed = scl_means ~ 1, random = ~ genotype , weights = 1,
           family = (asr_gaussian(dispersion = 1)),
           data = dataset,trace=F)
  }
  
  results2stg$qh2.2st.scl=qH2(stage2.scl)
  
  results2stg
  
  
  

}

#---------------------------Genotypic variance estimation---------------------

qH2 <- function(model,ref='genotype_'){
  varcomp <- summary(model)$varcomp
  
  BLUPl <- summary(model, coef=TRUE)$coef.random
  BLUPl <-as.data.frame(BLUPl)
  BLUPl <- BLUPl[grep(ref,rownames(BLUPl)),]
  
  BLUPl$PEV <- BLUPl$std.error^2
  
  PEV1<-BLUPl$PEV
  
  1-mean(PEV1)/varcomp['genotype','component']
}

#out <- ls_estimations[['output_7-0.8-0.9-100.rds']]

consolidated <- NULL

summary_analysis_iter_g_var <- for(i in 1:length(ls_estimations) ){ #seq(length(ls_estimations))
  
  print(i)
  
  out <- ls_estimations[[i]]
  
  #out <-ls_estimations[[which(names(ls_estimations)=='output_8-0.8-0.9-10.rds')]]
  
  ng <- out$pars$n.rep.genotype.env
  
  sigma_e <- out$pars$sigma.plot
  
  sigma_g <- out$pars$sigma.genotype
  
  
  
  #if(round(det(th_hess),3)==0){print('no_sing');return(NULL)}
  
  # --------------------- Thurstonian ---------------------
  
  th_hess <- out$th_mus$hessian
  
  Th_estimates = c(muG1=0,out$th_mus$par)
  
  th_inv_hess <- try(solve(th_hess),TRUE)
  
  sing_th <- "try-error" %in% class(th_inv_hess)
  
  # Calculate the g-inverse
  
  if(sing_th){
    th_inv_hess <- ginv(th_hess)
    colnames(th_inv_hess) <- colnames(th_hess) 
    rownames(th_inv_hess) <- colnames(th_hess)
  }
  
  # Create cov matrix
  
  vcov_th <- cbind(muG1=0,rbind(muG1=0,th_inv_hess))
  
  # Caculate qErrors
  
  quasiVar_th <- try( qvcalc(vcov_th),TRUE)
  
  # sumzero means
  
  Th_worthval <- matrix(Th_estimates,ncol=1)
  
  k <- nrow(Th_worthval)
  
  ma <- matrix(-1,k,k)
  
  C <- ma + diag(k)*k 
  
  Th_estimates_sum_zero <- 1/k*C %*% Th_worthval
  
  # sumzero se
  
  Th_sumZerocov <- 1/k^2*C %*% vcov_th %*% t(C)
  
  Th_estimates_se_table <- data.frame(estimate=Th_estimates,
                                   SE=sqrt(diag(vcov_th)), # Can give negative values
                                   est_smZ =Th_estimates_sum_zero,
                                   SE_smZ =sqrt(diag( Th_sumZerocov ))
                                   )
  
  Th_estimates_se_table$quasiSE <- if("try-error" %in% class(quasiVar_th)){
    NA
    }else{quasiVar_th$qvframe[,'quasiSE']}
    


  # --------------------- Plackett-Luce --------------------- 
  

  log_worths_hess <- out$log_worths$hessian
  
  PL_estimates = c(muG1=0,out$log_worths$par)
  
  pl_inv_hess <- try(solve(log_worths_hess),TRUE)
  
  sing_pl <- "try-error" %in% class(pl_inv_hess)
  
  # Calculate the g-inverse
  
  if(sing_pl){
    pl_inv_hess <- ginv(log_worths_hess)
    colnames(pl_inv_hess) <- colnames(log_worths_hess) 
    rownames(pl_inv_hess) <- colnames(log_worths_hess)
  }

  
  # Create cov matrix
  vcov_pl <- cbind(muG1=0,rbind(muG1=0,pl_inv_hess))
  
  # Caculate qErrors
  
  quasiVar_pl <- try(qvcalc(vcov_pl),TRUE)
  
  # Sumzer means
  
  PL_worthval <- matrix(PL_estimates,ncol=1)
  
  k <- nrow(PL_worthval)
  
  ma <- matrix(-1,k,k)
  
  C <- ma + diag(k)*k 
  
  PL_estimates_sum_zero <- 1/k*C %*% PL_worthval
  
  # sumzero se
  
  PL_sumZerocov <- 1/k^2*C %*% vcov_pl %*% t(C)
  
  
  PL_estimates_se_table <- data.frame(estimate=PL_estimates,
                                      SE=sqrt(diag(vcov_pl)),
                                      est_smZ =  PL_estimates_sum_zero,
                                      SE_smZ=sqrt(diag( PL_sumZerocov )))  
  
  PL_estimates_se_table$quasiSE <- if("try-error" %in% class(quasiVar_pl)){
    NA
  }else{quasiVar_pl$qvframe[,'quasiSE']}

  
  # -------------------Placket-Luce + Thurstonian-----------------------
  # cont means constrast, for ranking models will be always the same

  Th_estimates_se_table$genotype <- row.names(Th_estimates_se_table)
  PL_estimates_se_table$genotype<- row.names(PL_estimates_se_table)
  
  rank_analysis <- rbind(data.frame(genotype = row.names(Th_estimates_se_table),
                                    est_means=Th_estimates_se_table$estimate,
                                    est_cont=Th_estimates_se_table$estimate,
                                    est_smZ = Th_estimates_se_table$est_smZ,
                                    SE=Th_estimates_se_table$SE,
                                    SE_cont = Th_estimates_se_table$SE,
                                    SE_smZ = Th_estimates_se_table$SE_smZ,
                                    quasiSE = Th_estimates_se_table$quasiSE,
                                    PEV= PEV(dataset = Th_estimates_se_table,means = 'estimate'),
                                    model='Th',Singular=sing_th),
                         data.frame(genotype = row.names(PL_estimates_se_table),
                                    est_means=PL_estimates_se_table$estimate,
                                    est_cont=PL_estimates_se_table$estimate,
                                    est_smZ = PL_estimates_se_table$est_smZ,
                                    SE=PL_estimates_se_table$SE,
                                    SE_cont = PL_estimates_se_table$SE,
                                    SE_smZ = PL_estimates_se_table$SE_smZ,
                                    quasiSE = PL_estimates_se_table$quasiSE,
                                    PEV= PEV(PL_estimates_se_table,means = 'estimate'),
                                    model='PL',Singular=sing_pl
                                    ))
  
  rank_analysis$genotype <- gsub('^mu','',rank_analysis$genotype)
  
  truevals <-  unique(out$dataset[c('genotype','value.genotype')])
  
  rank_analysis <- rank_analysis %>% left_join(truevals,by='genotype')
  
  final_ds <- data.frame( rank_analysis)
  
  # Metric data
  
  #metric_analysis <- data.frame(data.frame(emmeans(out$mod.blues,'genotype'))[c('genotype','emmean',"SE")],model='MM')
  
  # Quasi-standar errors metric data
  
  blues <- fixef(out$mod.blues)
  
  sum_blues <- summary(out$mod.blues)
  
  con_coef <- sum_blues$coefficients
  
  L <- diag(1,nrow = length(blues))  # Contrast matrix
  
  L[,1] <- 1
  
  VCOV <- L %*% as.matrix(vcov(out$mod.blues)) %*% t(L)
  
  bluesNams <- gsub('genotype','',names(blues))
  
  bluesNams[1] <- "G1"
  
  qSErrors <- qvcalc(VCOV,labels =  bluesNams,estimate = L%*%blues)
  
  qSErrors <- as.data.frame(qSErrors$qvframe)
  
  qSErrors <- data.frame(genotype=row.names(qSErrors),qSErrors)
  

  metric_analysis <- data.frame(
    genotype = qSErrors$genotype,
    est_means=qSErrors$estimate[,1],
    est_cont=blues,
    est_smZ = NA,
    SE=qSErrors$SE,
    SE_cont = con_coef[,2],
    SE_smZ = NA,
    quasiSE = qSErrors$quasiSE,
    PEV= PEV(qSErrors,means = 'estimate'),
    model='MM',Singular=FALSE,scl_means=qSErrors$estimate[,1])
    
    
  names(metric_analysis) <- str_replace_all( names(metric_analysis), c('estimate'='est_means'))
  
  metric_analysis <- metric_analysis %>% left_join(truevals,by='genotype')
  
  refvar <- metric_analysis[metric_analysis$genotype=='G1','est_means'][1]
  
  final_ds <- final_ds %>% mutate(scl_means=refvar+est_means*sigma_e) 
  
  full_dataset <- rbind(final_ds,metric_analysis)
  
  #----------------------Variance components estimation------------------------
  
  # metric reference
  
  var_est <- data.frame(summary(out$mod.blups)$varcor)
  
  h2_broad_metric1 <- var_est[var_est$grp=='genotype','vcov']/(var_est[var_est$grp=='genotype','vcov']+var_est[var_est$grp=='Residual','vcov'])
  
  h2_broad_metric2 <- var_est[var_est$grp=='genotype','vcov']/(var_est[var_est$grp=='genotype','vcov']+var_est[var_est$grp=='Residual','vcov']/unique(out$pars$n.rep.genotype.env))
  
  # summary(lmer(trait.value.obs~1+(1|genotype)+(1|env.farm.vec),data=out$dataset))

  arSmodel <- asreml(fixed = trait.value.obs~1,random = ~genotype +env.farm.vec,
                     data=out$dataset,trace=F)
  
  h2_metric_cullis <- qH2(arSmodel)
  print(i)
  h2 <- do.call(rbind,lapply(split(full_dataset,full_dataset$model) ,genotypic_sd_2stages,ref="G1"))
  
  consolidated0 <- data.frame(iter=out$iter,ng=ng,sigma_e=sigma_e,sigma_g=sigma_g,estsigma_g_metric =  sqrt(var_est[var_est$grp=='genotype','vcov']),
                            estsigma_e_metric =  sqrt(var_est[var_est$grp=='Residual','vcov']),  
                            h2_broad_metric1 = h2_broad_metric1,
                            h2_broad_metric2 = h2_broad_metric2,h2_metric_cullis,full_dataset)
  
  
  consolidated0 <- consolidated0 %>% left_join(h2,'model')
  
  consolidated <- rbind(consolidated,consolidated0)
  
  write.csv(consolidated,'Results/summary_analysis_iter_g_var_table_V8.csv', #v3 second round v4 four round where I included more sample size
            row.names = F)
  
}

summary_analysis_iter_g_var_table <- do.call(rbind,summary_analysis_iter_g_var)

table(summary_analysis_iter_g_var_table$sing_pl,summary_analysis_iter_g_var_table$sing_pl)




