
# Hugo Dorado
# Estimation of genetic parameters based on rankings and metric data
# Created: 06-11/2024

library(tidyverse)
library(PlackettLuce)
library(ClimMobTools)
library(mvtnorm)
library(gosset)

setwd('Paper_1_Estimation_of_genetic_parameters')

source('Scripts/0. Ranking_pars_FUN.R')

#------------------------------------------------------------------------------
#------------------------------ Cassava ---------------------------------------

set.seed(123)

cassava_data <- read.csv('Processed_data/on-farm-data-cassava-nigeria.csv')

zone_state <- unique(read.csv('Processed_data/cassava_metric_dataset_post_AE.csv')[c("state","Zone")])

cassava_data <- cassava_data %>% left_join(zone_state,by=c('state'))

## Rankings analysis

# full dataset fresh yield------------------------------------------------------


cassava_data_tricot <- cassava_data[which(!is.na(cassava_data$rank_freshrootsyield) & cassava_data$rank_freshrootsyield %in% c(1,2,3)& !is.na(cassava_data$yieldton_ha) ),]

##remove farms with missing plots

cassava_data_tricot$farm <- as.character(cassava_data_tricot$farm)

tap<-with(cassava_data_tricot,tapply(plot, farm,function(x) length(unique(x))))
tap<-tap[tap==3]
cassava_data_tricot <- cassava_data_tricot[cassava_data_tricot$farm %in% names(tap),]


cassava_data_tricot$genotype <- as.character(cassava_data_tricot$genotype)

cassava_data_tricot <- cassava_data_tricot %>% mutate(ID=as.integer(as.numeric(factor(farm))))

ls_cassava_data_tricot <- split(cassava_data_tricot,paste(cassava_data_tricot$Zone,sep="_"))

log_worth_location <- lapply(ls_cassava_data_tricot,function(x){
  cassava_data_tricot2 <- x[c("farm","genotype","yieldton_ha")]
  
  length(unique(cassava_data_tricot2$farm))
  
  cassava_data_tricot2$yieldton_ha
  
  cassava_data_tricot2 <- cassava_data_tricot2 %>% group_by(farm) %>% group_split()  %>%  
    map( ~ .x %>%  mutate(order = base::rank(-1*yieldton_ha,ties.method = 'random'))) %>% bind_rows
  
  
  tricot.rank.data <-rank_numeric(data = cassava_data_tricot2[,-3],
                                  items = "genotype",
                                  input = "order", 
                                  id = "farm",ascending = FALSE)
  rankings <- tricot.rank.data
  
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
  
  outputs <- list()
  
  outputs[['rankings']] <- rankings
  
  outputs[['dataset']] <- x
  
  outputs[['th_mus']] <- th_mus
  
  outputs[['log_worths']] <- log_worths
  
  outputs[['th_time']] <- th_time
  
  outputs[['pl_time']] <- pl_time
  
  outputs[['n']] <- nrow(x)
  
  outputs[['loc']] <- unique(x$Zone)
  
  outputs[['ref']] <- names(matMu)[[1]]
  
  saveRDS(outputs,paste0('Observed_data_outputs/output_cassava_',unique(x$Zone),'_yield_ranked.rds'))
  
}
)



# set.seed(123)



#------------------------------------------------------------------------------
#----------------------------- Groundnut --------------------------------------

set.seed(123)

groundnut_data <- read.csv("Processed_data/on-farm-data-groundnut-tanzania.csv")

# overall

groundnut_data_yld <- groundnut_data[which(!is.na(groundnut_data$rank_yield) & groundnut_data$rank_yield %in% c(1,2,3)& !is.na(groundnut_data$yieldton_ha)) ,]

groundnut_data_yld <- groundnut_data_yld[!is.na(groundnut_data_yld$ecoregion),]

groundnut_data_yld$farm <- as.character(groundnut_data_yld$farm)

tap<-with(groundnut_data_yld,tapply(plot, farm,function(x) length(unique(x))))

tap<-tap[tap==3]

groundnut_data_yld <- groundnut_data_yld[groundnut_data_yld$farm %in% names(tap),]

groundnut_data_yld$genotype <- as.character(groundnut_data_yld$genotype)

groundnut_data_yld_tricot <- groundnut_data_yld %>% mutate(ID=as.integer(as.numeric(factor(farm))))

ls_groundnut_data_yld_tricot <- split(groundnut_data_yld_tricot,paste0(groundnut_data_yld_tricot$ecoregion,"-",groundnut_data_yld_tricot$year))

#ls_groundnut_data_yld_tricot <- ls_groundnut_data_yld_tricot[-c(2,3,6)]

# Lake-2021, `Western-2020`, `Western-2021`,`Southern-2020`, Southern-2021`, `Southern highlands-2023`

log_worth_location <- lapply(ls_groundnut_data_yld_tricot[c('Central-2020', 'Central-2022','Southern-2023','Southern highlands-2022')],function(x){
  print(unique(paste(x$ecoregion,x$year)))
  
  groundnut_data_yld_tricot2 <- x[c("farm","genotype","yieldton_ha")]
  
  groundnut_data_yld_tricot2 <- groundnut_data_yld_tricot2 %>% group_by(farm) %>% group_split()  %>%  
    map( ~ .x %>%  mutate(order = base::rank(-1*yieldton_ha,ties.method = 'random'))) %>% bind_rows
  
  
  tricot.rank.data <-rank_numeric(data = groundnut_data_yld_tricot2[,-3],
                                  items = "genotype",
                                  input = "order", 
                                  id = "farm",ascending = FALSE)
  rankings <- tricot.rank.data
  
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
  
  outputs <- list()
  
  outputs[['rankings']] <- rankings
  
  outputs[['dataset']] <- x
  
  outputs[['th_mus']] <- th_mus
  
  outputs[['log_worths']] <- log_worths
  
  outputs[['th_time']] <- th_time
  
  outputs[['pl_time']] <- pl_time
  
  outputs[['n']] <- nrow(x)
  
  outputs[['loc']] <- unique(x$Zone)
  
  outputs[['ref']] <- names(matMu)[[1]]
  
  saveRDS(outputs,paste0('Observed_data_outputs/output_groundnut_',unique(paste(x$ecoregion,x$year)),'_order.rds'))
}
)


#------------------------------------------------------------------------------
#------------------------------- Maize ----------------------------------------

# what should we do to compare genotypes tested with different maturity

set.seed(123)

maize_data <- read.csv('Processed_data/on-farm-data-maize-kenia.csv')

maize_gen_info <- readxl::read_xlsx('Processed_data/genotype-information-maize-kenya_GMK_HD.xlsx')

maize_data <- maize_data %>% left_join(maize_gen_info[c('Genotype','Maturity')],by=c('genotype'='Genotype'))

## Rankings analysis

# full dataset fresh yield------------------------------------------------------

maize_data_tricot <- maize_data[which(!is.na(maize_data$rank_yield) & maize_data$rank_yield %in% c(1,2,3) & maize_data$yield_farmer ),]

summary(maize_data_tricot$yield_farmer)

##remove farms with missing plots

maize_data_tricot$farm <- as.character(maize_data_tricot$farm)

tap<-with(maize_data_tricot,tapply(plot, farm,function(x) length(unique(x))))
tap<-tap[tap==3]
maize_data_tricot <- maize_data_tricot[maize_data_tricot$farm %in% names(tap),]


maize_data_tricot$genotype <- as.character(maize_data_tricot$genotype)

maize_data_tricot <- maize_data_tricot %>% mutate(ID=as.integer(as.numeric(factor(farm))))

ls_maize_data_tricot <- split(maize_data_tricot,paste(maize_data_tricot$zone,sep="_"))

log_worth_location <- lapply(ls_maize_data_tricot,function(x){
  
  maize_data_tricot2 <- x[c("farm","genotype","yield_farmer")]
  
  length(unique(maize_data_tricot2$farm))
  
  maize_data_tricot2 <- maize_data_tricot2 %>% group_by(farm) %>% group_split()  %>%  
    map( ~ .x %>%  mutate(order = base::rank(-1*yield_farmer,ties.method = 'random'))) %>% bind_rows
  
  tricot.rank.data <-rank_numeric(data = maize_data_tricot2[-3],
                                  items = "genotype",
                                  input = "order", 
                                  id = "farm",ascending = FALSE)
  rankings <- tricot.rank.data
  
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
  
  outputs <- list()
  
  outputs[['rankings']] <- rankings
  
  outputs[['dataset']] <- x
  
  outputs[['th_mus']] <- th_mus
  
  outputs[['log_worths']] <- log_worths
  
  outputs[['th_time']] <- th_time
  
  outputs[['pl_time']] <- pl_time
  
  outputs[['n']] <- nrow(x)
  
  outputs[['loc']] <- unique(x$Zone)
  
  outputs[['ref']] <- names(matMu)[[1]]
  
  saveRDS(outputs,paste0('Observed_data_outputs/output_maize_',unique(x$zone),'_order.rds'))
  
}
)

