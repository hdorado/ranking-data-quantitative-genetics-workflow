
# Data analysis observed datasets
# 19-10-2025
# Hugo Dorado

library(tidyverse)
library(readxl)
library(gosset)
library(PlackettLuce)
library(predictmeans)
library(asreml)

source("C:/Users/dorad002/R_Projects/GxExY_two_steps_ranking_data/Scripts2/master_funs.R")
setwd('Paper_1_Estimation_of_genetic_parameters')

# setwd('/home/joost/Hugo_Projects/Observed_data_reanalysis/')

# Cassava

dataset_cassava <- read.csv('Processed_data/On-farm-datasets/on-farm-data-cassava-nigeria.csv')

dataset <- dataset_cassava

# Apply filters to dataset

# We did not have authorization for Imo and discart Akwa Ibom

dataset <- dataset[!(dataset$state %in% c("Imo","Akwa Ibom")),] 

table(dataset$genotype)

dataset <- dataset[dataset$genotype!="Local",]

dataset <- dataset[!is.na(dataset$yieldton_ha),]

dataset <- dataset %>% group_by(farm) %>% mutate(valid_rank = sum(!is.na(plot))==3 & sum(!is.na(rank_freshrootsyield)) ==3 & length(unique(plot))==3 & length(unique(rank_freshrootsyield)) == 3 & plot %in% c('A','B','C') & rank_freshrootsyield %in% c(1,2,3))

dataset <- dataset %>% filter(valid_rank==TRUE)

dataset$year
dataset$state

length(unique(dataset$genotype))

nlc <- length(unique(dataset$genotype))*1/sum(1/table(dataset$genotype))

ranked_data <- rank_numeric(data = dataset,
                            items = 'genotype',
                            input = 'rank_freshrootsyield', 
                            id = 'farm',ascending = FALSE)

cassava_rankings <- ranked_data

dataset <-  dataset %>% group_by(farm) %>% mutate(ranking2 = round(rank(-yieldton_ha,ties.method = 'first'),0))

ranked_data2 <- rank_numeric(data = dataset,
                            items = 'genotype',
                            input = 'ranking2', 
                            id = 'farm',ascending = FALSE)

cassava_rankings2 <- ranked_data2

# metric data analysis

dataset$genotype <- factor(dataset$genotype)

dataset$farm <- factor(dataset$farm)

metric_mod_blups <- asreml(yieldton_ha~1,random = ~farm+genotype,data = dataset)

metric_mod_blues <- lmer(yieldton_ha~genotype+(1|farm),data = dataset)

blues <- fixef(metric_mod_blues)

pm <- predictmeans(metric_mod_blues,'genotype',plot = F)

pm$mean_table

# extra qsivar metric models

L <- diag(1,nrow = nrow(pm$mean_table))  # Contrast matrix

L[,1] <- 1

vcv <- vcov(metric_mod_blues)

VCOV <- L %*% as.matrix(vcv) %*% t(L)

bluesNams <- colnames(vcv)

bluesNams[1] <- levels(dataset$genotype)[1]

qSErrors <- qvcalc(VCOV,labels =  bluesNams,estimate = L%*%blues)

qSErrors <- as.data.frame(qSErrors$qvframe)

varcomp <- as.data.frame(summary(metric_mod_blues)$varcor)

sigma_e <- varcomp[varcomp$grp=='Residual',]$sdcor

tb_gen_rep <- table(dataset$genotype)

mean_table_cassava <- data.frame(qSErrors,sigma_e=sigma_e,mean_rep=mean(tb_gen_rep),max(tb_gen_rep),min(tb_gen_rep),sd(tb_gen_rep),nlc=nlc)

# Groundnut 

dataset_groundnut <- read.csv('Processed_data/On-farm-datasets/on-farm-data-groundnut-tanzania.csv')

dataset <- dataset_groundnut


dataset <- dataset[!is.na(dataset$yieldton_ha),]

dataset <- dataset %>% group_by(farm) %>% mutate(valid_rank = sum(!is.na(plot))==3 & sum(!is.na(rank_yield)) ==3 & length(unique(plot))==3 & length(unique(rank_yield)) == 3 & plot %in% c('A','B','C') & rank_yield %in% c(1,2,3))

dataset <- dataset %>% filter(valid_rank==TRUE)

unique(dataset$year)

length(unique(dataset$genotype))

nlc <- length(unique(dataset$genotype))*1/sum(1/table(dataset$genotype))

ranked_data <- rank_numeric(data = dataset,
                                 items = 'genotype',
                                 input = 'rank_yield', 
                                 id = 'farm',ascending = FALSE)

groundnut_rankings <- ranked_data

dataset <-  dataset %>% group_by(farm) %>% mutate(ranking2 = round(rank(-yieldton_ha,ties.method = 'first'),0))

ranked_data2 <- rank_numeric(data = dataset,
                             items = 'genotype',
                             input = 'ranking2', 
                             id = 'farm',ascending = FALSE)

groundnut_rankings2 <- ranked_data2

# metric data analysis

dataset$genotype <- factor(dataset$genotype)

dataset$farm <- factor(dataset$farm)

metric_mod_blups <- asreml(yieldton_ha~1,random = ~farm+genotype,data = dataset)

metric_mod_blues <- lmer(yieldton_ha~genotype+(1|farm),data = dataset)

blues <- fixef(metric_mod_blues)

pm <- predictmeans(metric_mod_blues,'genotype',plot = F)

pm$mean_table

# extra qsivar metric models

L <- diag(1,nrow = nrow(pm$mean_table))  # Contrast matrix

L[,1] <- 1

vcv <- vcov(metric_mod_blues)

VCOV <- L %*% as.matrix(vcv) %*% t(L)

bluesNams <- colnames(vcv)

bluesNams[1] <- levels(dataset$genotype)[1]

qSErrors <- qvcalc(VCOV,labels =  bluesNams,estimate = L%*%blues)

qSErrors <- as.data.frame(qSErrors$qvframe)

varcomp <- as.data.frame(summary(metric_mod_blues)$varcor)

sigma_e <- varcomp[varcomp$grp=='Residual',]$sdcor

tb_gen_rep <- table(dataset$genotype)

mean_table_groundnut <- data.frame(qSErrors,sigma_e=sigma_e,mean_rep=mean(tb_gen_rep),max(tb_gen_rep),min(tb_gen_rep),sd(tb_gen_rep),nlc=nlc)

#  maize

dataset_maize <- read.csv('Processed_data/On-farm-datasets/on-farm-data-maize-kenia.csv')

geno_info <- read_xlsx('Processed_data/Copy of genotype-information-maize-kenya_GMK_HD.xlsx')

dataset_maize <- dataset_maize %>% left_join(geno_info,by=c('genotype'='Genotype'))



# Early maize

dataset <- dataset_maize[dataset_maize$Maturity=='Early',]

dataset <- dataset[!is.na(dataset$yield_farmer),]

dataset <- dataset %>% group_by(farm) %>% mutate(valid_rank = sum(!is.na(plot))==3 & sum(!is.na(rank_yield)) ==3 & length(unique(plot))==3 & length(unique(rank_yield)) == 3 & plot %in% c('A','B','C') & rank_yield %in% c(1,2,3))

dataset <- dataset %>% filter(valid_rank==TRUE)

length(unique(dataset$genotype))

nlc <- length(unique(dataset$genotype))*1/sum(1/table(dataset$genotype))

ranked_data <- rank_numeric(data = dataset,
                            items = 'genotype',
                            input = 'rank_yield', 
                            id = 'farm',ascending = FALSE)


maize_early_rankings <- ranked_data

dataset <-  dataset %>% group_by(farm) %>% mutate(ranking2 = round(rank(-yield_farmer,ties.method = 'first'),0))

ranked_data2 <- rank_numeric(data = dataset,
                             items = 'genotype',
                             input = 'ranking2', 
                             id = 'farm',ascending = FALSE)

maize_early_rankings2 <- ranked_data2

# metric data analysis

dataset$genotype <- factor(dataset$genotype)

dataset$farm <- factor(dataset$farm)

metric_mod_blups <- asreml(yield_farmer~1,random = ~farm+genotype,data = dataset)

metric_mod_blues <- lmer(yield_farmer~genotype+(1|farm),data = dataset)

blues <- fixef(metric_mod_blues)

pm <- predictmeans(metric_mod_blues,'genotype',plot = F)

pm$mean_table

# extra qsivar metric models

L <- diag(1,nrow = nrow(pm$mean_table))  # Contrast matrix

L[,1] <- 1

vcv <- vcov(metric_mod_blues)

VCOV <- L %*% as.matrix(vcv) %*% t(L)

bluesNams <- colnames(vcv)

bluesNams[1] <- levels(dataset$genotype)[1]

qSErrors <- qvcalc(VCOV,labels =  bluesNams,estimate = L%*%blues)

qSErrors <- as.data.frame(qSErrors$qvframe)

varcomp <- as.data.frame(summary(metric_mod_blues)$varcor)

sigma_e <- varcomp[varcomp$grp=='Residual',]$sdcor

tb_gen_rep <- table(dataset$genotype)

mean_table_maize_early <- data.frame(qSErrors,sigma_e=sigma_e,mean_rep=mean(tb_gen_rep),max(tb_gen_rep),min(tb_gen_rep),sd(tb_gen_rep),nlc)



# Intermediate maize

dataset <- dataset_maize[dataset_maize$Maturity=='Intermediate',]

dataset <- dataset[!is.na(dataset$yield_farmer),]

dataset <- dataset %>% group_by(farm) %>% mutate(valid_rank = sum(!is.na(plot))==3 & sum(!is.na(rank_yield)) ==3 & length(unique(plot))==3 & length(unique(rank_yield)) == 3 & plot %in% c('A','B','C') & rank_yield %in% c(1,2,3))

dataset <- dataset %>% filter(valid_rank==TRUE)

length(unique(dataset$genotype))

ranked_data <- rank_numeric(data = dataset,
                            items = 'genotype',
                            input = 'rank_yield', 
                            id = 'farm',ascending = FALSE)

maize_intermedia_rankings <- ranked_data

dataset <-  dataset %>% group_by(farm) %>% mutate(ranking2 = round(rank(-yield_farmer,ties.method = 'first'),0))

ranked_data2 <- rank_numeric(data = dataset,
                             items = 'genotype',
                             input = 'ranking2', 
                             id = 'farm',ascending = FALSE)

maize_intermedia_rankings2 <- ranked_data2

nlc <- length(unique(dataset$genotype))*1/sum(1/table(dataset$genotype))

# metric data analysis

dataset$genotype <- factor(dataset$genotype)

dataset$farm <- factor(dataset$farm)

metric_mod_blups <- asreml(yield_farmer~1,random = ~farm+genotype,data = dataset)

metric_mod_blues <- lmer(yield_farmer~genotype+(1|farm),data = dataset)

blues <- fixef(metric_mod_blues)

pm <- predictmeans(metric_mod_blues,'genotype',plot = F)

pm$mean_table

# extra qsivar metric models

L <- diag(1,nrow = nrow(pm$mean_table))  # Contrast matrix

L[,1] <- 1

vcv <- vcov(metric_mod_blues)

VCOV <- L %*% as.matrix(vcv) %*% t(L)

bluesNams <- colnames(vcv)

bluesNams[1] <- levels(dataset$genotype)[1]

qSErrors <- qvcalc(VCOV,labels =  bluesNams,estimate = L%*%blues)

qSErrors <- as.data.frame(qSErrors$qvframe)

varcomp <- as.data.frame(summary(metric_mod_blues)$varcor)

sigma_e <- varcomp[varcomp$grp=='Residual',]$sdcor

tb_gen_rep <- table(dataset$genotype)

mean_table_maiz_intermedia <- data.frame(qSErrors,sigma_e=sigma_e,mean_rep=mean(tb_gen_rep),max(tb_gen_rep),min(tb_gen_rep),sd(tb_gen_rep),nlc=nlc)

# save rankings

rankings_datasets <- list(cassava_rankings=cassava_rankings,
        groundnut_rankings=groundnut_rankings,
        maize_early_rankings=maize_early_rankings,
        maize_intermedia_rankings=maize_intermedia_rankings
        )

saveRDS(rankings_datasets,file = 'Processed_data/rankings_datasets.rds')

rankings2_datasets <- list(cassava_rankings2=cassava_rankings2,
                          groundnut_rankings2=groundnut_rankings2,
                          maize_early_rankings2=maize_early_rankings2,
                          maize_intermedia_rankings2=maize_intermedia_rankings2
)

saveRDS(rankings2_datasets,file = 'Processed_data/rankings_datasets2.rds')


metric_tables <- list(mean_table_cassava=mean_table_cassava,mean_table_groundnut=mean_table_groundnut,
                      mean_table_maiz_intermedia=mean_table_maiz_intermedia,mean_table_maize_early=mean_table_maize_early)


saveRDS(metric_tables,file = 'Processed_data/metric_tables.rds')

# data analysis


source("master_funs.R")

rankings_datasets <- readRDS('Processed_data/rankings_datasets.rds')

rnks <- rankings_datasets$cassava_rankings

th_mod_out <- thurstonianMod(rnks = rnks,ref = sort(colnames(as.matrix(rnks)))[1])

cassava_rank_results <- th_mod_out

  
rnks <- rankings_datasets$groundnut_rankings

th_mod_out <- thurstonianMod(rnks = rnks,ref = sort(colnames(as.matrix(rnks)))[1])

groundnut_rank_results <- th_mod_out


rnks <- rankings_datasets$maize_early_rankings

th_mod_out <- thurstonianMod(rnks = rnks,ref = sort(colnames(as.matrix(rnks)))[1])

maize_early_rank_results <- th_mod_out

rnks <- rankings_datasets$maize_intermedia_rankings

th_mod_out <- thurstonianMod(rnks = rnks,ref = sort(colnames(as.matrix(rnks)))[1])

maize_interm_rank_results <- th_mod_out

#ls_results <- list(groundnut_rank_results=groundnut_rank_results,maize_early_rank_results=maize_early_rank_results,maize_interm_rank_results=maize_interm_rank_results)

#ls_results <- readRDS('Processed_data/results.rds')

#ls_results <- list(cassava_rank_results=cassava_rank_results,groundnut_rank_results=ls_results$groundnut_rank_results,maize_early_rank_results=ls_results$maize_early_rank_results,maize_interm_rank_results=ls_results$maize_interm_rank_results)

# 

#saveRDS(ls_results,'Processed_data/observed_data_paired.rds')

#------------------------------------------------------------------------------
#---- ranking 2

source("master_funs.R")

rankings_datasets <- readRDS('Processed_data/rankings_datasets2.rds')

rnks <- rankings_datasets$cassava_rankings2

th_mod_out <- thurstonianMod(rnks = rnks,ref = sort(colnames(as.matrix(rnks)))[1])

cassava_rank_results2 <- th_mod_out


rnks <- rankings_datasets$groundnut_rankings2

th_mod_out <- thurstonianMod(rnks = rnks,ref = sort(colnames(as.matrix(rnks)))[1])

groundnut_rank_results2 <- th_mod_out


rnks <- rankings_datasets$maize_early_rankings2

th_mod_out <- thurstonianMod(rnks = rnks,ref = sort(colnames(as.matrix(rnks)))[1])

maize_early_rank_results2 <- th_mod_out

rnks <- rankings_datasets$maize_intermedia_rankings2

th_mod_out <- thurstonianMod(rnks = rnks,ref = sort(colnames(as.matrix(rnks)))[1])

maize_interm_rank_results2 <- th_mod_out

ls_results2 <- list(cassava_rank_results2=cassava_rank_results2,groundnut_rank_results2=groundnut_rank_results2,maize_early_rank_results2=maize_early_rank_results2,maize_interm_rank_results2=maize_interm_rank_results2)

saveRDS('Processed_data/results.rds')

#ls_results <- readRDS('Processed_data/results.rds')

#ls_results <- list(cassava_rank_results=cassava_rank_results,groundnut_rank_results=ls_results$groundnut_rank_results,maize_early_rank_results=ls_results$maize_early_rank_results,maize_interm_rank_results=ls_results$maize_interm_rank_results)

# 

#saveRDS(ls_results,'Processed_data/observed_data_paired.rds')


#-------------------------------------------------------------------------------
# Here the real data analysis

ls_results <- readRDS('Processed_data/observed_data_paired.rds')

metric_table <- readRDS('Processed_data/metric_tables.rds')

# Cassava


rank_data <- ls_results$cassava_rank_result

metric_tab <- metric_table$mean_table_cassava

ut_means <- th_means_se(rank_data)

ut_means <- ut_means$qvframe

row.names(ut_means)[1] <- names(rank_data$est_mu)[1]

ut_means$genotype <- gsub('mu_','',row.names(ut_means))

metric_tab$genotype <- gsub('genotype','',row.names(metric_tab))

full_table <- inner_join(ut_means,metric_tab,by='genotype',suffix = c('rank','metric'))

# metric

weights_metric <- 1/full_table[,'quasiVarmetric']

full_table$genotype <- factor(full_table$genotype)

LM.met <- asreml(estimatemetric~1,random =~genotype , 
             family = (asr_gaussian(dispersion = 1)),weights =weights_metric,
             data=full_table)

LM.summ.met <-summary(LM.met)

nlc <- unique(full_table$nlc)

covMat.met <- as.data.frame(LM.summ.met$varcomp)

var_g_met <- covMat.met[row.names(covMat.met)=='genotype'  ,]$component

sd_e <- unique(full_table$sigma_e)

sd_g_met <-  sqrt(var_g_met)/(sd_e)

sd_g_met

qh2_met <- qH2(LM.met)

qh2_met

h2_met <- var_g_met/(var_g_met+sd_e^2/nlc)

h2_met

# ranking

full_table$weights_ranks <- 1/full_table[,'quasiVarrank']


LM.rank <- asreml(estimaterank~1,random = ~ genotype , 
                 family = (asr_gaussian(dispersion = 1)),weights =weights_ranks,
                 data=full_table,trace=F,na.action = na.method(x='omit'))

LM.summ.rank <-summary(LM.rank)

covMat.rank <- as.data.frame(LM.summ.rank$varcomp)

sd_g_rank <- sqrt(covMat.rank[row.names(covMat.rank)=='genotype'  ,]$component)

sd_g_rank

qh2_rank <- qH2(LM.rank)

qh2_rank 

h2_rank <- sd_g_rank^2/(sd_g_rank^2+1/nlc)

h2_rank

two_trait_mod <- rbind(
data.frame(estimates=full_table$estimatemetric,genotype=full_table$genotype,qVas = full_table$quasiVarmetric ,Type='Metric'),
data.frame(estimates=full_table$estimaterank,genotype=full_table$genotype,qVas = full_table$quasiVarrank, Type = 'Ranking' )
)

two_trait_mod$weight <-  1/two_trait_mod$qVas

two_trait_mod$genotype <- factor(two_trait_mod$genotype)

two_trait_mod$Type <- factor(two_trait_mod$Type)

model <- asreml(
  fixed = estimates~ Type,    random = ~us(Type):id(genotype) ,  weights =weight ,       # trait effect in fixed
  family =  (asr_gaussian(dispersion = 1)),                     # unstructured genetic covariance                # independent residuals
  data = two_trait_mod
)

varComp <- summary(model)$varcomp

gen_cor <- varComp['Type:genotype!Type_Ranking:Metric',]$comp/(sqrt(varComp['Type:genotype!Type_Metric:Metric',]$comp)*sqrt(varComp['Type:genotype!Type_Ranking:Ranking',]$comp))

cassava_data <- c(sd_g_met=sd_g_met,sd_g_rank=sd_g_rank,qh2_rank=qh2_rank,
                  qh2_met=qh2_met,h2_met=h2_met,h_rank=h2_rank,
                  medrep=unique(full_table$mean_rep),
                  minrep=unique(full_table$min.tb_gen_rep.),
                  maxrep=unique(full_table$max.tb_gen_rep.),
                  sdrep=unique(full_table$sd.tb_gen_rep.),gen_cor=gen_cor
                  )

# Groundnut 

rank_data <- ls_results$groundnut_rank_results

metric_tab <- metric_table$mean_table_groundnut

ut_means <- th_means_se(rank_data)

ut_means <- ut_means$qvframe

row.names(ut_means)[1] <- names(rank_data$est_mu)[1]

ut_means$genotype <- gsub('mu_','',row.names(ut_means))

metric_tab$genotype <- gsub('genotype','',row.names(metric_tab))

full_table <- inner_join(ut_means,metric_tab,by='genotype',suffix = c('rank','metric'))

# metric

weights_metric <- 1/full_table[,'quasiVarmetric']

full_table$genotype <- factor(full_table$genotype)

LM.met <- asreml(estimatemetric~1,random =~genotype , 
                 family = (asr_gaussian(dispersion = 1)),weights =weights_metric,
                 data=full_table)

LM.summ.met <-summary(LM.met)

covMat.met <- as.data.frame(LM.summ.met$varcomp)

var_g <- covMat.met[row.names(covMat.met)=='genotype'  ,]$component

sd_g_met <- sqrt(var_g) / unique(full_table$sigma_e)

sd_g_met

h2_met <- var_g/(var_g+unique(full_table$sigma_e)^2/unique(full_table$nlc))


qh2_met <- qH2(LM.met)

qh2_met




# ranking

full_table$weights_ranks <- 1/full_table[,'quasiVarrank']


LM.rank <- asreml(estimaterank~1,random = ~ genotype , 
                  family = (asr_gaussian(dispersion = 1)),weights =weights_ranks,
                  data=full_table,trace=F,na.action = na.method(x='omit'))

LM.summ.rank <-summary(LM.rank)

covMat.rank <- as.data.frame(LM.summ.rank$varcomp)

sd_g_rank <- sqrt(covMat.rank[row.names(covMat.rank)=='genotype'  ,]$component)

sd_g_rank

qh2_rank <- qH2(LM.rank)

qh2_rank 


h2_rank <- sd_g_rank^2/(sd_g_rank^2+1/unique(full_table$nlc))

h2_rank

two_trait_mod <- rbind(
  data.frame(estimates=full_table$estimatemetric,genotype=full_table$genotype,qVas = full_table$quasiVarmetric ,Type='Metric'),
  data.frame(estimates=full_table$estimaterank,genotype=full_table$genotype,qVas = full_table$quasiVarrank, Type = 'Ranking' )
)

two_trait_mod$weight <-  1/two_trait_mod$qVas

two_trait_mod$genotype <- factor(two_trait_mod$genotype)

two_trait_mod$Type <- factor(two_trait_mod$Type)

model <- asreml(
  fixed = estimates~ Type,    random = ~us(Type):id(genotype) ,  weights =weight ,       # trait effect in fixed
  family =  (asr_gaussian(dispersion = 1)),                     # unstructured genetic covariance                # independent residuals
  data = two_trait_mod
)

varComp <- summary(model)$varcomp

gen_cor <- varComp['Type:genotype!Type_Ranking:Metric',]$comp/(sqrt(varComp['Type:genotype!Type_Metric:Metric',]$comp)*sqrt(varComp['Type:genotype!Type_Ranking:Ranking',]$comp))



groudnut_data <- c(sd_g_met=sd_g_met,sd_g_rank=sd_g_rank,qh2_rank=qh2_rank,
                   qh2_met=qh2_met,h2_met=h2_met,h2_rank,medrep=unique(full_table$mean_rep),
                   minrep=unique(full_table$min.tb_gen_rep.),
                   maxrep=unique(full_table$max.tb_gen_rep.),
                   sdrep=unique(full_table$sd.tb_gen_rep.),gen_cor=gen_cor)


# Early maize

rank_data <- ls_results$maize_early_rank_results

metric_tab <- metric_table$mean_table_maize_early

ut_means <- th_means_se(rank_data)

ut_means <- ut_means$qvframe

row.names(ut_means)[1] <- names(rank_data$est_mu)[1]

ut_means$genotype <- gsub('mu_','',row.names(ut_means))

metric_tab$genotype <- gsub('genotype','',row.names(metric_tab))

full_table <- inner_join(ut_means,metric_tab,by='genotype',suffix = c('rank','metric'))

# metric

weights_metric <- 1/full_table[,'quasiVarmetric']

full_table$genotype <- factor(full_table$genotype)

LM.met <- asreml(estimatemetric~1,random =~genotype , 
                 family = (asr_gaussian(dispersion = 1)),weights =weights_metric,
                 data=full_table)

LM.summ.met <-summary(LM.met)

covMat.met <- as.data.frame(LM.summ.met$varcomp)

var_g <- covMat.met[row.names(covMat.met)=='genotype'  ,]$component

sd_g_met <- sqrt(var_g) / unique(full_table$sigma_e)

sd_g_met

h2_met <- var_g/(var_g+unique(full_table$sigma_e)^2/unique(full_table$nlc))

h2_met

qh2_met <- qH2(LM.met)

qh2_met

# ranking

full_table$weights_ranks <- 1/full_table[,'quasiVarrank']


LM.rank <- asreml(estimaterank~1,random = ~ genotype , 
                  family = (asr_gaussian(dispersion = 1)),weights =weights_ranks,
                  data=full_table,trace=F,na.action = na.method(x='omit'))

LM.summ.rank <-summary(LM.rank)

covMat.rank <- as.data.frame(LM.summ.rank$varcomp)

sd_g_rank <- sqrt(covMat.rank[row.names(covMat.rank)=='genotype'  ,]$component)

sd_g_rank

h2_rank <- sd_g_rank^2/(sd_g_rank^2+1/unique(full_table$nlc))

h2_rank

qh2_rank <- qH2(LM.rank)

qh2_rank 

two_trait_mod <- rbind(
  data.frame(estimates=full_table$estimatemetric,genotype=full_table$genotype,qVas = full_table$quasiVarmetric ,Type='Metric'),
  data.frame(estimates=full_table$estimaterank,genotype=full_table$genotype,qVas = full_table$quasiVarrank, Type = 'Ranking' )
)

two_trait_mod$weight <-  1/two_trait_mod$qVas

two_trait_mod$genotype <- factor(two_trait_mod$genotype)

two_trait_mod$Type <- factor(two_trait_mod$Type)

model <- asreml(
  fixed = estimates~ Type,    random = ~us(Type):id(genotype) ,  weights =weight ,       # trait effect in fixed
  family =  (asr_gaussian(dispersion = 1)),                     # unstructured genetic covariance                # independent residuals
  data = two_trait_mod
)

varComp <- summary(model)$varcomp

gen_cor <- varComp['Type:genotype!Type_Ranking:Metric',]$comp/(sqrt(varComp['Type:genotype!Type_Metric:Metric',]$comp)*sqrt(varComp['Type:genotype!Type_Ranking:Ranking',]$comp))


early_maize_data <- c(sd_g_met=sd_g_met,sd_g_rank=sd_g_rank,qh2_rank=qh2_rank,
                      qh2_met=qh2_met,h2_met=h2_met,h2_rank=h2_rank,medrep=unique(full_table$mean_rep),
                      minrep=unique(full_table$min.tb_gen_rep.),
                      maxrep=unique(full_table$max.tb_gen_rep.),sdrep=unique(full_table$sd.tb_gen_rep.),gen_cor=gen_cor)


# intermerdiate maize maize

rank_data <- ls_results$maize_interm_rank_results

metric_tab <- metric_table$mean_table_maiz_intermedia

ut_means <- th_means_se(rank_data)

ut_means <- ut_means$qvframe

row.names(ut_means)[1] <- names(rank_data$est_mu)[1]

ut_means$genotype <- gsub('mu_','',row.names(ut_means))

metric_tab$genotype <- gsub('genotype','',row.names(metric_tab))

full_table <- inner_join(ut_means,metric_tab,by='genotype',suffix = c('rank','metric'))




# metric

weights_metric <- 1/full_table[,'quasiVarmetric']

full_table$genotype <- factor(full_table$genotype)

LM.met <- asreml(estimatemetric~1,random =~genotype , 
                 family = (asr_gaussian(dispersion = 1)),weights =weights_metric,
                 data=full_table)


LM.met <- update(LM.met)

LM.summ.met <-summary(LM.met)

covMat.met <- as.data.frame(LM.summ.met$varcomp)

var_g <- covMat.met[row.names(covMat.met)=='genotype'  ,]$component

sd_g_met <- sqrt(var_g) / unique(full_table$sigma_e)

sd_g_met

h2_met <- var_g/(var_g+unique(full_table$sigma_e)^2/unique(full_table$nlc))

h2_met

qh2_met <- qH2(LM.met)

qh2_met

# ranking

full_table$weights_ranks <- 1/full_table[,'quasiVarrank']


LM.rank <- asreml(estimaterank~1,random = ~ genotype , 
                  family = (asr_gaussian(dispersion = 1)),weights =weights_ranks,
                  data=full_table,trace=F,na.action = na.method(x='omit'))

LM.summ.rank <-summary(LM.rank)

covMat.rank <- as.data.frame(LM.summ.rank$varcomp)

sd_g_rank <- sqrt(covMat.rank[row.names(covMat.rank)=='genotype'  ,]$component)

sd_g_rank

h2_rank <- sd_g_rank^2/(sd_g_rank^2+1/unique(full_table$nlc))

h2_rank

qh2_rank <- qH2(LM.rank)

qh2_rank 

two_trait_mod <- rbind(
  data.frame(estimates=full_table$estimatemetric,genotype=full_table$genotype,qVas = full_table$quasiVarmetric ,Type='Metric'),
  data.frame(estimates=full_table$estimaterank,genotype=full_table$genotype,qVas = full_table$quasiVarrank, Type = 'Ranking' )
)

two_trait_mod$weight <-  1/two_trait_mod$qVas

two_trait_mod$genotype <- factor(two_trait_mod$genotype)

two_trait_mod$Type <- factor(two_trait_mod$Type)

model <- asreml(
  fixed = estimates~ Type,    random = ~us(Type):id(genotype) ,  weights =weight ,       # trait effect in fixed
  family =  (asr_gaussian(dispersion = 1)),                     # unstructured genetic covariance                # independent residuals
  data = two_trait_mod
)

varComp <- summary(model)$varcomp

gen_cor <- varComp['Type:genotype!Type_Ranking:Metric',]$comp/(sqrt(varComp['Type:genotype!Type_Metric:Metric',]$comp)*sqrt(varComp['Type:genotype!Type_Ranking:Ranking',]$comp))


intermed_maize_data <- c(sd_g_met=sd_g_met,sd_g_rank=sd_g_rank,
                         qh2_rank=qh2_rank,qh2_met=qh2_met,h2_met=h2_met,h2_rank=h2_rank,
                         medrep=unique(full_table$mean_rep),minrep=unique(full_table$min.tb_gen_rep.),maxrep=unique(full_table$max.tb_gen_rep.),sdrep=unique(full_table$sd.tb_gen_rep.),gen_cor=gen_cor)

rbind(
cassava_data,
groudnut_data,
early_maize_data,
intermed_maize_data
)


######################################################################
#########################################################
#################ranking 2

ls_results2 <- readRDS('Processed_data/results2.rds')

# Here the real data analysis

#ls_results <- readRDS('Processed_data/observed_data_paired.rds')

metric_table <- readRDS('Processed_data/metric_tables.rds')

# Cassava


rank_data <- ls_results2$cassava_rank_results2

# ranking


ut_means <- th_means_se(rank_data)

ut_means <- ut_means$qvframe

row.names(ut_means)[1] <- names(rank_data$est_mu)[1]

ut_means$genotype <- gsub('mu_','',row.names(ut_means))


full_table <- ut_means

full_table$weights_ranks <- 1/full_table$quasiVar 

full_table$genotype <- factor(full_table$genotype)

LM.rank <- asreml(estimate~1,random = ~ genotype , 
                  family = (asr_gaussian(dispersion = 1)),weights =weights_ranks,
                  data=full_table,trace=F,na.action = na.method(x='omit'))

LM.summ.rank <-summary(LM.rank)

covMat.rank <- as.data.frame(LM.summ.rank$varcomp)

sd_g_rank <- sqrt(covMat.rank[row.names(covMat.rank)=='genotype'  ,]$component)

sd_g_rank

h2_rank <- sd_g_rank^2/(sd_g_rank^2+1/unique(metric_table$mean_table_cassava$nlc))

h2_rank

qh2_rank <- qH2(LM.rank)

qh2_rank 


# Groundnut 

rank_data <- ls_results2$groundnut_rank_results2

# ranking


ut_means <- th_means_se(rank_data)

ut_means <- ut_means$qvframe

row.names(ut_means)[1] <- names(rank_data$est_mu)[1]

ut_means$genotype <- gsub('mu_','',row.names(ut_means))


full_table <- ut_means

full_table$weights_ranks <- 1/full_table$quasiVar 

full_table$genotype <- factor(full_table$genotype)

LM.rank <- asreml(estimate~1,random = ~ genotype , 
                  family = (asr_gaussian(dispersion = 1)),weights =weights_ranks,
                  data=full_table,trace=F,na.action = na.method(x='omit'))

LM.summ.rank <-summary(LM.rank)

covMat.rank <- as.data.frame(LM.summ.rank$varcomp)

sd_g_rank <- sqrt(covMat.rank[row.names(covMat.rank)=='genotype'  ,]$component)

sd_g_rank

h2_rank <- sd_g_rank^2/(sd_g_rank^2+1/unique(metric_table$mean_table_groundnut$nlc))

h2_rank

qh2_rank <- qH2(LM.rank)

qh2_rank 
# Early maize

rank_data <- ls_results2$maize_early_rank_results2

# ranking


ut_means <- th_means_se(rank_data)

ut_means <- ut_means$qvframe

row.names(ut_means)[1] <- names(rank_data$est_mu)[1]

ut_means$genotype <- gsub('mu_','',row.names(ut_means))


full_table <- ut_means

full_table$weights_ranks <- 1/full_table$quasiVar 

full_table$genotype <- factor(full_table$genotype)

LM.rank <- asreml(estimate~1,random = ~ genotype , 
                  family = (asr_gaussian(dispersion = 1)),weights =weights_ranks,
                  data=full_table,trace=F,na.action = na.method(x='omit'))

LM.summ.rank <-summary(LM.rank)

covMat.rank <- as.data.frame(LM.summ.rank$varcomp)

sd_g_rank <- sqrt(covMat.rank[row.names(covMat.rank)=='genotype'  ,]$component)

sd_g_rank

h2_rank <-sd_g_rank^2/(sd_g_rank^2+1/unique(metric_table$mean_table_maize_early$nlc))

h2_rank

qh2_rank <- qH2(LM.rank)

qh2_rank 

# intermerdiate maize maize

rank_data <- ls_results2$maize_interm_rank_results2

# ranking


ut_means <- th_means_se(rank_data)

ut_means <- ut_means$qvframe

row.names(ut_means)[1] <- names(rank_data$est_mu)[1]

ut_means$genotype <- gsub('mu_','',row.names(ut_means))


full_table <- ut_means

full_table$weights_ranks <- 1/full_table$quasiVar 

full_table$genotype <- factor(full_table$genotype)

LM.rank <- asreml(estimate~1,random = ~ genotype , 
                  family = (asr_gaussian(dispersion = 1)),weights =weights_ranks,
                  data=full_table,trace=F,na.action = na.method(x='omit'))

LM.summ.rank <-summary(LM.rank)

covMat.rank <- as.data.frame(LM.summ.rank$varcomp)


h2_rank

sd_g_rank <- sqrt(covMat.rank[row.names(covMat.rank)=='genotype'  ,]$component)

sd_g_rank

h2_rank <- sd_g_rank^2/(sd_g_rank^2+1/unique(metric_table$mean_table_maiz_intermedia$nlc))

h2_rank

qh2_rank <- qH2(LM.rank)

qh2_rank 

##------------------------------sweeet potato data analysis--------------------

dataset_sweet_potato <- read.csv('Processed_data/On-farm-datasets/sweet_potato_processed.csv')



table(dataset_sweet_potato$genotype,dataset_sweet_potato$rank)

presence <- dataset_sweet_potato %>% select(id, genotype)%>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = genotype, values_from = presence, values_fill = 0)

# Remove the farm column
mat <- as.matrix(presence[,-1])


# split by project

pro1 <- dataset_sweet_potato[dataset_sweet_potato$projectid=='2ea3466e3686',]



table(pro1$genotype,pro1$rank)

presence <- pro1 %>% select(id, genotype)%>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = genotype, values_from = presence, values_fill = 0)

# Remove the farm column
mat <- as.matrix(presence[,-1])

pro2 <- dataset_sweet_potato[dataset_sweet_potato$projectid=='c8d22e12ddd0',]

table(pro2$genotype,pro2$rank)

presence <- pro2 %>% select(id, genotype)%>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = genotype, values_from = presence, values_fill = 0)

# Remove the farm column
mat <- as.matrix(presence[,-1])

# they are not connected, I had to fixed



# Cross-product gives co-occurrence counts
co_occurrence <- t(mat) %*% mat
co_occurrence

# select only pro1

dataset_sweet_potato <- pro2


dataset <- dataset_sweet_potato

nlc <- length(unique(dataset$genotype))*1/sum(1/table(dataset$genotype))

1.96*sd(table(dataset$genotype))/sqrt(8)

ranked_data <- rank_numeric(data = dataset,
                            items = 'genotype',
                            input = 'rank', 
                            id = 'id',ascending = FALSE)


ranked_data2 <- rank_numeric(data = dataset,
                             items = 'genotype',
                             input = 'rank2', 
                             id = 'id',ascending = FALSE)

# metric data analysis

dataset$genotype <- factor(dataset$genotype)

dataset$farm <- factor(dataset$id)

metric_mod_blups <- asreml(yield ~1,random = ~farm+genotype,data = dataset)

metric_mod_blues <- lmer(yield~genotype+(1|farm),data = dataset)

blues <- fixef(metric_mod_blues)

pm <- predictmeans(metric_mod_blues,'genotype',plot = F)

pm$mean_table

# extra qsivar metric models

L <- diag(1,nrow = nrow(pm$mean_table))  # Contrast matrix

L[,1] <- 1

vcv <- vcov(metric_mod_blues)

VCOV <- L %*% as.matrix(vcv) %*% t(L)

bluesNams <- colnames(vcv)

bluesNams[1] <- levels(dataset$genotype)[1]

qSErrors <- qvcalc(VCOV,labels =  bluesNams,estimate = L%*%blues)

qSErrors <- as.data.frame(qSErrors$qvframe)

varcomp <- as.data.frame(summary(metric_mod_blues)$varcor)

sigma_e <- varcomp[varcomp$grp=='Residual',]$sdcor

tb_gen_rep <- table(dataset$genotype)

mean_table_sweet_potato <- data.frame(qSErrors,sigma_e=sigma_e,mean_rep=mean(tb_gen_rep),max(tb_gen_rep),min(tb_gen_rep),sd(tb_gen_rep),nlc=nlc)


# --------------------------------data analysis-----------------------------

#rankings_datasets <- readRDS('Processed_data/rankings_datasets.rds')

rnks <- ranked_data

th_mod_out <- thurstonianMod(rnks = rnks,ref = 'Local_check')

sweet_potato_rank_results <- th_mod_out

# rankings 2

rnks <- ranked_data2

th_mod_out <- thurstonianMod(rnks = rnks,ref = 'Local_check')

sweet_potato_rank_results2 <- th_mod_out

# intermerdiate maize maize

rank_data <- sweet_potato_rank_results

rank_data2 <- sweet_potato_rank_results2

metric_tab <- mean_table_sweet_potato

# ranking 1

ut_means <- th_means_se(rank_data)

ut_means <- ut_means$qvframe

row.names(ut_means)[1] <- names(rank_data$est_mu)[1]

ut_means$genotype <- gsub('mu_','',row.names(ut_means))

metric_tab$genotype <- gsub('genotype','',row.names(metric_tab))

full_table <- inner_join(ut_means,metric_tab,by='genotype',suffix = c('rank','metric'))

# ranking 2

ut_means2 <- th_means_se(rank_data2)

ut_means2 <- ut_means2$qvframe

row.names(ut_means2)[1] <- names(rank_data2$est_mu)[1]

ut_means2$genotype <- gsub('mu_','',row.names(ut_means2))


# metric

weights_metric <- 1/full_table[,'quasiVarmetric']

full_table$genotype <- factor(full_table$genotype)

LM.met <- asreml(estimatemetric~1,random =~genotype , 
                 family = (asr_gaussian(dispersion = 1)),weights =weights_metric,
                 data=full_table)



#LM.met <- update(LM.met)

LM.summ.met <-summary(LM.met)

covMat.met <- as.data.frame(LM.summ.met$varcomp)

var_g <- covMat.met[row.names(covMat.met)=='genotype'  ,]$component

sd_g_met <- sqrt(var_g) / unique(full_table$sigma_e)

sd_g_met

h2_met <- var_g/(var_g+unique(full_table$sigma_e)^2/unique(full_table$nlc))

h2_met

qh2_met <- qH2(LM.met)

qh2_met

# ranking

full_table$weights_ranks <- 1/full_table[,'quasiVarrank']


LM.rank <- asreml(estimaterank~1,random = ~ genotype , 
                  family = (asr_gaussian(dispersion = 1)),weights =weights_ranks,
                  data=full_table,trace=F,na.action = na.method(x='omit'))

LM.summ.rank <-summary(LM.rank)

covMat.rank <- as.data.frame(LM.summ.rank$varcomp)

sd_g_rank <- sqrt(covMat.rank[row.names(covMat.rank)=='genotype'  ,]$component)

sd_g_rank

h2_rank <- sd_g_rank^2/(sd_g_rank^2+1/unique(full_table$nlc))

h2_rank

qh2_rank <- qH2(LM.rank)

qh2_rank 

# ranking 2

ut_means2$weights_ranks <- 1/ut_means2[,'quasiVar']

ut_means2$genotype <- factor(ut_means2$genotype)

LM.rank <- asreml(estimate~1,random = ~ genotype , 
                  family = (asr_gaussian(dispersion = 1)),weights =weights_ranks,
                  data=ut_means2,trace=F,na.action = na.method(x='omit'))

LM.summ.rank <-summary(LM.rank)

covMat.rank <- as.data.frame(LM.summ.rank$varcomp)

sd_g_rank <- sqrt(covMat.rank[row.names(covMat.rank)=='genotype'  ,]$component)

sd_g_rank

h2_rank2 <- sd_g_rank^2/(sd_g_rank^2+1/unique(full_table$nlc))

h2_rank2

qh2_rank2 <- qH2(LM.rank)

qh2_rank2 

# genetic correlation

two_trait_mod <- rbind(
  data.frame(estimates=full_table$estimatemetric,genotype=full_table$genotype,qVas = full_table$quasiVarmetric ,Type='Metric'),
  data.frame(estimates=full_table$estimaterank,genotype=full_table$genotype,qVas = full_table$quasiVarrank, Type = 'Ranking' )
)

two_trait_mod$weight <-  1/two_trait_mod$qVas

two_trait_mod$genotype <- factor(two_trait_mod$genotype)

two_trait_mod$Type <- factor(two_trait_mod$Type)

model <- asreml(
  fixed = estimates~ Type,    random = ~us(Type):id(genotype) ,  weights =weight ,       # trait effect in fixed
  family =  (asr_gaussian(dispersion = 1)),                     # unstructured genetic covariance                # independent residuals
  data = two_trait_mod
)

varComp <- summary(model)$varcomp

gen_cor <- varComp['Type:genotype!Type_Ranking:Metric',]$comp/(sqrt(varComp['Type:genotype!Type_Metric:Metric',]$comp)*sqrt(varComp['Type:genotype!Type_Ranking:Ranking',]$comp))


sweet_potatodata <- c(sd_g_met=sd_g_met,sd_g_rank=sd_g_rank,
                         qh2_rank=qh2_rank,qh2_met=qh2_met,h2_met=h2_met,h2_rank=h2_rank,
                         medrep=unique(full_table$mean_rep),minrep=unique(full_table$min.tb_gen_rep.),maxrep=unique(full_table$max.tb_gen_rep.),sdrep=unique(full_table$sd.tb_gen_rep.),gen_cor=gen_cor)


sweet_potatodata


rbind(
cassava_data,
groudnut_data,
early_maize_data,
intermed_maize_data,
sweet_potatodata
)






