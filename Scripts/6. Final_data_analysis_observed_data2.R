
# Data analysis
# all crops
# Hugo Dorado
# 11/15/2024

library(asreml)
library(tidyverse)
library(lme4)

setwd('Paper_1_Estimation_of_genetic_parameters')

source("Scripts/0. Ranking_pars_FUN.R")

qH2 <- function(model){
  varcomp <- summary(model)$varcomp
  
  BLUPl <- summary(model, coef=TRUE)$coef.random
  BLUPl <-as.data.frame(BLUPl)
  BLUPl <- BLUPl[grep('genotype_',rownames(BLUPl)),]
  
  BLUPl$PEV <- BLUPl$std.error^2
  
  PEV1<-BLUPl$PEV
  
  1-mean(PEV1)/varcomp['genotype','component']
}

# ----------------------------------- Cassava ----------------------------------

blues_metric_data <- read.csv("Results/Observed_dataset_results/blues_cassava.csv")

yield_ranked <- read.csv('Results/Observed_dataset_results/cassava_worth_values.csv')

fulldataset <- full_join(blues_metric_data,yield_ranked,by="genotype",
                         suffix = c('_blues',"_yieldrank"))

cor(fulldataset$estimate_blues,fulldataset$estimate_yieldrank,use = "pairwise.complete.obs")

fulldataset %>% ggplot(aes(x=estimate_blues,y=estimate_yieldrank)) + 
  theme_bw() + geom_text(aes(label = genotype)) + geom_smooth(method = 'lm',se=F) +
  xlab('Genotypic means ton/ha') +ylab('Genotypic worth values')

## Metric data

# genotypic variance

fulldataset$yieldrank_weights <- 1/fulldataset$SE_yieldrank^2

fulldataset$yieldmetric_weights <- 1/fulldataset$SE_blues^2

fulldataset$yieldrank_weights_q <- 1/fulldataset$quasiSE_yieldrank^2

fulldataset$genotype <- factor(fulldataset$genotype)
  
blups <- asreml(fixed = estimate_blues ~ 1, random = ~ genotype , weights = yieldmetric_weights,
       family = (asr_gaussian(dispersion = 1)),
       data = fulldataset,trace=F,na.action = na.method(x='omit'))

summary(blups)$varcomp

# 


# heritabilitey

qH2(blups)

## Ranking data se

fulldataset$genotype <- factor(fulldataset$genotype)

blups <- asreml(fixed = estimate_yieldrank ~ 1, random = ~ genotype , weights = yieldrank_weights,
                family = (asr_gaussian(dispersion = 1)),
                data = fulldataset[-1,],trace=F,na.action = na.method(x='omit'))

summary(blups)$varcomp

# heritabilitey

qH2(blups)

## Ranking data qse

fulldataset$genotype <- factor(fulldataset$genotype)

blups <- asreml(fixed = estimate_yieldrank ~ 1, random = ~ genotype , weights = yieldrank_weights_q,
                family = (asr_gaussian(dispersion = 1)),
                data = fulldataset,trace=F,na.action = na.method(x='omit'))

summary(blups)$varcomp

# heritabilitey

qH2(blups)

# ------------------------------Groudnut------------------------------

blues_metric_data <- read.csv("Results/Observed_dataset_results/blues_groundnut.csv")

yield_ranked <- read.csv('Results/Observed_dataset_results/groundnut_worth_values.csv')

fulldataset <- full_join(blues_metric_data,yield_ranked,by="genotype",
                         suffix = c('_blues',"_yieldrank"))

cor(fulldataset$estimate_blues,fulldataset$estimate_yieldrank,use = "pairwise.complete.obs")

fulldataset %>% ggplot(aes(x=estimate_blues,y=estimate_yieldrank)) + 
  theme_bw() + geom_text(aes(label = genotype)) + geom_smooth(method = 'lm',se=F) +
  xlab('Genotypic means ton/ha') +ylab('Genotypic worth values')

## Metric data

# genotypic variance

fulldataset$yieldrank_weights <- 1/fulldataset$SE_yieldrank^2

fulldataset$yieldmetric_weights <- 1/fulldataset$SE_blues^2

fulldataset$yieldrank_weights_q <- 1/fulldataset$quasiSE_yieldrank^2

fulldataset$genotype <- factor(fulldataset$genotype)

blups <- asreml(fixed = estimate_blues ~ 1, random = ~ genotype , weights = yieldmetric_weights,
                family = (asr_gaussian(dispersion = 1)),
                data = fulldataset,trace=F,na.action = na.method(x='omit'))

summary(blups)$varcomp

# heritabilitey

qH2(blups)

## Ranking data

fulldataset$genotype <- factor(fulldataset$genotype)

blups <- asreml(fixed = estimate_yieldrank ~ 1, random = ~ genotype , weights = yieldrank_weights,
                family = (asr_gaussian(dispersion = 1)),
                data = fulldataset[-1,],trace=F,na.action = na.method(x='omit'))

summary(blups)$varcomp

# heritabilitey

qH2(blups)


## Ranking data

fulldataset$genotype <- factor(fulldataset$genotype)

blups <- asreml(fixed = estimate_yieldrank ~ 1, random = ~ genotype , weights = yieldrank_weights_q,
                family = (asr_gaussian(dispersion = 1)),
                data = fulldataset,trace=F,na.action = na.method(x='omit'))

summary(blups)$varcomp

# heritabilitey

qH2(blups)



#---------------------------------Maize-------------------------------

# Early

yield_ranked <- read.csv("Results/Observed_dataset_results/maize_early_worth_values.csv")

blues_metric_data <- read.csv('Results/Observed_dataset_results/blues_maize_Early.csv')

fulldataset <- full_join(blues_metric_data,yield_ranked,by="genotype",
                         suffix = c('_blues',"_yieldrank"))

cor(fulldataset$estimate_blues,fulldataset$estimate_yieldrank,use = "pairwise.complete.obs")

fulldataset %>% ggplot(aes(x=estimate_blues,y=estimate_yieldrank)) + 
  theme_bw() + geom_text(aes(label = genotype)) + geom_smooth(method = 'lm',se=F) +
  xlab('Genotypic means ton/ha') +ylab('Genotypic worth values')



## Metric data

# genotypic variance

fulldataset$yieldrank_weights <- 1/fulldataset$SE_yieldrank^2

fulldataset$yieldmetric_weights <- 1/fulldataset$SE_blues^2

fulldataset$yieldrank_weights_q <- 1/fulldataset$quasiSE_yieldrank^2

fulldataset$genotype <- factor(fulldataset$genotype)

blups <- asreml(fixed = estimate_blues ~ 1, random = ~ genotype , weights = yieldmetric_weights,
                family = (asr_gaussian(dispersion = 1)),
                data = fulldataset,trace=F,na.action = na.method(x='omit'))

summary(blups)$varcomp

# heritabilitey

qH2(blups)

## Ranking data

fulldataset$genotype <- factor(fulldataset$genotype)

blups <- asreml(fixed = estimate_yieldrank ~ 1, random = ~ genotype , weights = yieldrank_weights,
                family = (asr_gaussian(dispersion = 1)),
                data = fulldataset[-13,],trace=F,na.action = na.method(x='omit'))

summary(blups)$varcomp

# heritabilitey

qH2(blups)


## Ranking data

fulldataset$genotype <- factor(fulldataset$genotype)

blups <- asreml(fixed = estimate_yieldrank ~ 1, random = ~ genotype , weights = yieldrank_weights_q,
                family = (asr_gaussian(dispersion = 1)),
                data = fulldataset,trace=F,na.action = na.method(x='omit'))

summary(blups)$varcomp

# heritabilitey

qH2(blups)


# Intermediate---

blues_metric_data <- read.csv("Results/Observed_dataset_results/blues_maize_Intermediate.csv")

yield_ranked <- read.csv('Results/Observed_dataset_results/maize_Intermediate_worth_values.csv')

fulldataset <- full_join(blues_metric_data,yield_ranked,by="genotype",
                         suffix = c('_blues',"_yieldrank"))

cor(fulldataset$estimate_blues,fulldataset$estimate_yieldrank,use = "pairwise.complete.obs")

fulldataset %>% ggplot(aes(x=estimate_blues,y=estimate_yieldrank)) + 
  theme_bw() + geom_text(aes(label = genotype)) + geom_smooth(method = 'lm',se=F) +
  xlab('Genotypic means ton/ha') +ylab('Genotypic worth values')
  
## Metric data

# genotypic variance

fulldataset$yieldrank_weights <- 1/fulldataset$SE_yieldrank^2

fulldataset$yieldmetric_weights <- 1/fulldataset$SE_blues^2

fulldataset$yieldrank_weights_q <- 1/fulldataset$quasiSE_yieldrank^2

fulldataset$genotype <- factor(fulldataset$genotype)

blups <- asreml(fixed = estimate_blues ~ 1, random = ~ genotype , weights = yieldmetric_weights,
                family = (asr_gaussian(dispersion = 1)),
                data = fulldataset,trace=F,na.action = na.method(x='omit'))

summary(blups)$varcomp


saas <- asreml(fulldataset)


# heritabilitey

qH2(blups)

## Ranking data

fulldataset$genotype <- factor(fulldataset$genotype)

blups <- asreml(fixed = estimate_yieldrank ~ 1, random = ~ genotype , weights = yieldrank_weights,
                family = (asr_gaussian(dispersion = 1)),
                data = fulldataset[-1,],trace=F,na.action = na.method(x='omit'))

summary(blups)$varcomp

# heritabilitey

qH2(blups)


## Ranking data

fulldataset$genotype <- factor(fulldataset$genotype)

blups <- asreml(fixed = estimate_yieldrank ~ 1, random = ~ genotype , weights = yieldrank_weights_q,
                family = (asr_gaussian(dispersion = 1)),
                data = fulldataset,trace=F,na.action = na.method(x='omit'))

summary(blups)$varcomp

# heritabilitey

qH2(blups)




# Intermediate ---------- discarted -------------------------

blues_metric_data <- read.csv("Results/Observed_dataset_results/blues_maize_Late.csv")

yield_ranked <- read.csv('Results/Observed_dataset_results/maize_Late_worth_values.csv')

fulldataset <- full_join(blues_metric_data,yield_ranked,by="genotype",
                         suffix = c('_blues',"_yieldrank"))

cor(fulldataset$estimate_blues,fulldataset$estimate_yieldrank,use = "pairwise.complete.obs")

fulldataset %>% ggplot(aes(x=estimate_blues,y=estimate_yieldrank)) + 
  theme_bw() + geom_text(aes(label = genotype)) + geom_smooth(method = 'lm',se=F) +
  xlab('Genotypic means ton/ha') +ylab('Genotypic worth values')

## Metric data

# genotypic variance

fulldataset$yieldrank_weights <- 1/fulldataset$SE_yieldrank^2

fulldataset$yieldmetric_weights <- 1/fulldataset$SE_blues^2

fulldataset$genotype <- factor(fulldataset$genotype)

blups <- asreml(fixed = estimate_blues ~ 1, random = ~ genotype , weights = yieldmetric_weights,
                family = (asr_gaussian(dispersion = 1)),
                data = fulldataset,trace=F,na.action = na.method(x='omit'))

summary(blups)$varcomp

# heritabilitey

qH2(blups)

## Ranking data

fulldataset$genotype <- factor(fulldataset$genotype)

blups <- asreml(fixed = estimate_yieldrank ~ 1, random = ~ genotype , weights = yieldrank_weights,
                family = (asr_gaussian(dispersion = 1)),
                data = fulldataset[-1,],trace=F,na.action = na.method(x='omit'))

summary(blups)$varcomp

# heritabilitey

qH2(blups)

