
# Data processing thurstonian and PL experiments
# Reading files collected from the servers and adding mixed models
# Hugo Dorado
# 6-18/2024

library(tidyverse)
library(stringi)
library(asreml)
library(PlackettLuce)
library(emmeans)
library(HelpersMG)
library(complexlm)
library(asreml)
library(lme4)

# Data astringi# Data analyisis - sample size

setwd('/home/joost/Hugo_Projects/Sim_est_genotypic_values/Genotypic_values_performance3/')

nam.fls <- list.files()

paths <- list.files(full.names = T)

ls_estimations <- lapply(paths, readRDS)

names(ls_estimations) <- nam.fls

#-------------------------------------------------------------------------------

# Check that all simulated tests are available


final_results <- do.call(rbind,lapply(ls_estimations,function(wa){wa$pars})) %>% 
  dplyr::select(sigma.genotype,sigma.plot,n.rep.genotype.env) %>% 
  group_by(sigma.genotype,sigma.plot,n.rep.genotype.env) %>%
  summarise(n=n())

table(final_results$n) # 6*6*6


#dataset <- out$dataset

ls_estimations <- lapply(ls_estimations, function(out){

  out$mod.blues <- lmer(trait.value.obs~genotype+(1|farm.vec),data=out$dataset)

  out$mod.blups <- lmer(trait.value.obs~(1|genotype)+(1|farm.vec),data=out$dataset)
  
  out
}
)

#names(ls_estimations) <- nam.fls

saveRDS(ls_estimations,'/home/joost/Hugo_Projects/Sim_est_genotypic_values/Processed_data/output_ls_estimations2.rds')

#saveRDS(ls_estimations,'Results/se_simulations.rds')




