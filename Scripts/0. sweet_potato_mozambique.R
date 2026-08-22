
library(tidyverse)
library(PlackettLuce)
library(gosset)

setwd('Paper_1_Estimation_of_genetic_parameters')

source('Scripts/0. Ranking_pars_FUN.R')

sweet_potato <- read.csv('Processed_data/On-farm-datasets/sweet_potato.csv',row.names = 1)

sweet_potato <- data.frame(id=1:nrow(sweet_potato),sweet_potato)
           
names(sweet_potato)[7:8] <- paste0('rnk',names(sweet_potato)[7:8] )

names(sweet_potato) <- str_replace_all(names(sweet_potato) ,c('_a'='_A',"_b"='_B','_c'='_C'))

var_metric <- sweet_potato %>% dplyr::select(id,package_item_A:package_item_C,crw_A:crw_C)

names(var_metric) <- gsub("package_item", 'packageitem',names(var_metric))
library(stringr)



var_metric <- var_metric %>% 
  pivot_longer(
    cols = !id, 
    names_to = c(".value", "plot"), 
    names_sep = "_", 
    values_drop_na = F
  )

table(sweet_potato$rnkcrw_pos==sweet_potato$rnkcrw_neg)

var_ranking <- sweet_potato %>% group_by(id) %>% 
  mutate(rnkcrw_mid=  c("A","B","C")[!( c("A","B","C") %in% c(rnkcrw_pos,rnkcrw_neg))]  ) %>% 
  select(id,rnkcrw_pos,rnkcrw_mid,rnkcrw_neg) %>% pivot_longer(!id) %>%
  separate(name,c('var', 'rank'),"_") %>% mutate(rank=str_replace_all(rank,c('pos'='1','mid'='2','neg'='3')))





full_sweet_potato <- var_metric %>% left_join(sweet_potato[c('id','project_id','project_Code')],by=c('id')) %>%   
  left_join(var_ranking,by = c('id','plot'='value')) %>%
  select(-var) 
  
full_sweet_potato <- full_sweet_potato[,c('id','packageitem','project_id' ,'project_Code','crw','rank')]

full_sweet_potato$crw

full_sweet_potato <- full_sweet_potato %>% group_by(id) %>% mutate(rank2=round(rank(-crw,ties.method = 'first'),0))

table(full_sweet_potato$rank,full_sweet_potato$rank2)

names(full_sweet_potato) <- c('id','genotype', 'projectid','projectcode', 'yield','rank', 'rank2')

write.csv(full_sweet_potato,'Processed_data/On-farm-datasets/sweet_potato_processed.csv')

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


