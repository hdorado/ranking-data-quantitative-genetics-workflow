
# Ordered use the metric data reported in the dataset and calculate a ranking 
# The actual ranking reported isn't account for the analysis.

library(lme4)
library(emmeans)
library(stringr)
library(tidyverse)
library(ggpmisc)

setwd('Paper_1_Estimation_of_genetic_parameters/')

nams <- list.files('Results/Observed_data_outputs_ranked/')

observed_dataset <- lapply(list.files('Results/Observed_data_outputs_ranked/',full.names = T),readRDS)

names(observed_dataset) <- nams

# exploratory data analysis

# Cassava----------------------------------------------------------------------

cassava_datasets <-  do.call(rbind,
  lapply(observed_dataset[grep('cassava',nams)],function(obs){
  
    obs$dataset 
  })
)

table(cassava_datasets$genotype,cassava_datasets$Zone)

cassava_datasets$planting_date <- as.Date(cassava_datasets$planting_date)

cassava_datasets %>% group_by(Zone) %>% summarise(Min = min(planting_date),
                                                  Max = max(planting_date)
                                                  )


cassava_datasets <- cassava_datasets %>% group_by(farm) %>% group_split()  %>%  
  map( ~ .x %>%  mutate(order = base::rank(-1*yieldton_ha,ties.method = 'random'))) %>% bind_rows

table(cassava_datasets$rank_freshrootsyield,cassava_datasets$order,cassava_datasets$Zone)
  
table(cassava_datasets$genotype,cassava_datasets$rank_freshrootsyield)

data.frame(
cassava_datasets %>% group_by(Zone,genotype) %>% 
  summarise(rank_1=sum(rank_freshrootsyield==1),
                          rank_2=sum(rank_freshrootsyield==2),
                                        rank_3=sum(rank_freshrootsyield==3))

)

# Groundnut-------------------------------------------------------------------


groundnut_datasets <-  do.call(rbind,
                             lapply(observed_dataset[grep('groundnut',nams)],function(obs){
                               
                               obs$dataset 
                             })
)

table(groundnut_datasets$genotype,paste0(groundnut_datasets$ecoregion,groundnut_datasets$year))

groundnut_datasets$plantingdate <- as.Date(groundnut_datasets$plantingdate)

groundnut_datasets %>% group_by(ecoregion,year) %>% summarise(Min = min(plantingdate,na.rm = T),
                                                  Max = max(plantingdate,na.rm=T)
)


groundnut_datasets <- groundnut_datasets %>% group_by(farm) %>% group_split()  %>%  
  map( ~ .x %>%  mutate(order = base::rank(-1*yieldton_ha,ties.method = 'random'))) %>% bind_rows

table(groundnut_datasets$rank_yield,groundnut_datasets$order,paste0(groundnut_datasets$ecoregion,groundnut_datasets$year))


data.frame(
  groundnut_datasets %>% group_by(ecoregion,year,genotype) %>% 
    summarise(rank_1=sum(rank_yield==1),
              rank_2=sum(rank_yield==2),
              rank_3=sum(rank_yield==3))
  
)

# maize-------------------------------------------------------------------


maize_datasets <-  do.call(rbind,
                               lapply(observed_dataset[grep('maize',nams)],function(obs){
                                 
                                 obs$dataset 
                               })
)

table(maize_datasets$genotype,paste0(maize_datasets$zone))

maize_datasets$plantingdate <- as.Date(maize_datasets$plantingdate)

maize_datasets %>% group_by(zone) %>% summarise(Min = min(plantingdate,na.rm = T),
                                                              Max = max(plantingdate,na.rm=T)
)


maize_datasets <- maize_datasets %>% group_by(farm) %>% group_split()  %>%  
  map( ~ .x %>%  mutate(order = base::rank(-1*yield_farmer,ties.method = 'random'))) %>% bind_rows

table(maize_datasets$rank_yield,maize_datasets$order,paste0(maize_datasets$zone))


data.frame(
  maize_datasets %>% group_by(zone,genotype) %>% 
    summarise(rank_1=sum(rank_yield==1),
              rank_2=sum(rank_yield==2),
              rank_3=sum(rank_yield==3))
  
)






# -------------------------------Cassava ---------------------------------------



yield_nam <- 'yieldton_ha'
genotype_nam <- 'genotype'
farm_nam <- 'farm'

cassava_results <-  do.call(rbind,lapply(observed_dataset[grep('cassava',nams)],function(exp){
  
    dataset <-  exp$dataset
    
    selec_var <- c(yield_nam,genotype_nam,farm_nam)
    
    ds <- dataset[selec_var]
    
    names(ds) <- c('yield','genotype','farm')
    
    model <- lmer(yield~genotype+(1|farm),data = ds)
    
    sigma_vals <- data.frame(VarCorr(model))
    
    sigma_e <- sigma_vals[sigma_vals$grp=='Residual',]$sdcor
    
    genotypic_means <- emmeans(model,'genotype')
    
    genotypic_means <- data.frame(genotypic_means)[c('genotype',  'emmean', 'SE' )]
    
    rank_analysis <- rbind(data.frame(genotype=exp$ref,est_means=0,model='Th'),
                           data.frame(genotype = names(exp$th_mus$par),
                                      est_means=exp$th_mus$par,model='Th'),
                           data.frame(genotype=exp$ref,est_means=0,model='PL'),
        data.frame(genotype = names(exp$log_worths$par),est_means=exp$log_worths$par,
                   model='PL'))
    
    rank_analysis$genotype <- gsub('^mu','',rank_analysis$genotype)
    
    rank_analysis <- rank_analysis %>% mutate(est_means_sca = genotypic_means[genotypic_means$genotype == gsub('^mu','',exp$ref) , ]$emmean + est_means*sigma_e)
    
    rank_analysis <- rank_analysis %>% left_join(genotypic_means,by='genotype')
    
    data.frame(Location = exp$loc,rank_analysis)
  }))


ga <- cassava_results %>% ggplot(aes(x=emmean ,y= est_means_sca)) + geom_point(aes(colour = model)) +
  theme_bw() + labs(color = 'Model') +  geom_abline(intercept = 0,slope=1) + 
  stat_poly_eq(aes(label =  paste(..eq.label.., ..rr.label.., 'rho','`=`',round(sqrt(r.squared),2),sep = "~~~~"),colour = model),formula = y ~ x,parse = TRUE,, size = 2)+
  geom_smooth(method = 'lm',aes(colour=model),se=F,linewidth=0.7)+xlab('Estimated genotypic means for yield t/ha')+ theme(plot.title = element_text(hjust = 0.5)) +
  ylab('Estimated genotypic means for yield ranked') +facet_wrap(~Location,scales='free') +ggtitle(label ='Cassava')

ggsave('Results/cassava_genotypic_means_ranked.png',ga,height = 8,width = 9.5)

# -------------------------------Groundnut---------------------------------------

yield_nam <- 'yieldton_ha'
genotype_nam <- 'genotype'
farm_nam <- 'farm'

groundnut_results <-  do.call(rbind,lapply(observed_dataset[grep('groundnut',nams)],function(exp){
  
  dataset <-  exp$dataset
  
  selec_var <- c(yield_nam,genotype_nam,farm_nam)
  
  ds <- dataset[selec_var]
  
  names(ds) <- c('yield','genotype','farm')
  
  model <- lmer(yield~genotype+(1|farm),data = ds)
  
  sigma_vals <- data.frame(VarCorr(model))
  
  sigma_e <- sigma_vals[sigma_vals$grp=='Residual',]$sdcor
  
  genotypic_means <- emmeans(model,'genotype')
  
  genotypic_means <- data.frame(genotypic_means)[c('genotype',  'emmean', 'SE' )]
  
  rank_analysis <- rbind(data.frame(genotype=exp$ref,est_means=0,model='Th'),
                         data.frame(genotype = names(exp$th_mus$par),
                                    est_means=exp$th_mus$par,model='Th'),
                         data.frame(genotype=exp$ref,est_means=0,model='PL'),
                         data.frame(genotype = names(exp$log_worths$par),est_means=exp$log_worths$par,
                                    model='PL'))
  
  rank_analysis$genotype <- gsub('^mu','',rank_analysis$genotype)
  
  rank_analysis <- rank_analysis %>% mutate(est_means_sca = genotypic_means[genotypic_means$genotype == gsub('^mu','',exp$ref) , ]$emmean + est_means*sigma_e)
  
  rank_analysis <- rank_analysis %>% left_join(genotypic_means,by='genotype')
  
  data.frame(Location = unique(exp$dataset$ecoregion), year = unique(exp$dataset$year),rank_analysis)
}))


ga <- groundnut_results %>% ggplot(aes(x=emmean ,y= est_means_sca)) + geom_point(aes(colour = model)) +
  theme_bw() + labs(color = 'Model') +  geom_abline(intercept = 0,slope=1) + 
  stat_poly_eq(aes(label =  paste(..eq.label.., ..rr.label.., 'rho','`=`',round(sqrt(r.squared),2),sep = "~~~~"),colour = model),formula = y ~ x,parse = TRUE,, size = 2)+
  geom_smooth(method = 'lm',aes(colour=model),se=F,linewidth=0.7)+xlab('Estimated genotypic means for yield ton/ha')+ theme(plot.title = element_text(hjust = 0.5)) +
  ylab('Estimated genotypic means for yield ranked') +facet_wrap(~Location+year,scales='free') +ggtitle(label ='Groundnut')

ggsave('Results/groundnut_genotypic_means_ranked.png',ga,height = 8,width = 9.5)


# -------------------------------Maize---------------------------------------

yield_nam <- 'yield_farmer'
genotype_nam <- 'genotype'
farm_nam <- 'farm'

maize_results <-  do.call(rbind,lapply(observed_dataset[grep('maize',nams)],function(exp){
  
  dataset <-  exp$dataset
  
  selec_var <- c(yield_nam,genotype_nam,farm_nam)
  
  ds <- dataset[selec_var]
  
  names(ds) <- c('yield','genotype','farm')
  
  ds$farm <- factor(ds$farm)
  
  model <- lmer(yield~genotype+(1|farm),data = ds)
  
  sigma_vals <- data.frame(VarCorr(model))
  
  sigma_e <- sigma_vals[sigma_vals$grp=='Residual',]$sdcor
  
  genotypic_means <- emmeans(model,'genotype')
  
  genotypic_means <- data.frame(genotypic_means)[c('genotype',  'emmean', 'SE' )]
  
  rank_analysis <- rbind(data.frame(genotype=exp$ref,est_means=0,model='Th'),
                         data.frame(genotype = names(exp$th_mus$par),
                                    est_means=exp$th_mus$par,model='Th'),
                         data.frame(genotype=exp$ref,est_means=0,model='PL'),
                         data.frame(genotype = names(exp$log_worths$par),est_means=exp$log_worths$par,
                                    model='PL'))
  
  rank_analysis$genotype <- gsub('^mu','',rank_analysis$genotype)
  
  rank_analysis <- rank_analysis %>% mutate(est_means_sca = genotypic_means[genotypic_means$genotype == gsub('^mu','',exp$ref) , ]$emmean + est_means*sigma_e)
  
  rank_analysis <- rank_analysis %>% left_join(genotypic_means,by='genotype')
  
  data.frame(Location = unique(exp$dataset$zone), rank_analysis)
}))


ga <- maize_results %>% ggplot(aes(x=emmean ,y= est_means_sca)) + geom_point(aes(colour = model)) +
  theme_bw() + labs(color = 'Model') +  geom_abline(intercept = 0,slope=1) + 
  stat_poly_eq(aes(label =  paste(..eq.label.., ..rr.label.., 'rho','`=`',
                                  round(sqrt(r.squared),2),sep = "~~~~"),
                   colour = model),formula = y ~ x,parse = TRUE, size = 2)+
  geom_smooth(method = 'lm',aes(colour=model),se=F,linewidth=0.7)+xlab('Estimated genotypic means for yield ton/ha')+ theme(plot.title = element_text(hjust = 0.5)) +
  ylab('Estimated genotypic means for yield ranked') +facet_wrap(~Location,scales='free') + 
  ggtitle(label ='Maize')


ggsave('Results/maize_genotypic_means_ranked.png',ga,height = 8,width = 9.5)


# calculate, statistics of overprefered varieties and tables with number of observations by varieties

