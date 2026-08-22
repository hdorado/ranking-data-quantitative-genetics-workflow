
# Final results
# Hugo Andres Dorado
# 2025-03-28

library(tidyverse)
library(ggpmisc)
library(ggpubr)
library(ggiraph)
library(ggrepel)
library(gtools)
library(ggpp)

setwd('Paper_1_Estimation_of_genetic_parameters/')

custom_colors2 <- c("Th" = "#00CCBE",  # Deep Red
                    "PL" = "#FEA800",
                    "MM"= "#F8766D") 

#-----------------------------------------------------------

data_30_exmple <- read.csv('Results/Plot_inputs/dataset_summary_out_30_ref_case.csv')




data_30_exmple <- data_30_exmple %>% group_by(iter,sigma_g,sigma_e,ng,model) %>% 
  mutate(ref_est_mean = est_means[genotype=="G1"])

data_30_exmple <- data_30_exmple %>% mutate(est_means = ifelse(model=="MM",
                                                               (est_means-ref_est_mean)/sigma_e,est_means),
                                            SE = ifelse(model=="MM",SE/sigma_e,SE))


data_30_exmple$model <- factor(data_30_exmple$model,c('Th','PL','MM'))

plot_val <- data_30_exmple %>% filter(iter==3 & model=='Th') %>% 
  ggplot(aes(x=scl_true_means,y=est_means)) + 
  # geom_text_repel(  #+ geom_text(aes(label = genotype))
  #  aes(label=genotype,color=model),
  #  box.padding = unit(0.2, "lines"),size=3,max.overlaps = 25
  #)+ 
  #facet_wrap(~iter,scales='free') + 
  geom_abline(intercept = 0,slope = 1)+ 
  #geom_smooth(aes(color=model),se=F,method = 'lm')+  
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th rank-based", PL="PL rank-based",MM="MM metric-based"))+
  ylab('Estimated genotypic means') + geom_point(aes(color=model))+
  xlab("True genotypic means") + 
  scale_x_continuous(breaks = c(-1,0,1, 2),limits = c(-1,2.3))+
  scale_y_continuous(breaks = c(-1,0,1, 2, 3)) +  
  geom_errorbar(aes(ymax = est_means+SE,ymin=est_means-SE,colour = model),
                width=0.05, alpha=0.7, linewidth=0.5)+theme_bw() +theme(legend.position = "none")

ggsave('Results/Paper_figures_V5/worth_values_estimated_vs_true_EUCARPIA.png',plot_val,height = 3.2,width = 4)


plot_serror <- data_30_exmple %>% filter(iter==3 ,genotype!='G1',model=="Th")  %>% 
  ggplot(aes(x=est_means,y=SE,colour = model)) +geom_point() +
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th rank-based", PL="PL rank-based",MM="MM metric-based"))+
  ylab('Estimated genotypic means') + geom_point(aes(color=model))+ theme_bw()+
  xlab('Estimated genotypic means') + ylab('Std. error of estimated genotypic means')+theme(legend.position = "none")


ggsave('Results/Paper_figures_V5/stand_err_vs_genVals_EUCARPIA.png',plot_serror,height = 3.2,width = 4)

# Scaled case

data_exmple_scl <- data_30_exmple %>% filter(iter==3 ,genotype!='G1')

data_exmple_scl <- data_exmple_scl %>% mutate(est_means2= est_means, est_means2= ifelse(model=="PL", 
                                                                                        est_means*0.8372322,
                                                                                        est_means ),SE_2=SE,SE_2=ifelse(model=="PL", 
                                                                                                                        SE_2*0.8372322,
                                                                                                                        SE_2 ))


data_exmple_scl %>% ggplot(aes(x=est_means2,y=scl_true_means)) + geom_point(aes(color=model))


data_exmple_scl$model


data_exmple_scl$model <- factor(data_exmple_scl$model,levels = c("Th","PL","MM"))

plot_serror <- data_exmple_scl %>% filter(iter==3 ,genotype!='G1',model %in% c("Th", "MM"))  %>% 
  ggplot(aes(x=est_means2,y=SE_2,colour = model)) +geom_point() +
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th rank-based",PL="PL rank-based (scl)",MM= "MM quant.-based"))+ theme_bw()+
  xlab('Estimated genotypic means') + ylab('Std. error of estimated genotypic means')


ggsave('Results/Paper_figures_V5/stand_err_vs_genVals2_EUCARPIA.png',plot_serror,
       height = 3.2,width = 5)


se_dataset <- read.csv("Results/Plot_inputs/summary_analysis_simulation_SE.csv")

se_dataset$genotype <- factor(se_dataset$genotype,levels=mixedsort((unique(se_dataset$genotype))))

ge_val_out <- se_dataset%>% filter(model!="MM" &  genotype != "G1")  %>% ggplot(aes(x=genotype,y=est_means)) +  
  scale_color_manual(name= "Model",values=custom_colors2,
                     labels = c(Th="Th rank-based", PL="PL rank-based",MM="MM metric-based"))+
  geom_boxplot(aes(color=model))+ylab("Estimated genotypic values") +theme_bw() +xlab('Genotype')

#ggsave('Results/Paper_figures_V5/appendix_geno_values_dist.png',ge_val_out,height = 3.2,width = 8)

rmvIt <- se_dataset$iter[which(se_dataset$model=='PL'& se_dataset$est_means>4.3 |se_dataset$model=='PL'& se_dataset$est_means < -3 )]


#  Adjust scale

scale_pl_th_ajs <- se_dataset %>% select(iter,genotype,model,est_means) %>% 
  pivot_wider(names_from = model,values_from = est_means) %>% 
  group_by() %>% mutate(slope=coefficients(lm(Th~-1+PL))[1])


se_dataset <- se_dataset %>% left_join(scale_pl_th_ajs[c('iter','genotype','slope')],by=c('iter', 'genotype'))

se_dataset <- se_dataset %>% filter(model!="MM"  ) %>% mutate(est_means2 = ifelse(model=="PL",est_means*slope ,est_means),
                                                              SE2 = ifelse(model=="PL",SE*slope ,SE))

se_dataset_summ <- se_dataset %>% filter(genotype!='G1' ) %>% filter(model!="MM"  )  %>% filter(iter!=25  ) %>%
  group_by(ng,genotype,sigma_e,sigma_g,model) %>% 
  summarise(mean_SE=mean(SE2),SE_mean=sd(SE2)/sqrt(n()),sd_estimate=sd(est_means2),
            median_SE=median(SE2),sd_me_estimate=mad(est_means2),n=n())



#se_dataset_summ <- mutate(correct_mean = ifelse(model=='PL',mean_SE,),correct_se=ifelse(model=='PL',,))

se_dataset_summ$model <- factor(se_dataset_summ$model,levels = c("Th","PL"))

se_val_plot <- se_dataset_summ %>%  ggplot(aes(x=sd_estimate,y=mean_SE)) + geom_point(aes(color=model)) +
  geom_abline(intercept=0,slope=1)+theme_bw() +  
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th rank-based", PL="PL rank-based (scl)",MM="MM quant.-based")) +
  # geom_errorbar(aes(ymax = mean_SE+SE_mean,ymin =mean_SE-SE_mean ,color=model),width=0.005, alpha=0.7, linewidth=0.5)+
  xlab("Std. dev. of estimated genotypic means") + ylab("Mean of std. errors of estimated \ngenotypic means")+
  ylim(0.35,0.5)+xlim(0.35,0.5)

se_val_plot

#ggsave('Results/Paper_figures_V5/standr_erro_val2.png',se_val_plot,height = 3.2,width =5)

data_30_exmple$genotype <- factor(data_30_exmple$genotype,levels=mixedsort((unique(data_30_exmple$genotype))))

library(stringi)
library(stringr)

one_run <- data_30_exmple %>% filter(iter==3 & model!='MM') %>% 
  dplyr::select(genotype,est_means,SE,quasiSE) %>% 
  pivot_longer(!c(iter:est_means),names_to = 'variance_type',values_to = 'value') %>%
  
  mutate(variance_type=str_replace_all(variance_type,c( 'quasiSE'="Based on quasi-standard errors",'SE'='Based on conventional standard errors')))



est_true <- one_run %>% 
  ggplot(aes(x=genotype,y=est_means,color=model)) + geom_point(position=position_dodge(width=0.5)) + 
  scale_color_manual(name= "Model",values=custom_colors2,
                     labels = c(Th="Th", PL="PL", MM="MM"))+
  ylab('Estimated genotypic values') +
  xlab("True genotypic values") + 
  #scale_x_continuous(breaks = c(-1,0,1, 2),limits = c(-1,2.3))+
  #scale_y_continuous(breaks = c(-1,0,1, 2, 3)) +  
  geom_errorbar(aes(ymax = est_means + 1.96 * value, ymin = est_means - 1.96 * value),
                width=0.3, alpha=0.7, linewidth=0.6, position=position_dodge(width=0.5))+
  facet_grid(variance_type~.)+
  theme_bw()  


#ggsave('Results/Paper_figures_V5/ques_SE.png',est_true,height = 5,width = 8)

# missing values and indeterminate cases analysis

extreme_cases <- read.csv('Results/Plot_inputs/outlayers_inderterminates.csv')

extreme_cases$model <- factor(extreme_cases$model ,levels = c("Th","PL"))

4*14*11*30

extreme_cases$ng2 <- factor(extreme_cases$ng,levels = c(3,10,30,100),
                            labels = c('r = 3', 
                                       # expression(r[g] == 5), 
                                       'r = 10', 
                                       'r = 30', 
                                       'r = 100'))

custom_colors2 <- c("Th" = "#00CCBE",  # Deep Red
                    "PL" = "#FEA800",
                    "MM"= "#F8766D") 



extreme_cases <- extreme_cases %>% pivot_longer(-c(X:n,ng2),names_to='cases',values_to ='percentaje')

extreme_cases$cases2 <- "Extreme estimated genotypic values"

extreme_cases$cases2[extreme_cases$cases=='perc_sing'] <- "Undefined estimated standard errors"

extreme_cases<- extreme_cases[extreme_cases$cases2 != 'Extreme estimated genotypic values',]

extreme_cases_plot <- extreme_cases %>% filter(model=='Th')%>%ggplot(aes(x=sigma_std,y=percentaje)) +
  geom_point(aes(color=model)) +geom_line(aes(color=model)) + 
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th rank-based", PL="PL rank-based"))+
  facet_grid(.~ng2)+
  theme_bw()+
  ylab("Percentage of trials with at least one \nundefined estimate of variance") + xlab(expression(paste("True genetic sd ", (sigma[g]), "")))+theme(legend.position= 'none')


ggsave('Results/Paper_figures_V5/extreme_cases_plot_EUCARPIA.png',extreme_cases_plot,height = 3.2,width = 9)


#--------------------------


# equivalence

PL_th_equivalence <- read.csv('Results/Plot_inputs/equva_factor_ds_val.csv')

PL_th_equivalence <-PL_th_equivalence[round(PL_th_equivalence$Ratio,8) %in% c(0.18181818,0.71428571 ,0.50000000,2.25000000,4.50000000),] 

PL_th_equivalence$Ratio <- factor(round(PL_th_equivalence$Ratio,2))


PL_th_equivalence <- PL_th_equivalence %>% ggplot(aes(x=ng,y=slope_mean))+ 
  geom_point(aes(color=factor(Ratio))) + geom_line(aes(color=factor(Ratio))) +
  theme_bw()+ylab("Regression slope (Thurstonian vs. Plackett-Luce estimates)")+
  labs(color= expression(paste("True genetic sd ", sigma[g], "")))+
  xlab(expression(paste("Number of replicates", (r), "")))

#ggsave('Results/Paper_figures_V5/PL_th_equivalence.png',PL_th_equivalence,height = 4.5 ,width =6.5 )


# Magical factor

mag_fact2 <- read.csv('Results/Plot_inputs/mag_fact2.csv')

mag_fact2b <- mag_fact2[round(mag_fact2$ratio,8) %in% c(0.18181818,0.71428571,0.5,2.25,4.50000000),] 

mag_fact2b <- mag_fact2b %>% mutate(ratio2=round(ratio,2))


mag_fact2b <- mag_fact2b %>% group_by(ratio2,model,ng) %>% 
  summarise(medn=median(mu_mag_fact_median_qseQse,na.rm = T)) 


mag_fact2b <- mag_fact2b %>% mutate(model2 = ifelse(model=="PL", "Plackett-luce","Thurstonian"))

mag_fact2b$model2 <- factor(mag_fact2b$model2,levels = c("Thurstonian","Plackett-luce"))

ggmagFac <- mag_fact2b  %>%   #[round(mag_fact2$ratio,8) %in% c(0.18181818,0.71428571,0.5,2.25,4.50000000),] %>% 
  ggplot(aes(x=ng,y=medn)) +geom_point(aes(color=factor(ratio2))) + 
  labs(color= expression(paste("True genetic sd ", (sigma[g]), "")))+
  geom_line(aes(color=factor(ratio2)))+theme_bw()+ylim(0,5)+ 
  facet_grid(.~model2) + 
  ylab('Ratio of medians of quasi-standard errors:\nrank- vs metric-based')+
  xlab(expression(paste("Number of replicates", (r), "")))


#ggsave("Results/Paper_figures_V5/magica_factor.png",ggmagFac,height = 3.5,width = 8)


#mag_fact2 <- mag_fact2 %>% filter(ng!=5) 



mag_fact2b %>% ggplot(aes(x=ng,y=mu_mag_fact_mean_SEQse))+ 
  geom_point(aes(color=factor(ratio2))) + 
  geom_line(aes(color=factor(ratio2)))+theme_bw()+#ylim(c(0,5))+#geom_smooth(aes(color=factor(ratio)),se=F)+
  ylab("median(qSE_rank)/(median(metric_SE)/sigma_e)")+
  xlab(expression(paste("True genetic sd (", widetilde(sigma)[g], ")")))


#--------------------------Predictive performance-------------------------------

index_summary <- read.csv("Results/Plot_inputs/index_summary_RMSE_CCOR_R2_V3_fulldataset.csv")

performance_Values <- index_summary %>% pivot_longer(-c(iter,ng,sigma_e,sigma_g,model,Ratio),
                                                     names_to = 'variable',values_to = 'value') %>%
  group_by( ng, sigma_e, sigma_g, model, Ratio, variable) %>% 
  summarise(mu_value=mean(value),SE=1.96*sd(value)/sqrt(n())) #Na.rm to include missing values

index_summary %>% filter(sigma_e==2.0 & sigma_g== 0.1 & model=='MM') %>% group_by(ng) %>% summarise(mean(cor))


custom_colors2 <- c("Th" = "#00CCBE",  # Deep Red
                    "PL" = "#FEA800",
                    "MM"= "#F8766D") 

performance_Values$variable <- factor(performance_Values$variable,levels = c("cor","Bias","RMSE"))

levels(performance_Values$variable) <- c(expression(paste("Kendall correlation coef. (", tau,")")),
                                         #expression(r[g]==50),
                                         "Bias","RMSE")

performance_Values$ng2 <- factor(performance_Values$ng,levels = c(3,5,10,30,100),
                                 labels = c(expression(r == 3), #r[g] 
                                            expression(r == 5), 
                                            expression(r == 10), 
                                            expression(r == 30), 
                                            expression(r == 100)))

#-------------------------------------------------------------------------------
# GV - Person correlation + kendall-index vs radio of variation 

performance_Values$model <- factor(performance_Values$model,levels = c("Th","PL","MM"))

pv <- performance_Values %>% 
  filter(ng %in% c(3,10,30,100) & model %in% c('Th','MM') &
           (variable %in% c('Bias','RMSE',
                            "paste(\"Kendall correlation coef. (\", tau, \")\")")
           )) %>%
  ggplot(aes(x=Ratio,y=mu_value,group=interaction(model,ng))) +
  geom_point(aes(colour = model),alpha=0.6,size=1.3) +
  geom_vline(xintercept = 0.8888889,linetype= 'dashed')+
  geom_errorbar(aes(ymax = mu_value+SE,ymin=mu_value-SE,colour = model),
                width=0.05, alpha=0.8, linewidth=0.5)+
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th rank-based", PL="PL rank-based (scl)",MM="MM quant.-based"))+
  facet_grid(variable~ng2,labeller =label_parsed,scales='free')+theme_bw() + #,switch="both"
  ylab('') + geom_line(aes(colour = model),alpha=0.7)+
  xlab(expression(paste("True genetic sd ", (sigma[g]), "")))


pv

ggsave('Results/Paper_figures_V5/Figure_1_GV_rho_c_index_vs_ratio_EUCARPIA.png',pv,height = 6,width = 9)


#-------------------------------------------------------------------------------

# Estimation genotypic variances new version


genotypic_se_g0 <- read.csv('Results/Plot_inputs/estimation_genotypic_variance.csv')  # UPDATE THIS FILE FROM THE SOURCE

head(genotypic_se_g0)

# paste valid cases

# full_ds_valid_cases <- read.csv('Results/full_experiment_valid_novalidcases.csv')

# genotypic_se_g0 <- genotypic_se_g0 %>% mutate( Experiment0 = paste(iter,sigma_g,sigma_e,ng,sep='-') )

# genotypic_se_g0 <-  genotypic_se_g0 %>% left_join(unique(full_ds_valid_cases[c('Experiment0','Exp_valid')]))

# write.csv(genotypic_se_g0,'Results/Plot_inputs/estimation_genotypic_variance.csv',row.names=F)

genotypic_se_g0 <- read.csv('Results/Plot_inputs/estimation_genotypic_variance.csv')  # UPDATE THIS FILE FROM THE SOURCE



custom_colors2 <- c("Th" = "#00CCBE",  # Deep Red
                    "PL" = "#FEA800",
                    "MM"= "#F8766D") 


genotypic_se_g <- genotypic_se_g0 %>% mutate(true_sigma_sc=sigma_g/sigma_e) %>% 
  group_by(ng, true_sigma_sc, model) %>% filter(!is.na(sg_g_estimate_scl2)) %>% 
  summarise(mu_sd_g_wQV = mean(sg_g_estimate_scl2,na.rm=T),
            SE_mu_sd_g_wQV = sd(sg_g_estimate_scl2,na.rm=T)/sqrt(n()) ,n=n())



genotypic_se_g <- genotypic_se_g %>% filter(ng %in% c(3,10,30,100)  ) 

genotypic_se_g <- genotypic_se_g %>% mutate(ng2=factor(paste0('r=',ng),
                                                       levels=c('r=3','r=10',
                                                                'r=30','r=100')))


genotypic_se_g$ng3 <- factor(genotypic_se_g$ng,levels = c(3,10,30,100),
                             labels = c(expression(r == 3), 
                                        expression(r == 10), 
                                        expression(r == 30), 
                                        expression(r == 100)))

genotypic_se_g$model <- factor(genotypic_se_g$model,levels = c("Th","PL",'MM'))

gene_vars <- genotypic_se_g %>% filter(model!="PL",ng!=10) %>% ggplot(aes(y=mu_sd_g_wQV ,x=true_sigma_sc,color=model ))+ 
  geom_point(alpha=0.3,show.legend = F)+
  # geom_errorbar(aes(ymin=mu_sd_g_wQV_scaled-1.96*SE_mu_sd_g_wQV_scaled,ymax=mu_sd_g_wQV_scaled+1.96*SE_mu_sd_g_wQV_scaled,colour = model),width=0.1, alpha=0.7, linewidth=0.3)+
  geom_line(aes(color=model),alpha=0.8,linewidth=0.5)+
  geom_abline(intercept =  0,slope=1,linetype=1)+
  #scale_shape_manual(name='Number \nof replicates',values = c(19,17,15,3,8,7))+
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th rank-based", PL="PL rank-based (scl)",MM="MM quant.-based"))+
  scale_linetype(name=expression(atop("Residual ",paste('stad. dev. ',sigma[e]))))+
  xlab(expression(paste("True genetic sd ", (sigma[g]), "")))+
  geom_vline(xintercept = 0.8,linetype =2)+
  #facet_wrap(~ng)+
  ylab(expression(paste("Estimated genetic sd ", (hat(sigma[g])), "")))+#+ggtitle(label = 'include out')+
  facet_wrap(ng3~.,labeller =label_parsed)+theme_bw()+ylim(0,5.5)

ggsave("Results/Paper_figures_V5/genotypic_var_estm_2_EUCARPIA.png",gene_vars,width =8.5 ,height = 2.7)

# Estimation of heritabilities new version

genotypic_se_g0 <- read.csv('Results/Plot_inputs/estimation_genotypic_variance.csv')  # UPDATE THIS FILE FROM THE SOURCE

genotypic_se_g0 <- genotypic_se_g0 %>% filter(Exp_valid=="Valid")

genotypic_se_g0$heritability_fnl[!is.na(genotypic_se_g0$heritability_fnl) & (genotypic_se_g0$heritability_fnl<=0)] <- NA

genotypic_se_g0$h2_metric_cullis[!is.na(genotypic_se_g0$h2_metric_cullis) & (genotypic_se_g0$h2_metric_cullis<=0)] <- NA

genotypic_se_g0$heritability_cullisEstQV[!is.na(genotypic_se_g0$heritability_cullisEstQV) & (genotypic_se_g0$heritability_cullisEstQV<=0)] <- NA

genotypic_se_g0 <-  genotypic_se_g0 %>% filter(Exp_valid=="Valid")


heritability <- genotypic_se_g0 %>% group_by(ng, sigma_e, sigma_g ,model) %>% 
  summarise(mu_cullis_ref=mean(h2_metric_cullis,na.rm = T),
            mu_cullis=mean(heritability_fnl,na.rm=T),
            mu_cullis2=mean(heritability_cullisEstQV,na.rm = T))



heritability %>% ggplot(aes(x=mu_cullis_ref ,y=mu_cullis2) ) + 
  facet_wrap(~ng) + geom_point(aes(color=model))+
  geom_abline(intercept = 0,slope=1)

heritability$ng3 <-  factor(heritability$ng,levels = c(3,10,30,100),
                            labels = c(expression(r == 3), 
                                       expression(r == 10), 
                                       expression(r == 30), 
                                       expression(r == 100)))

heritability <- heritability %>% ungroup()%>% select(ng,ng3,model,mu_cullis_ref,mu_cullis2) %>% 
  group_by(ng,model,mu_cullis_ref,ng3) %>% summarise(mu_mu=mean(mu_cullis2)) 




heritability <- heritability %>% filter(ng %in% c(3,10,30,100))

heritability_mean <- expand.grid(ng=unique(genotypic_se_g0$ng),
                                 sigma_g=unique(genotypic_se_g0$sigma_g),
                                 sigma_e=unique(genotypic_se_g0$sigma_e)) %>%
  mutate(model="Expected H2",
         mu_cullis_ref=sigma_g^2/(sigma_g^2+sigma_e^2/ng),
         mu_mu=sigma_g^2/(sigma_g^2+1.4*sigma_e^2/ng))
#ng=unique(heritability$ng3) ) 


heritability_mean <- heritability_mean %>% filter(ng %in% c(3,10,30,100))

heritability_mean$ng <-  factor(heritability_mean$ng,levels = c(3,10,30,100),
                                labels = c(expression(r == 3), 
                                           expression(r == 10), 
                                           expression(r == 30), 
                                           expression(r == 100)))



ggplot(heritability_mean,aes(mu_cullis_ref,mu_mu)) + geom_point()+facet_grid(.~ng)

heritability$model <-  factor(heritability$model,levels = c("Th", "PL", "MM"))            
heritability_plot <-heritability  %>% filter(ng!=10,model !='PL') %>%
  
  ggplot(aes(x=mu_cullis_ref,y= mu_mu)) +
  #geom_line(data=heritability_mean,aes(x=mu_cullis_ref,y= mu_mu)) +
  geom_point(aes(color=model),alpha=0.3) +
  geom_abline(intercept =  0,slope=1,linetype=2) +
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th rank-based", PL="PL rank-based",MM="MM quant.-based"))+
  xlab("Heritability, one-stage analysis calculation")+
  ylab("Heritability, two-stage analysis calculation")+ 
  facet_wrap(~ng3,labeller = label_parsed) + theme_bw()
#geom_abline(aes(intercept = 0, slope = 1, linetype = "Expected H²"), color = "black") 

ggsave('Results/Paper_figures_V5/cullis_heritability_EUCARPIA.png',heritability_plot,width =8.5 ,height = 2.9)

# Observed dataset

observed_dataset <- read.csv("Results/Plot_inputs/genotypic_means_observed_dataset.csv")

observed_dataset$crop[observed_dataset$crop == "Intermerdia-maturing maize"] <- "Intermedia-maturing maize"


observed_dataset$crop <- factor(observed_dataset$crop,levels = c("Cassava","Groundnut",
                                                                 "Early-maturing maize",
                                                                 
                                                                 "Intermedia-maturing maize"))



rg_labels <- data.frame(
  crop = c("Cassava", "Groundnut", "Early-maturing maize", "Intermedia-maturing maize"),
  rg_value = c(0.99, 0.97, 0.99, 0.97),
  x = c(0),  # Adjust x position (relative)
  y = c(4)   # Adjust y position (relative)
)

obs_dataset_results <- observed_dataset %>% ggplot(aes(x=estimate_blues,y=estimate_yieldrank)) +  
  geom_smooth(method = 'lm',se=F,color='black',linewidth=0.6) + facet_wrap(~crop,scales = 'free')+
  xlab(expression(paste('Estimated genotypic means quant.-based (t  ha'^"-1",")"))) +ylab('Estimated genotypic means ranked-based')+
  geom_point(color='azure4')  + #geom_text_repel(  #+ geom_text(aes(label = genotype))
  #aes(label=genotype),
  #box.padding = unit(0.2, "lines"),size=2,max.overlaps = 25
  
  #)+ 
  stat_cor(
    cor.coef.name = 'tau',
    method = "kendall",
    label.x.npc = 0.001, # Adjust placement to the right
    label.y.npc = 0.99 ,# Adjust placement to the top
    p.accuracy = 0.001, r.accuracy = 0.01
  )+
  geom_text_npc(
    data = rg_labels,
    aes(label = paste("r[g]==", rg_value)),
    parse = TRUE,
    npcx = 0.045,
    npcy = 0.87 
  ) +
  theme_bw()

obs_dataset_results


ggsave('Results/Paper_figures_V5/observed_datasets.png',obs_dataset_results,height = 6,width = 9)








