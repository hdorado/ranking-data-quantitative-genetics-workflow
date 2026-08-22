
# Only paper figures
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

plot_val <- data_30_exmple %>% filter(iter==3 & model!='MM') %>% 
  ggplot(aes(x=scl_true_means,y=est_means)) + 
 # geom_text_repel(  #+ geom_text(aes(label = genotype))
  #  aes(label=genotype,color=model),
  #  box.padding = unit(0.2, "lines"),size=3,max.overlaps = 25
  #)+ 
  #facet_wrap(~iter,scales='free') + 
  geom_abline(intercept = 0,slope = 1)+ 
  geom_smooth(aes(color=model),se=F,method = 'lm')+  
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based", PL="PL Utility-based",MM="MM metric-based"))+
  ylab('Estimated genotypic means') + geom_point(aes(color=model))+
  xlab("True genotypic means") + 
  scale_x_continuous(breaks = c(-1,0,1, 2),limits = c(-1,2.3))+
  scale_y_continuous(breaks = c(-1,0,1, 2, 3)) +  
  #geom_errorbar(aes(ymax = est_means+SE,ymin=est_means-SE,colour = model),
  #width=0.05, alpha=0.7, linewidth=0.5)+
  theme_bw() 


ggsave('Results/Final_figures/tiff/Fig2.tiff',plot_val,height = 3.2,width = 5)

ggsave('Results/Final_figures/pdf/Fig2.pdf',plot_val,height = 3.2,width = 5)

# option 2, shorter

plot_val <- data_30_exmple %>% filter(iter==3 & model!='MM') %>% 
  ggplot(aes(x=scl_true_means,y=est_means)) + 
 # geom_text_repel(  #+ geom_text(aes(label = genotype))
  #  aes(label=genotype,color=model),
  #  box.padding = unit(0.2, "lines"),size=3,max.overlaps = 25
  #)+ 
  #facet_wrap(~iter,scales='free') + 
  geom_abline(intercept = 0,slope = 1)+ 
  geom_smooth(aes(color=model),se=F,method = 'lm')+  
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based", PL="PL Utility-based",MM="MM metric-based"))+
  ylab('Estimated genotypic means') + geom_point(aes(color=model))+
  xlab("True genotypic means") + 
  scale_x_continuous(breaks = c(-1,0,1, 2),limits = c(-1,2.3))+
  scale_y_continuous(breaks = c(-1,0,1, 2, 3)) +  
  #geom_errorbar(aes(ymax = est_means+SE,ymin=est_means-SE,colour = model),
  #width=0.05, alpha=0.7, linewidth=0.5)+
  theme_bw() +  theme(
  legend.position = c(0.03, 0.97),
  legend.justification = c(0, 1),
  legend.background = element_rect(
    fill = "white",
    colour = "black"
  )#,
  #legend.key.size = unit(0.30, "cm"),
  #legend.spacing.y = unit(0.05, "cm"),
  #legend.margin = margin(2, 2, 2, 2),
  #legend.text = element_text(size = 7.5),
  #legend.title = element_text(size = 8.5)
)

ggsave('Results/Final_figures/pdf/Fig2.pdf',plot_val,width = 3.3,height = 3)


# Scaled case

data_exmple_scl <- data_30_exmple %>% filter(iter==3 ,genotype!='G1')

data_exmple_scl <- data_exmple_scl %>% mutate(est_means2= est_means, est_means2= ifelse(model=="PL", 
                                                                                        est_means*0.8372322,
                                                                                        est_means ),SE_2=SE,SE_2=ifelse(model=="PL", 
                                                                                                                        SE_2*0.8372322,
                                                                                                                        SE_2 ))


data_exmple_scl$model <- factor(data_exmple_scl$model,levels = c("Th","PL","MM"))

plot_serror <- data_exmple_scl %>% filter(iter==3 ,genotype!='G1')  %>% 
  ggplot(aes(x=est_means2,y=SE_2,colour = model)) +geom_point() +
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based",PL="PL Utility-based (scl)",MM= "MM Quant.-based"))+ theme_bw()+
  xlab('Estimated genotypic means') + ylab('SE of estimated genotypic means') +ylim(c(0.3,0.5))+
theme(
  legend.position = c(0.03, 0.97),      # top-left inside plot
  legend.justification = c(0, 1),
  legend.background = element_rect(
    fill = "white",
    colour = "black"
  ),
  legend.key.size = unit(0.4, "cm"),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 9)
)


ggsave('Results/Final_figures/pdf/fig3.pdf',plot_serror,height = 3,width = 3.3)

se_dataset <- read.csv("Results/Plot_inputs/summary_analysis_simulation_SE.csv")

se_dataset$genotype <- factor(se_dataset$genotype,levels=mixedsort((unique(se_dataset$genotype))))

#  Adjust scale

scale_pl_th_ajs <- se_dataset %>% select(iter,genotype,model,est_means) %>% 
  pivot_wider(names_from = model,values_from = est_means) %>% 
  group_by() %>% mutate(slope=coefficients(lm(Th~-1+PL))[1])


se_dataset <- se_dataset %>% left_join(scale_pl_th_ajs[c('iter','genotype','slope')],by=c('iter', 'genotype'))

se_dataset <- se_dataset %>% filter(model!="MM"  ) %>% mutate(est_means2 = ifelse(model=="PL",est_means*slope ,est_means),
                                                SE2 = ifelse(model=="PL",SE*slope ,SE))

se_dataset_summ <- se_dataset %>% filter(genotype!='G1' ) %>% filter(model!="MM"  )  %>% filter(iter!=25  ) %>%
  group_by(ng,genotype,sigma_e,sigma_g,model) %>% 
  summarise(mean_SE=mean(SE2),SE_mean=1.96*sd(SE2)/sqrt(n()),sd_estimate=sd(est_means2),
            sd_SE1 = sqrt((n()-1)*sd(est_means2)^2/qchisq(0.025,n()-1,lower.tail = T)),
            sd_SE2 = sqrt((n()-1)*sd(est_means2)^2/qchisq(0.025,n()-1,lower.tail = F)),
            median_SE=median(SE2),sd_me_estimate=mad(est_means2),n=n())



#se_dataset_summ <- mutate(correct_mean = ifelse(model=='PL',mean_SE,),correct_se=ifelse(model=='PL',,))

se_dataset_summ$model <- factor(se_dataset_summ$model,levels = c("Th","PL"))

se_val_plot <- se_dataset_summ %>%  ggplot(aes(x=sd_estimate,y=mean_SE)) + geom_point(aes(color=model)) +
  geom_abline(intercept=0,slope=1)+theme_bw() +  
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based", PL="PL Utility-based (scl)",MM="MM quant.-based")) +

  xlab("SD of estimated genotypic means") + ylab("Mean of SEs of estimated genotypic means")+
  ylim(0.35,0.52)+xlim(0.35,0.5) +   theme(
    legend.position = c(0.04, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = "white",
      colour = "black"
    )#,
    #legend.key.size = unit(0.35, "cm"),
    #legend.spacing.y = unit(0.05, "cm"),
    #legend.margin = margin(2, 2, 2, 2),
    #legend.text = element_text(size = 8),
    #legend.title = element_text(size = 9)
  )

se_val_plot

ggsave('Results/Final_figures/tiff/Fig4.png',se_val_plot,height = 3.5,width =5)
ggsave('Results/Final_figures/pdf/fig4.pdf',se_val_plot,height = 3.5,width =3.3)

# missing values and indeterminate cases analysis

extreme_cases <- read.csv('Results/Plot_inputs/outlayers_inderterminates.csv')

extreme_cases$model <- factor(extreme_cases$model ,levels = c("Th","PL"))

4*14*11*30

extreme_cases$ng2 <- factor(extreme_cases$ng,levels = c(3,10,30,100),
                                 labels = c(expression(n[r] == 3), 
                                           # expression(r[g] == 5), 
                                           expression(n[r] == 10), 
                                           expression(n[r] == 30), 
                                           expression(n[r] == 100)))

custom_colors2 <- c("Th" = "#00CCBE",  # Deep Red
                    "PL" = "#FEA800",
                    "MM"= "#F8766D") 



extreme_cases <- extreme_cases %>% pivot_longer(-c(X:n,ng2),names_to='cases',values_to ='percentaje')

extreme_cases$cases2 <- "Extreme estimated genotypic values"

extreme_cases$cases2[extreme_cases$cases=='perc_sing'] <- "Undefined estimated standard errors"

extreme_cases<- extreme_cases[extreme_cases$cases2 != 'Extreme estimated genotypic values',]

extreme_cases_plot <- extreme_cases %>% ggplot(aes(x=sigma_std,y=percentaje)) +
  geom_point(aes(color=model)) +geom_line(aes(color=model)) + 
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based", PL="PL Utility-based"))+
  facet_grid(.~ng2, labeller = label_parsed)+
  theme_bw()+
  ylab("Percentage of trials") + xlab(expression(paste("True genetic SD ", (sigma[geno]), "")))


ggsave('Results/Final_figures/png/Fig6.png',extreme_cases_plot,height = 3.2,width = 9, dpi=600)


ggsave('Results/Final_figures/pdf/Fig6.pdf',extreme_cases_plot,height = 3.2,width = 9)

#--------------------------


# equivalence

PL_th_equivalence <- read.csv('Results/Plot_inputs/equva_factor_ds_val.csv')

PL_th_equivalence <-PL_th_equivalence[round(PL_th_equivalence$Ratio,8) %in% c(0.18181818,0.71428571 ,0.50000000,2.25000000,4.50000000),] 

PL_th_equivalence$Ratio <- factor(round(PL_th_equivalence$Ratio,2))


PL_th_equivalence2 <- PL_th_equivalence %>% ggplot(aes(x=ng,y=slope_mean))+ 
  geom_point(aes(color=factor(Ratio))) + geom_line(aes(color=factor(Ratio))) +
  theme_bw()+ylab("Regression slope (TH vs. PL estimates)")+
  labs(color= expression(paste("True genetic SD (", sigma[geno], ")")))+
  xlab(expression(paste("Number of replicates (", n[r], ")")))

ggsave('Results/Paper_figures_V5/PL_th_equivalence.png',PL_th_equivalence2,height = 3.5 ,width =6.5 )


# Magical factor

mag_fact2 <- read.csv('Results/Plot_inputs/estimation_genotypic_variance.csv')

mag_fact2 <- mag_fact2 %>%  mutate(ratio=sigma_g/sigma_e,Valid_SE = !Singular & !invalid_SE )%>% filter(Valid_SE) 

mag_fact2b <- mag_fact2[round(mag_fact2$ratio,8) %in% c(0.18181818,0.71428571,0.5,2.25,1.5,3),] 

mag_fact2b <- mag_fact2b %>% select(iter:model,ratio,SE_median,slope_th_pl,Valid_SE)

MM_mag_fact2b <- mag_fact2b %>% filter(model=='MM')

mag_fact2b <- mag_fact2b %>% filter(model  !='MM') %>% 
  left_join(MM_mag_fact2b[c('iter',  'ng' ,'sigma_e', 'sigma_g' ,'SE_median')],by=c('iter',  'ng' ,'sigma_e', 'sigma_g' ), suffix = c('RK','MM')) %>%
  mutate(SE_medianMM=SE_medianMM/sigma_e)  %>% # standarize metric
  mutate(SE_medianRK=ifelse(model=="PL",SE_medianRK*slope_th_pl,SE_medianRK)) # standarize PL

mag_fact2b <- mag_fact2b %>% mutate(ratio2=round(ratio,2),ratio_se = SE_medianRK/SE_medianMM )

mag_fact2b <- mag_fact2b %>% group_by(ratio2,model,ng) %>% 
  summarise(medn=median(ratio_se,na.rm = T)) 


mag_fact2b <- mag_fact2b %>% mutate(model2 = ifelse(model=="PL", "Plackett-luce (scl)","Thurstonian"))

mag_fact2b$model2 <- factor(mag_fact2b$model2,levels = c("Thurstonian","Plackett-luce (scl)"))

ggmagFac <- mag_fact2b %>% filter(ng>=5) %>%   #[round(mag_fact2$ratio,8) %in% c(0.18181818,0.71428571,0.5,2.25,4.50000000),] %>% 
  ggplot(aes(x=ng,y=medn)) +geom_point(aes(color=factor(ratio2))) + 
  labs(color= expression(paste("True genetic SD ", (sigma[geno]), "")))+
  geom_line(aes(color=factor(ratio2)))+theme_bw()+ylim(0,5)+ 
 facet_grid(.~model2) + 
  ylab('Median of ratio of SE estimates\n ranking vs quantitative data')+
  xlab(expression(paste("Number of replicates (", n[r], ")")))

mean(mag_fact2b[mag_fact2b$ng>=30 & mag_fact2b$model== 'Th' & mag_fact2b$ratio2<2,]$medn)

ggsave("Results/Paper_figures_V5/magical_factor.png",ggmagFac,height = 4,width = 8)


#mag_fact2 <- mag_fact2 %>% filter(ng!=5) 


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
                                 labels = c(expression(n[r] == 3), #r[g] 
                                            expression(n[r] == 5), 
                                            expression(n[r] == 10), 
                                            expression(n[r] == 30), 
                                            expression(n[r] == 100)))

#-------------------------------------------------------------------------------
# GV - Person correlation + kendall-index vs radio of variation 

performance_Values$model <- factor(performance_Values$model,levels = c("Th","PL","MM"))

pv <- performance_Values %>% 
  filter(ng %in% c(3,10,30,100) & 
           (variable %in% c('Bias','RMSE',
                            "paste(\"Kendall correlation coef. (\", tau, \")\")")
             )) %>%
  ggplot(aes(x=Ratio,y=mu_value,group=interaction(model,ng))) +
  geom_point(aes(colour = model),alpha=0.6,size=1.3) +
  geom_vline(xintercept = 0.8888889,linetype= 'dashed')+
  geom_errorbar(aes(ymax = mu_value+SE,ymin=mu_value-SE,colour = model),
                width=0.05, alpha=0.8, linewidth=0.5)+
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based", PL="PL Utility-based (scl)",MM="MM Quant.-based"))+
  facet_grid(variable~ng2,labeller =label_parsed,scales='free')+theme_bw() + #,switch="both"
  ylab('') + geom_line(aes(colour = model),alpha=0.7)+
  xlab(expression(paste("True genetic SD ", (sigma[geno]), "")))


pv

ggsave('Results/Final_figures/png/Fig5.png',pv,height = 6,width = 9,dpi = 600)

ggsave('Results/Final_figures/pdf/fig5.pdf',pv,height = 6,width = 9)

#-------------------------------------------------------------------------------

# Estimation genotypic variances new version


genotypic_se_g0 <- read.csv('Results/Plot_inputs/estimation_genotypic_variance.csv')  # UPDATE THIS FILE FROM THE SOURCE

head(genotypic_se_g0)

genotypic_se_g0 <- genotypic_se_g0 %>% mutate(experiment = paste(iter, ng ,sigma_e, sigma_g,sep='-'))

genotypic_se_g0 <- genotypic_se_g0 %>% filter(ng!=5)

#----------------------- Counting missing cases--------------------------------

# analysis of missing estimates of genetic variance an heritability

genotypic_se_g0$sigma_g 

11*14*4*30

length(unique(genotypic_se_g0$experiment)) # full experiment

#

table(genotypic_se_g0$invalid_SE,genotypic_se_g0$invalid_qSE,genotypic_se_g0$model)
      
genotypic_se_g0 %>% group_by(model) %>%
  summarise(n=n(),invalid_SE_ag=sum(invalid_SE),invalid_qSE_ag=sum( !invalid_SE & invalid_qSE),n_estm = sum(is.na(sd_g_w.qV)),
                                                  n_hertab=sum(is.na(qh2.2st.w.qV)|qh2.2st.w.qV<0))

# Counting singular cases

table(genotypic_se_g0$Singular,genotypic_se_g0$invalid_SE,genotypic_se_g0$model)

genotypic_se_g0 <- genotypic_se_g0 %>% mutate(Valid_qSE = !Singular & !invalid_SE & !invalid_qSE)

genotypic_se_g0 <- genotypic_se_g0 %>% mutate(sd_g_w_.qv_filtered = ifelse(Valid_qSE,sd_g_w.qV,NA))

table(genotypic_se_g0$Valid_qSE,!is.na(genotypic_se_g0$sd_g_w.qV),genotypic_se_g0$model)

genotypic_se_g0[!genotypic_se_g0$Valid_qSE ,]

genotypic_se_g0 %>% group_by(model) %>%
  summarise(n=n(),Singular_ag=sum(Singular),invalid_SE_ag=sum(!Singular & invalid_SE),
            invalid_qSE_ag=sum( !Singular & !invalid_SE & invalid_qSE),
            nton__estm = sum(is.na(sd_g_w_.qv_filtered)),
            n_hertab=sum(is.na(qh2.2st.w.qV)|qh2.2st.w.qV<0))


# paste valid cases

# full_ds_valid_cases <- read.csv('Results/full_experiment_valid_novalidcases.csv')

# genotypic_se_g0 <- genotypic_se_g0 %>% mutate( Experiment0 = paste(iter,sigma_g,sigma_e,ng,sep='-') )

# genotypic_se_g0 <-  genotypic_se_g0 %>% left_join(unique(full_ds_valid_cases[c('Experiment0','Exp_valid')]))

# write.csv(genotypic_se_g0,'Results/Plot_inputs/estimation_genotypic_variance.csv',row.names=F)

#genotypic_se_g0 <- read.csv('Results/Plot_inputs/estimation_genotypic_variance.csv')  # UPDATE THIS FILE FROM THE SOURCE

genotypic_se_g0 <- genotypic_se_g0 %>% filter(ng %in% c(3,10,30,100)  ) 

genotypic_se_g0 <- genotypic_se_g0 %>% mutate(experiment = paste(iter, ng ,sigma_e, sigma_g,sep='-'))

length(unique(genotypic_se_g0$experiment))

custom_colors2 <- c("Th" = "#00CCBE",  # Deep Red
                    "PL" = "#FEA800",
                    "MM"= "#F8766D") 

# Sample label

genotypic_se_g0$ng3 <- factor(genotypic_se_g0$ng,levels = c(3,10,30,100),
                             labels = c(expression(n[r] == 3), 
                                        expression(n[r] == 10), 
                                        expression(n[r] == 30), 
                                        expression(n[r] == 100)))


length(unique(genotypic_se_g0$experiment))

# complete plot before aggregation


# all estimates of genetic varaince

bef_agg <- genotypic_var_estm_2_appx <- genotypic_se_g0 %>% filter(Valid_qSE) %>% mutate(true_sigma_sc=sigma_g/sigma_e) %>%  ggplot(aes(y=sg_g_estimate_scl2 ,x=true_sigma_sc,color=model ))+ 
  geom_point(alpha=0.3,show.legend = T)+
  # geom_errorbar(aes(ymin=mu_sd_g_wQV_scaled-1.96*SE_mu_sd_g_wQV_scaled,ymax=mu_sd_g_wQV_scaled+1.96*SE_mu_sd_g_wQV_scaled,colour = model),width=0.1, alpha=0.7, linewidth=0.3)+
  #geom_line(aes(color=model),alpha=0.8,linewidth=0.5)+
  geom_abline(intercept =  0,slope=1,linetype=1)+
  #scale_shape_manual(name='Number \nof replicates',values = c(19,17,15,3,8,7))+
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based", PL="PL Utility-based (scl)",MM="MM Quant.-based"))+
  scale_linetype(name=expression(atop("Residual ",paste('stad. dev. ',sigma[e]))))+
  xlab(expression(paste("True genetic SD ", (sigma[g]), "")))+
  geom_vline(xintercept = 0.8,linetype =2)+
  #facet_wrap(~ng)+
  ylab(expression(paste("Estimated genetic SD ", (hat(sigma[g])), "")))+#+ggtitle(label = 'include out')+
  facet_wrap(~ng3,labeller =label_parsed)+theme_bw()#+ylim(0,5.5)

bef_agg

ggsave("Results/Paper_figures_V5/full_genetic_variance_appx.png",bef_agg,width =7 ,height = 6)


####


genotypic_se_g <- genotypic_se_g0 %>% filter(Valid_qSE) %>% mutate(true_sigma_sc=sigma_g/sigma_e) %>% 
  group_by(ng, true_sigma_sc, model) %>%# filter(!is.na(sg_g_estimate_scl2)) %>% 
  summarise(mu_sd_g_wQV = mean(sg_g_estimate_scl2,na.rm=T),
            SE_mu_sd_g_wQV = sd(sg_g_estimate_scl2,na.rm=T)/sqrt(n()) ,n=n())


genotypic_se_g <- genotypic_se_g %>% filter(ng %in% c(3,10,30,100)  ) 


genotypic_se_g$ng3 <- factor(genotypic_se_g$ng,levels = c(3,10,30,100),
                              labels = c(expression(n[r] == 3), 
                                         expression(n[r] == 10), 
                                         expression(n[r] == 30), 
                                         expression(n[r] == 100)))


genotypic_se_g$model <- factor(genotypic_se_g$model,levels = c("Th","PL",'MM'))

# Aggreated

gene_vars <- genotypic_se_g%>% ggplot(aes(y=mu_sd_g_wQV ,x=true_sigma_sc,color=model ))+ 
  geom_point(alpha=0.3,show.legend = T)+
 # geom_errorbar(aes(ymin=mu_sd_g_wQV_scaled-1.96*SE_mu_sd_g_wQV_scaled,ymax=mu_sd_g_wQV_scaled+1.96*SE_mu_sd_g_wQV_scaled,colour = model),width=0.1, alpha=0.7, linewidth=0.3)+
  geom_line(aes(color=model),alpha=0.8,linewidth=0.5)+
  geom_abline(intercept =  0,slope=1,linetype=1)+
  #scale_shape_manual(name='Number \nof replicates',values = c(19,17,15,3,8,7))+
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based", PL="PL Utility-based (scl)",MM="MM Quant.-based"))+
  scale_linetype(name=expression(atop("Residual ",paste('stad. dev. ',sigma[e]))))+
  xlab(expression(paste("True genetic SD ", (sigma[g]), "")))+
  geom_vline(xintercept = 0.8,linetype =2)+
  #facet_wrap(~ng)+
  ylab(expression(paste("Estimated genetic SD ", (hat(sigma[g])), "")))+#+ggtitle(label = 'include out')+
  facet_wrap(ng3~.,labeller =label_parsed)+theme_bw()+ylim(0,5)

ggsave("Results/Paper_figures_V5/genotypic_var_estm_2.png",gene_vars,width =7 ,height = 6)

ggsave("Results/Final_figures/png/Fig7.png",gene_vars,width =6.8 ,height = 4,dpi = 600)

ggsave("Results/Final_figures/pdf/fig7.pdf",gene_vars,width =7 ,height = 6)

# Contraste#-----------------------------------------------------------------

?ggsave

genotypic_se_g <- genotypic_se_g0 %>% filter(Valid_qSE) %>% mutate(true_sigma_sc=sigma_g/sigma_e) %>% 
  group_by(ng, true_sigma_sc, model) %>%# filter(!is.na(sg_g_estimate_scl2)) %>% 
  summarise(mu_sd_g_wQV = mean(sg_g_estimate_scl2,na.rm=T),
            mu_sdSE_g_w= mean(sg_gse_estimate_scl2,na.rm=T) ,n=n())

genotypic_se_g <- genotypic_se_g %>% filter(ng %in% c(3,10,30,100)  ) 


genotypic_se_g$ng3 <- factor(genotypic_se_g$ng,levels = c(3,10,30,100),
                             labels = c(expression(n[r] == 3), 
                                        expression(n[r] == 10), 
                                        expression(n[r] == 30), 
                                        expression(n[r] == 100)))


genotypic_se_g$model <- factor(genotypic_se_g$model,levels = c("Th","PL",'MM'))

genotypic_se_g <- genotypic_se_g %>% pivot_longer(!c(ng,true_sigma_sc, model,n ,ng3 ),names_to='method',values_to='values')

library(stringr)


genotypic_se_g <- genotypic_se_g %>% mutate(method=str_replace_all(method,c('mu_sdSE_g_w'='Variances-based','mu_sd_g_wQV'='Quasi-variances-based')))

genotypic_se_g$method <- factor(genotypic_se_g$method , levels = c('Variances-based','Quasi-variances-based'))

gene_vars <- genotypic_se_g%>% filter( model!='MM',ng==30) %>% 
  ggplot(aes(y=values ,x=true_sigma_sc,color=model ))+ 
  geom_point(alpha=0.3,show.legend = T)+
  # geom_errorbar(aes(ymin=mu_sd_g_wQV_scaled-1.96*SE_mu_sd_g_wQV_scaled,ymax=mu_sd_g_wQV_scaled+1.96*SE_mu_sd_g_wQV_scaled,colour = model),width=0.1, alpha=0.7, linewidth=0.3)+
  geom_line(aes(color=model),alpha=0.8,linewidth=0.5)+
  geom_abline(intercept =  0,slope=1,linetype=1)+
  #scale_shape_manual(name='Number \nof replicates',values = c(19,17,15,3,8,7))+
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based", PL="PL Utility-based (scl)",MM="MM Quant.-based"))+
  scale_linetype(name=expression(atop("Residual ",paste('stad. dev. ',sigma[e]))))+
  xlab(expression(paste("True genetic SD ", (sigma[g]), "")))+xlim(0,2)+ylim(0,2)+
  geom_vline(xintercept = 0.8,linetype =2)+
  #facet_wrap(~ng)+
  ylab(expression(paste("Estimated genetic SD ", (hat(sigma[g])), "")))+#+ggtitle(label = 'include out')+
  facet_wrap(method~.,labeller =label_parsed)+theme_bw()#+ylim(0,5.5)

ggsave("Results/Paper_figures_V5/qv-vrs.png",gene_vars,width =7 ,height = 3)


#------------ Estimation of heritabilities new version-------------------------


genotypic_se_g0 <- genotypic_se_g0 %>% mutate(valid_sdg = Valid_qSE & !is.na(sd_g_w.qV))

genotypic_se_g0 <- genotypic_se_g0 %>% mutate(valid_sdg = Valid_qSE & !is.na(sd_g_w.qV))

genotypic_se_g0 %>% group_by(model) %>% summarise(n=n())

heritabilities <-  genotypic_se_g0 %>% filter(valid_sdg)

heritabilities <- heritabilities %>% select(iter:model,ng3,valid_sdg,h2_metric_cullis,h2_broad_metric1,h2_broad_metric2,heritability_cullisEstQV,h2_SEd2)

heritabilities <- heritabilities %>% mutate(heritability_valid = h2_metric_cullis>0 & heritability_cullisEstQV>0)

heritabilities %>% group_by(model) %>% summarise(n=n(),h2_n=sum(heritability_valid))


# genotypic_se_g0 <- genotypic_se_g0 %>% select(iter:model,Valid_qSE)

#new_calc_heritab <- read.csv('Results/Paper_figures_V5/new_calc_heritab.csv')

new_calc_heritab <- heritabilities 


custom_colors2 <- c("Th" = "#00CCBE",  # Deep Red
                    "PL" = "#FEA800",
                    "MM"= "#F8766D") 

# Sample label

new_calc_heritab <- new_calc_heritab %>% filter(ng %in% c(3,10,30,100)  ) 


new_calc_heritab$ng3 <- factor(new_calc_heritab$ng,levels = c(3,10,30,100),
                              labels = c(expression(n[r] == 3), 
                                         expression(n[r] == 10), 
                                         expression(n[r] == 30), 
                                         expression(n[r] == 100)))


#new_calc_heritab <- new_calc_heritab %>% select(ng:ng3) %>% unique()

new_calc_heritab %>% group_by(model) %>% summarise(n=n())

new_calc_heritab <- new_calc_heritab %>% filter(heritability_valid)

#  cullis heritability

# Mean of Cullis heritability estimates obtained from the two-stage analysis

plot_heritab_cullis <-  new_calc_heritab %>% 
  ggplot(aes(x=h2_metric_cullis,y=heritability_cullisEstQV,color=model)) +
   ylab('H2 Cullis estim. 2 stg') + xlab('H2 Cullis est. 1 stg') + geom_point(alpha=0.3,show.legend = T)+
geom_abline(intercept =  0,slope=1,linetype=1) +
facet_wrap(ng3~.,labeller =label_parsed)+
scale_color_manual(name= "Estimates",values=custom_colors2,
                   labels = c(Th="Th Utility-based", PL="PL Utility-based (scl)",MM="MM Quant.-based")) +
 # xlim(-1,1)+ylim(-2,1) +
  theme_bw() + xlab("Cullis heritability estimates (one-stage analysis)") + 
  ylab("Cullis heritability estimates (two-stage analysis)")

ggsave('Results/Paper_figures_V5/cullis_hertab.png',plot_heritab_cullis,width =8 ,height = 6)

# SE2 no aggreated

plot_heritab <-  new_calc_heritab %>% 
  ggplot(aes(x=h2_broad_metric2,y=h2_SEd2,color=model)) +
  ylab(expression("Entry-mean basis heritability estimate " ~ (hat(sigma)[g]^2 / (hat(sigma[g])^2 + var(hat(mu)[i])/ 2) ) )) +
  xlab(expression("Entry-mean basis heritability estimate " ~(hat(sigma)[g]^2 / (hat(sigma)[g]^2 + hat(sigma)[e]^2 / n[r])))) + geom_point(alpha=0.3,show.legend = T)+
  geom_abline(intercept =  0,slope=1,linetype=1) +
  facet_wrap(ng3~.,labeller =label_parsed)+
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based", PL="PL Utility-based (scl)",MM="MM Quant.-based")) +
  # xlim(-1,1)+ylim(-2,1) +
  theme_bw() #+# xlab("Cullis heritability estimates (one-stage analysis)") + 
  #ylab("Cullis heritability estimates (two-stage analysis)")
plot_heritab

ggsave('Results/Paper_figures_V5/h2_se2_hertab.png',plot_heritab,width =8 ,height = 6)

#  cullis heritability aggregated

plot_heritab_cullis <-  new_calc_heritab %>% 
  ggplot(aes(x=h2_metric_cullis,y=heritability_cullisEstQV,color=model)) +
  ylab('H2 Cullis estim. 2 stg') + xlab('H2 Cullis est. 1 stg') + geom_point(alpha=0.3,show.legend = T)+
  geom_abline(intercept =  0,slope=1,linetype=1) +
  facet_wrap(ng3~.,labeller =label_parsed)+
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based", PL="PL Utility-based (scl)",MM="MM Quant.-based")) +
  # xlim(-1,1)+ylim(-2,1) +
  theme_bw() + xlab("Cullis heritability estimates (one-stage analysis)") + 
  ylab("Cullis heritability estimates (two-stage analysis)")

ggsave('Results/Paper_figures_V5/cullis_hertab.png',plot_heritab_cullis,width =8 ,height = 6)

heritability_mean <- expand.grid(ng=unique(genotypic_se_g0$ng),
                                 sigma_g=unique(genotypic_se_g0$sigma_g),
                                 sigma_e=unique(genotypic_se_g0$sigma_e)) %>%
  mutate(model="Expected H2",
         mu_cullis_ref=sigma_g^2/(sigma_g^2+sigma_e^2/ng),
         mu_mu=sigma_g^2/(sigma_g^2+1.3*sigma_e^2/ng))
#ng=unique(heritability$ng3) )

heritability_mean$ng3 <- factor(heritability_mean$ng,levels = c(3,10,30,100),
                                labels = c(expression(n[r] == 3), 
                                           expression(n[r] == 10), 
                                           expression(n[r] == 30), 
                                           expression(n[r] == 100)))


plot_heritab_cullis_agg <-  new_calc_heritab %>% select(-iter) %>% 
  group_by( ng, sigma_e, sigma_g, model, ng3) %>%
  summarise(mean_h2=mean(heritability_cullisEstQV),mean_ref_h2=mean(h2_metric_cullis)) %>%
  ggplot(aes(x=mean_ref_h2,y=mean_h2,color=model)) +
  geom_line(data=heritability_mean,aes(x=mu_cullis_ref,y= mu_mu))+
  ylab('H2 Cullis estim. 2 stg') + xlab('H2 Cullis est. 1 stg') + geom_point(alpha=0.3,show.legend = T)+
  geom_abline(intercept =  0,slope=1,linetype=1) +
  facet_wrap(ng3~.,labeller =label_parsed)+
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based", PL="PL Utility-based (scl)",MM="MM Quant.-based")) +
  # xlim(-1,1)+ylim(-2,1) +
  theme_bw() + xlab("Cullis heritability estimates (one-stage analysis)") + 
  ylab("Cullis heritability estimates (two-stage analysis)")

plot_heritab_cullis_agg

ggsave('Results/Paper_figures_V5/em_hertab.png',plot_heritab_cullis_agg,width =8 ,height = 6)
ggsave('Results/Final_figures/pdf/fig8.pdf',plot_heritab_cullis_agg,width =8 ,height = 6)
ggsave('Results/Final_figures/png/Fig8.png',plot_heritab_cullis_agg,width =6.8 ,height = 4,dpi=600)


plot_heritab_cullis_agg <-  new_calc_heritab %>% select(-iter) %>% 
  group_by( ng, sigma_e, sigma_g, model, ng3) %>%
  summarise(mean_h2=mean(h2_SEd2),mean_ref_h2=mean(h2_broad_metric2)) %>%
  ggplot(aes(x=mean_ref_h2,y=mean_h2,color=model)) +
  geom_line(data=heritability_mean,aes(x=mu_cullis_ref,y= mu_mu))+
  ylab('H2 Cullis estim.') + xlab('H2 Cullis estm.') + geom_point(alpha=0.3,show.legend = T)+
  geom_abline(intercept =  0,slope=1,linetype=1) +
  facet_wrap(ng3~.,labeller =label_parsed)+
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based", PL="PL Utility-based (scl)",MM="MM Quant.-based")) +
  # xlim(-1,1)+ylim(-2,1) +
  theme_bw() +ylab('H2 entry mean se^2/2 estim.') + xlab('H2 entry mean est.') 

plot_heritab_cullis_agg

ggsave('Results/Paper_figures_V5/h2_se2_b_agg.png',plot_heritab_cullis_agg,width =8 ,height = 6)

# broad sense heritability

plot_he_em <- new_calc_heritab  %>% filter(mod_h2 != 'MM-qh2.2st.w.qV') %>% ggplot(aes(x = h2_broad_metric2,y = estim_broad_sense_h,color=model)) + 
  ylab('H2 Cullis estim. 2 stg') + xlab('H2 Cullis est. 1 stg') + geom_point(alpha=0.3,show.legend = T)+
  geom_abline(intercept =  0,slope=1,linetype=1) +
  facet_wrap(ng3~.,labeller =label_parsed)+
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based", PL="PL Utility-based (scl)",MM="MM Quant.-based")) +
  #xlim(-1,1)+ylim(-2,1) +
  theme_bw() + xlab("Entry-mean heritability estimates (one-stage analysis)") + 
  ylab("Entry-mean heritability estimates (two-stage analysis)")

ggsave('Results/Paper_figures_V5/em_hertab.png',plot_he_em,width =8 ,height = 6)

# broad sense heritability aggreated

plot_he_mu_em <-  new_calc_heritab %>% filter(mod_h2!='MM-qh2.2st.w.qV') %>% group_by(ng3,mod_h2, sigma_e, sigma_g ,model) %>% 
  summarise(mu_h2=mean(h2_broad_metric2,na.mr=T),mu_es_h2 =mean(estim_broad_sense_h,na.mr=T) ) %>% ggplot(aes(x = mu_h2,y = mu_es_h2,color=model)) + 
  ylab('H2 Cullis estim. 2 stg') + xlab('H2 Cullis est. 1 stg') + geom_point(alpha=0.3,show.legend = T)+
  geom_abline(intercept =  0,slope=1,linetype=1) +
  facet_wrap(ng3~.,labeller =label_parsed)+
  scale_color_manual(name= "Estimates",values=custom_colors2,
                     labels = c(Th="Th Utility-based", PL="PL Utility-based (scl)",MM="MM Quant.-based")) +
  #xlim(-1,1)+ylim(-2,1) +
  theme_bw() + xlab("Mean of Entry-mean heritability estimates (one-stage analysis)") + 
  ylab("Mean of Entry-mean heritability estimates (two-stage analysis)")

ggsave('Results/Paper_figures_V5/mean_em_hertab.png',plot_he_mu_em,width =8 ,height = 6)


#---------------------all code of heritability estimation----------------------

# genotypic_se_g0 <- read.csv('Results/Plot_inputs/estimation_genotypic_variance.csv')  # UPDATE THIS FILE FROM THE SOURCE
# 
# genotypic_se_g0 <- genotypic_se_g0 %>% filter(Exp_valid=="Valid")
# 
# genotypic_se_g0$heritability_fnl[!is.na(genotypic_se_g0$heritability_fnl) & (genotypic_se_g0$heritability_fnl<=0)] <- NA
# 
# genotypic_se_g0$h2_metric_cullis[!is.na(genotypic_se_g0$h2_metric_cullis) & (genotypic_se_g0$h2_metric_cullis<=0)] <- NA
# 
# 
# genotypic_se_g0$heritability_cullisEstQV[!is.na(genotypic_se_g0$heritability_cullisEstQV) & (genotypic_se_g0$heritability_cullisEstQV<=0)] <- NA
# 
# genotypic_se_g0 <-  genotypic_se_g0 %>% filter(Exp_valid=="Valid")
# 




# heritability <- genotypic_se_g0 %>% group_by(ng, sigma_e, sigma_g ,model) %>% 
#   summarise(mu_cullis_ref=mean(h2_metric_cullis,na.rm = T),
#             mu_cullis=mean(heritability_fnl,na.rm=T),
#             mu_cullis2=mean(heritability_cullisEstQV,na.rm = T))
# 
# 
# 
# heritability %>% ggplot(aes(x=mu_cullis_ref ,y=mu_cullis2) ) + 
#   facet_wrap(~ng) + geom_point(aes(color=model))+
#   geom_abline(intercept = 0,slope=1)
# 
# heritability$ng3 <-  factor(heritability$ng,levels = c(3,10,30,100),
#                      labels = c(expression(n[r] == 3), 
#                                 expression(n[r] == 10), 
#                                 expression(n[r] == 30), 
#                                 expression(n[r] == 100)))
# 
# heritability <- heritability %>% ungroup()%>% select(ng,ng3,model,mu_cullis_ref,mu_cullis2) %>% 
#   group_by(ng,model,mu_cullis_ref,ng3) %>% summarise(mu_mu=mean(mu_cullis2)) 
# 
# 
# 
# 
# heritability <- heritability %>% filter(ng %in% c(3,10,30,100))
# 
# heritability_mean <- expand.grid(ng=unique(genotypic_se_g0$ng),
#                                            sigma_g=unique(genotypic_se_g0$sigma_g),
#                                  sigma_e=unique(genotypic_se_g0$sigma_e)) %>%
#                                  mutate(model="Expected H2",
#                                         mu_cullis_ref=sigma_g^2/(sigma_g^2+sigma_e^2/ng),
#                                         mu_mu=sigma_g^2/(sigma_g^2+1.4*sigma_e^2/ng))
#                                  #ng=unique(heritability$ng3) ) 
# 
# 
# heritability_mean <- heritability_mean %>% filter(ng %in% c(3,10,30,100))
# 
# heritability_mean$ng3 <-  factor(heritability_mean$ng,levels = c(3,10,30,100),
#                                 labels = c(expression(n[r] == 3), 
#                                            expression(n[r] == 10), 
#                                            expression(n[r] == 30), 
#                                            expression(n[r] == 100)))
# 
# 
# 
# ggplot(heritability_mean,aes(mu_cullis_ref,mu_mu)) + geom_point()+facet_grid(.~ng)
# 
# heritability$model <-  factor(heritability$model,levels = c("Th", "PL", "MM"))            
# heritability_plot <-heritability  %>%
#   
#   ggplot(aes(x=mu_cullis_ref,y= mu_mu)) +
#   geom_line(data=heritability_mean,aes(x=mu_cullis_ref,y= mu_mu)) +
#   geom_point(aes(color=model),alpha=0.3) +
#   geom_abline(intercept =  0,slope=1,linetype=2) +
#   scale_color_manual(name= "Estimates",values=custom_colors2,
#                      labels = c(Th="Th Utility-based", PL="PL Utility-based",MM="MM Quant.-based"))+
#   xlab("Heritability, one-stage analysis calculation")+
#   ylab("Heritability, two-stage analysis calculation")+ 
#   facet_wrap(~ng3,labeller = label_parsed) + theme_bw()
# #geom_abline(aes(intercept = 0, slope = 1, linetype = "Expected H²"), color = "black") 
# 
# ggsave('Results/Paper_figures_V5/cullis_heritability.png',heritability_plot,width =8 ,height = 6)

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

observed_dataset$crop2 <- str_replace_all(observed_dataset$crop,c('Cassava'='Nigeria - Cassava','Groundnut'='Tanzania - Groundnut','Early-maturing maize'='Kenia - Early-maturing maize','Intermedia-maturing maize'= 'Kenia - Intermedia-maturing maize')) 

observed_dataset


rg_labels <- data.frame(
  crop2 = c("Nigeria - Cassava", "Tanzania - Groundnut", "Kenia - Early-maturing maize", "Kenia - Intermedia-maturing maize"),
  rg_value = c(0.99, 0.97, 0.99, 0.97),
  x = c(0),  # Adjust x position (relative)
  y = c(4)   # Adjust y position (relative)
)


obs_dataset_results <- observed_dataset %>% ggplot(aes(x=estimate_blues,y=estimate_yieldrank)) +  
  geom_smooth(method = 'lm',se=F,color='black',linewidth=0.6) + facet_wrap(~crop2,scales = 'free')+
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





ggsave('Results/Paper_figures_V5/observed_datasets_EUCARPIA.png',obs_dataset_results,height = 6,width = 9)








