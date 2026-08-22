# Script updated to remove only singular cases of cov matrix
# Hugo Dorado
# 07-01-2025
# 175 outlayers***

library(tidyverse)
library(dplyr)
library(Hmisc)
library(ggpmisc)
library(asreml)

# setwd('/home/joost/Hugo_Projects/Sim_est_genotypic_values/')

setwd("Paper_1_Estimation_of_genetic_parameters/")

dataset_summary0 <- read.csv('Results/summary_analysis_iter_g_var_table_V9.csv') # 

#dataset_summary <- dataset_summary %>% select()
dataset_summary %>% group_by(model) %>% summarise(n=n(),n_estm = sum(is.na(sd_g_w.qV)))

# Add reference variety and incoporate reference variables

dataset_summary <- dataset_summary0 %>% 
  mutate(Experiment0= paste(iter, sigma_g, sigma_e, ng,sep='-'),
         Experiment = paste(iter, ng, sigma_g,sigma_e,model,sep='-'))%>% 
  group_by(iter, ng, sigma_e, sigma_g, model ) %>% 
  mutate(ref_val=value.genotype[genotype=='G1'],ref_est_val=est_means[genotype=='G1'])

# Estandarize the metric values

dataset_summary <- dataset_summary %>% 
  mutate(scl_true_means=(value.genotype-ref_val)/sigma_e,
         scl_est_means=ifelse(model=="MM",
                              (est_means-ref_est_val)/sigma_e,est_means))

# Equivalence validation

equva_factor_ds0 <- dataset_summary %>% mutate(Ratio = sigma_g/sigma_e) %>%
  select(iter,ng, sigma_e, sigma_g,value.genotype,scl_true_means,model,Ratio,est_means,genotype,quasiSE) %>%
  pivot_wider(values_from =  c(est_means,quasiSE) ,names_from = c(model))


equva_factor_ds <- equva_factor_ds0 %>% group_by(iter,ng, sigma_e, sigma_g) %>% 
  mutate(outalayer_th=ifelse(sum(abs(est_means_Th)>4)==0,"Not",'Yes'))

equva_factor_ds <- equva_factor_ds %>% group_by(iter,ng, sigma_e, sigma_g) %>% 
  mutate(slope_th_pl=sum(est_means_Th*est_means_PL)/sum(est_means_PL^2))

equva_factor_ds <- equva_factor_ds %>% select(iter,ng,sigma_e,sigma_g,Ratio,outalayer_th,slope_th_pl) %>% unique()

equva_factor_ds_full <- equva_factor_ds %>% group_by(ng,  Ratio) %>% summarise(slope_mean= mean(slope_th_pl) ) 

equva_factor_ds_val <- equva_factor_ds %>% filter(outalayer_th=="Not")  %>% group_by(ng,  Ratio) %>% summarise(slope_mean= mean(slope_th_pl) )

equva_factor_ds_full[round(equva_factor_ds_full$Ratio,8) %in% c(0.18181818,0.71428571,0.5,2.25,4.50000000),] %>% 
  ggplot(aes(x=ng,y=slope_mean)) +geom_point(aes(color=factor(Ratio))) + geom_line(aes(color=factor(Ratio)))+theme_bw()

write.csv(equva_factor_ds_full,'Results/Plot_inputs/equva_factor_ds_val.csv',row.names = F)




equva_factor_ds_val[round(equva_factor_ds_val$Ratio,8) %in% c(0.18181818,0.71428571,0.5,2.25,4.50000000),] %>% 
  ggplot(aes(x=ng,y=slope_mean)) +geom_point(aes(color=factor(Ratio))) + geom_line(aes(color=factor(Ratio)))+theme_bw()

# Independent analysis of correlation

equva_factor_ds0 <- equva_factor_ds0 %>% filter(genotype!="G1") %>% mutate(Ratio=sigma_g/sigma_e) %>% group_by(iter,ng, sigma_e ,sigma_g,  Ratio) %>% 
  mutate(person_Th=cor(value.genotype,est_means_Th), person_Pl=cor(value.genotype,est_means_PL),
         RMSE_Th=sqrt(mean((est_means_Th-value.genotype))^2),RMSE_PL=sqrt(mean((est_means_PL-value.genotype))^2),n=n())

one_iter <- equva_factor_ds0 %>% dplyr::filter(iter==15,ng==10,Ratio==3.75) 


one_iter[abs(one_iter$est_means_Th)<4 , ]%>% 
  ggplot(aes(x=value.genotype)) +geom_point(aes(y=est_means_Th,color='Thurstonian'))+
  geom_point(aes(y=est_means_PL,color='PL')) + 
  geom_smooth(aes(y=est_means_Th,color='Thurstonian'),method = 'lm',se = F)#+
#  geom_smooth(aes(y=est_means_PL,color='PL'),method = 'lm',se = F) + theme_bw()#+ 
# geom_text(x = -Inf, y = Inf, aes(label = paste("RMSE: PL=",round(((RMSE_PL)),3),
#                                               "Th=",round(((RMSE_Th)),3),
#                                              "\nCOR: ",'PL=',round(((person_Pl)),3),' Th=',round(((person_Th)),3) ,sep=' ')),
#        hjust = 0, vjust = 1, size = 5)

one_iter%>% 
  ggplot(aes(x=value.genotype)) +geom_point(aes(y=est_means_Th,color='Thurstonian'))

one_iter2<-one_iter[abs(one_iter$est_means_Th)<4 , ]

cor(one_iter2$est_means_Th,one_iter2$value.genotype)


one_iter%>% 
  ggplot(aes(x=value.genotype)) +geom_point(aes(y=est_means_PL,color='PL'))

sort(unique(equva_factor_ds0$Ratio))

# 

dataset_RMSE_fixed <- equva_factor_ds0

dataset_RMSE_fixed <-  dataset_RMSE_fixed %>% group_by(iter,  ng, sigma_e, sigma_g, Ratio) %>% 
  mutate(slope_th_pl=sum(est_means_Th*est_means_PL)/sum(est_means_PL^2)) 

dataset_RMSE_fixed <- dataset_RMSE_fixed %>% mutate(est_means_PL_scl=est_means_PL*slope_th_pl)

data.frame(head(dataset_RMSE_fixed))

dataset_RMSE_fixed_long <- dataset_RMSE_fixed %>% select(iter, ng, sigma_e, sigma_g ,value.genotype,scl_true_means, Ratio, genotype,est_means_Th, est_means_PL_scl, est_means_MM) %>%
  pivot_longer(-(iter:genotype),names_to = "Model" ,names_prefix ="est_means_",
               values_to = "est_mean")

dataset_RMSE_fixed_long <- dataset_RMSE_fixed_long %>% filter(Model!= "MM") %>% group_by(iter,ng, sigma_e, sigma_g,Model,Ratio ) %>% summarise(RMSE=mean((est_mean-scl_true_means)^2))

dataset_RMSE_fixed_long <- dataset_RMSE_fixed_long %>% group_by(ng,sigma_e,sigma_g,Model,Ratio)  %>% summarise(meanRMSE=mean(RMSE))

dataset_RMSE_fixed_long %>% ggplot(aes(x=Ratio,y=meanRMSE)) + geom_point(aes(color=Model)) +
  geom_line(aes(color=Model)) + facet_wrap(~ng,scales='free')


# Transformed

# Unbiased estimates

dataset_summary_out_30_ref_case <- dataset_summary %>% filter(sigma_e==0.9,sigma_g==0.8,ng==30)

dataset_summary_out_30_ref_case %>% ggplot(aes(x=scl_true_means,y=est_means)) + geom_point(aes(color=model)) + 
  facet_wrap(~iter,scales='free')+ geom_abline(intercept=0,slope=1)+  geom_smooth(aes(color=model),se=F,method = 'lm')

dataset_summary_out_30_ref_case %>% filter(iter==3) %>% ggplot(aes(x=scl_true_means,y=est_means)) + 
  geom_text(aes(color=model,label=genotype) ) + 
  facet_wrap(~iter,scales='free') + 
  geom_abline(intercept = 0,slope = 1,linetype='dashed') + theme_bw() + 
  geom_smooth(aes(color=model),se=F,method = 'lm')+  
  ylab('Estimated genotypic values') + 
  xlab("True genotypic values")

write.csv(dataset_summary_out_30_ref_case,'Results/Plot_inputs/dataset_summary_out_30_ref_case.csv')

# Magical factor


dataset_summary_MF <- dataset_summary %>% filter(!is.na(qh2.2st.w.qV),qh2.2st.w.qV>0) %>% select("iter","ng","model","genotype","sigma_g","sigma_e","est_means","scl_means","scl_est_means","SE","quasiSE")


dataset_summary_MF$genotype <- factor(dataset_summary_MF$genotype)

dataset_summary_MF_rank  <- dataset_summary_MF %>% filter(model!='MM')

dataset_summary_MF_metric  <- dataset_summary_MF %>% filter(model=='MM')

dataset_magical_factor <- dataset_summary_MF_rank %>% 
  left_join(dataset_summary_MF_metric[,-3],by=c('iter','ng','genotype','sigma_g','sigma_e'),suffix = c("_rank","_metric"))

dataset_magical_factor <- dataset_magical_factor %>% filter(genotype != "G1")
  
agg_vals <- dataset_magical_factor %>% group_by(iter,ng,model,sigma_g,sigma_e) %>%
  summarise(
    mean_SE_Rank=mean(SE_rank),
    median_qSE_Rank=median(SE_rank),
    mean_qSE_Rank=mean(quasiSE_rank),
    median_qSE_Rank=median(quasiSE_rank),
    mean_SE_metric=mean(SE_metric)/unique(sigma_e),
    median_qSE_metric=median(SE_metric)/unique(sigma_e),
    mean_qSE_metric=mean(quasiSE_metric)/unique(sigma_e),
    median_qSE_metric=median(quasiSE_metric)/unique(sigma_e))

mag_fact <- agg_vals %>% mutate(mag_fact_mean_SEQse=mean_qSE_Rank/mean_SE_metric,
                    mag_fact_median_SEQse=median_qSE_Rank/median_qSE_metric,
                    mag_fact_mean_qseQse=mean_qSE_Rank/mean_qSE_metric,
                    mag_fact_median_qseQse=median_qSE_Rank/median_qSE_metric)



mag_fact<- mag_fact %>% mutate(ratio=sigma_g/sigma_e)


mag_fact2 <- mag_fact %>% group_by(iter,ratio,ng,model,sigma_g,sigma_e) %>% 
  summarise(mu_mag_fact_mean_SEQse=mean(mag_fact_mean_SEQse,na.rm=T),mu_mag_fact_median_SEQse=mean(mag_fact_median_SEQse,na.rm=T),
            mu_mag_fact_mean_qseQse=mean(mag_fact_mean_qseQse,na.rm=T),mu_mag_fact_median_qseQse=mean(mag_fact_median_qseQse,na.rm=T))

write.csv(mag_fact2,'Results/Plot_inputs/mag_fact2.csv',row.names = F)

mag_fact2  %>%  #[round(mag_fact2$ratio,8) %in% c(0.18181818,0.71428571,0.5,2.25,4.50000000),] %>% 
  ggplot(aes(x=ratio,y=mu_mag_fact_median_SEQse)) +geom_point(aes(color=factor(ng))) + 
  geom_line(aes(color=factor(ng)))+theme_bw()+ylim(c(0,5))+geom_smooth(aes(color=factor(ng)),se=F)



mag_fact2  %>%  #[round(mag_fact2$ratio,8) %in% c(0.18181818,0.71428571,0.5,2.25,4.50000000),] %>% 
  ggplot(aes(x=ratio,y=mu_mag_fact_mean_SEQse)) +geom_point(aes(color=factor(ng))) + 
  geom_line(aes(color=factor(ng)))+theme_bw()+ylim(c(0,5))+geom_smooth(aes(color=factor(ng)),se=F)+
  ylab("median(qSE_rank)/(median(metric_SE)/sigma_e)")




sa <- mag_fact2 %>% filter(ratio<1,ng==100) 

median(sa$mu_mag_fact_median_SEQse,na.rm=T)^2

# Later analysis


sa <- dataset_summary_MF %>% filter(iter==7,ng==3,model=="Th",sigma_g==1.5,sigma_e==3.5)

genotype <- factor(sa$genotype)
est_values <- sa$est_means
weights <- 1/(sa$quasiSE)^2

dataset <- data.frame(genotype,est_means=est_values,w2=weights)

ds <- dataset[c('est_means','genotype','w2')]

as <- asreml(fixed = est_means ~ 1, random = ~ genotype , weights = w3,
       family = (asr_gaussian(dispersion = 1)),
       data = dataset,trace=F)
qH2(as)

as <- asreml(fixed = est_means ~ 1, random = ~ genotype , weights = w2,
             family = (asr_gaussian(dispersion = 1)),
             data = ds,trace=F)

asreml(fixed = est_means ~ 1, random = ~ genotype , weights = w2,
       family = (asr_gaussian(dispersion = 1)),
       data = dataset,trace=F)



as <- asreml(fixed = est_means ~ 1, random = ~ genotype , weights = w2,
             family = (asr_gaussian(dispersion = 1)),
             data = as,trace=F)

PEV <- function(genotype,est_values,weights){
  require(asreml)
  ds <- data.frame(genotype,est_values,weights)
  ds$genotype <- factor(ds$genotype )
  asmod <- asreml(fixed = est_values ~ 1, random = ~ genotype , weights = weights,
         family = (asr_gaussian(dispersion = 1)),
         data = ds,trace=F)
  
  BLUPl <- summary(asmod, coef=TRUE)$coef.random
  BLUPl <-as.data.frame(BLUPl)
  BLUPl <- BLUPl[grep('genotype_',rownames(BLUPl)),]
  
  BLUPl$std.error^2
  
}

# runningZeroConstant variance


ls_dataset_summary_MF<- split(dataset_summary_MF,paste(dataset_summary_MF$iter,dataset_summary_MF$ng,dataset_summary_MF$model,dataset_summary_MF$sigma_g,dataset_summary_MF$sigma_e,sep = '-'))

ls_dataset_summary_MF2 <- list()

for(i in 1:length(ls_dataset_summary_MF)){
  print(i)
   addvar <- ls_dataset_summary_MF[[i]]
   addvar$PEV <- PEV(genotype=addvar$genotype,est_values=addvar$est_means,w=1/addvar$quasiSE^2)
   ls_dataset_summary_MF2 <- list(ls_dataset_summary_MF2,addvar)
}


dataset_summary_MF2 <- dataset_summary_MF %>% group_by(iter,ng,model,sigma_g,sigma_e) %>%
  mutate(PEV=PEV(genotype=genotype,est_values=est_means,w=1/quasiSE^2))



# Analysis of standard errors


se_dataset <- read.csv("Results/Plot_inputs/summary_analysis_simulation_SE.csv")

rmvIt <- se_dataset$iter[which(se_dataset$model=='PL'& se_dataset$est_means>4.3 |se_dataset$model=='PL'& se_dataset$est_means < -3 )]

se_dataset_summ <- se_dataset %>% filter( !(iter %in% rmvIt) ) %>% 
  group_by(iter,ng,sigma_e,sigma_g,model) %>% mutate(ref_gen=value.genotype[genotype=="G1"],true_worth=(value.genotype-ref_gen)/sigma_e) %>%
  group_by(ng,genotype,sigma_e,sigma_g,model) %>%  mutate(n=n()) %>%
  summarise(mean_SE=mean(SE),sd_estimate=sd(est_means),sd_estimate2 = sqrt(sum((est_means-true_worth)^2)/(unique(n)-1)),
            median_SE=median(SE),sd_me_estimate=mad(est_means),n=n()) %>% filter(genotype!='G1'& model!="MM" )


se_dataset %>% ggplot(aes(x=genotype,y=est_means)) + geom_boxplot() +
  facet_wrap(~model)

se_dataset %>% filter(genotype=="G20",model!="MM") %>% ggplot(aes(x=est_means)) +geom_histogram()+
  facet_wrap(~model)

se_dataset %>% filter(genotype=="G20",model!="MM") %>% 
  mutate(abs_dif = abs(est_means-median(est_means)))  %>% ggplot(aes(x=abs_dif)) +geom_histogram()+
  facet_wrap(~model) + geom_vline(aes(xintercept=median(abs_dif)))+ 
  geom_vline(aes(xintercept=mad(est_means,constant = 1)))




se_dataset_summ %>% ggplot(aes(x=sd_estimate,y=mean_SE)) + geom_point(aes(color=model)) +
  geom_abline(intercept=0,slope=1,linetype='dashed')+geom_text(aes(label=genotype))




se_dataset_summ %>% ggplot(aes(x=sd_estimate2,y=mean_SE)) + geom_point(aes(color=model)) +
  geom_abline(intercept=0,slope=1,linetype='dashed')

se_dataset_summ %>% ggplot(aes(x=median_SE,y=sd_me_estimate)) + geom_point(aes(color=model)) +
  geom_abline(intercept=0,slope=1,linetype='dashed')+geom_text(aes(label=genotype))




hist(se_dataset[se_dataset$model=="Th" & se_dataset$genotype=="G20",]$est_means)

hist(se_dataset[se_dataset$model=="Th" & se_dataset$genotype=="G20",]$SE)

# Analising outlayers

dataset_summary_out <- dataset_summary %>%  filter(model!="MM") %>% 
  mutate(Estimates_val=ifelse(model=="Th","Th worth values","PL Log worth values"))

dataset_summary_out$Outalyer_value <- abs(dataset_summary_out$est_means) > 4


dataset_summary_out_mod <- dataset_summary_out %>% group_by(iter,ng,sigma_e,sigma_g,model,Singular,
                                 Experiment0,Experiment) %>% summarise(indeterminate = sum(Outalyer_value)>0)

dim(dataset_summary_out_mod)

dataset_summary_out_mod <- dataset_summary_out_mod %>% filter(ng %in% c(3,10,30,100))

#4*14*11*30*2

dataset_summary_out_mod2 <- dataset_summary_out_mod %>% mutate(sigma_std=sigma_g/sigma_e) %>%
  group_by(ng,model,sigma_std) %>% summarise( tot_singular = sum(Singular),tot_indet = sum(indeterminate), n=n())


dataset_summary_out_mod2 <- dataset_summary_out_mod2 %>% group_by(ng,model,sigma_std) %>%
  mutate(perc_sing = tot_singular/n*100,perc_indet=tot_indet/n*100)


sum(dataset_summary_out_mod2$n)

write.csv(dataset_summary_out_mod2,'Results/Plot_inputs/outlayers_inderterminates.csv')

bx <- dataset_summary_out %>%  ggplot(aes(x='',y=est_means)) +
  geom_boxplot() + facet_wrap(~Estimates_val)+ ylim(c(-8,8))+theme_bw()+ylab('Estimates')



ggsave('bx.png',bx,heigh=4,width=5.5)

sing_data <- unique( dataset_summary_out %>% select(iter, ng, model,sigma_e, sigma_g,model ,Singular ) ) %>% mutate(Ratio = sigma_g/sigma_e)



sing_data_agg <- sing_data %>% group_by( ng, sigma_e, sigma_g,model,Ratio ) %>% summarise(tot_sing=sum(Singular),tot=n())

sing_data2_full <- split(sing_data_agg,with(sing_data_agg,paste(ng,model))) 

ls_sing_data2_full <- lapply(sing_data2_full,function(x){data<- x[order(x$Ratio,decreasing = F),]; data$c_tot_sing <- cumsum(data$tot_sing); data$c_tot <- cumsum(data$tot) ;data })

cum_dis <- do.call(rbind,ls_sing_data2_full)

gssc <- cum_dis %>% mutate(Percentaje=c_tot_sing/(770*30)*100) %>% filter(ng!=5)  %>%ggplot(aes(x=Ratio,y=Percentaje,color=model,group=paste(ng,model))) +geom_point(aes(shape=factor(ng))) + geom_line() + ylab("Percetaje of trials with singular variance")+theme_bw()

ggsave('gssc.png',height = 4,width = 5)


# %>% mutate(tot_sing=cumsum(tot_sing), ctot=cumsum(sing_data1$tot))


sing_data1 <- sing_data_agg %>% filter(model=='Th' ) %>% arrange(Ratio)

sing_data1$ctot <- cumsum(sing_data1$tot)

sing_data1$ctot_sing <- cumsum(sing_data1$tot_sing)

sing_data2 <- sing_data_agg %>% filter(model=='PL' ) %>% arrange(Ratio)

sing_data2$ctot <- cumsum(sing_data2$tot)

sing_data2$ctot_sing <- cumsum(sing_data2$tot_sing)

sing_data2_full <- rbind(sing_data1,sing_data2)

sing_data2_full <- sing_data2_full %>% mutate(Percentaje = ctot_sing/ctot) %>% filter(ng %in% c(3,10,30,100))

sing_data2_full %>% group_by(ng,model,Ratio) %>% summarise( Pect_missing= sum(Percentaje)) %>% ggplot(aes(x=Ratio,y=Pect_missing,color=model)) +geom_point() + geom_line() + facet_wrap(~ng)




data.frame(sing_data %>% group_by(model) %>% summarise(tot_sim=sum(Singular_count) ,tot=sum(n) , perc=tot_sim/tot*100))[,c(1,4)]

boxplot(dataset_summary$est_means)

sing_data2 <- sing_data %>% mutate(Ratio = sigma_g/sigma_e) %>% select(ng,Ratio,model,Singular_count,n)

sing_data2 <- sing_data2 %>% arrange(Ratio) %>%  group_by(model,Ratio,n) %>% mutate(cum_sum=cumsum(Singular_count),sumn=cumsum(n))

# Evaluating the standard errrors
ref_case <- dataset_summary %>% filter(ng==30,iter==10,sigma_g==0.8,sigma_e==0.9) 

ref_case %>% filter(model=='MM' & genotype!='G11') %>% ggplot(aes(x=SE,y=quasiSE)) + geom_point()+
  theme_bw() + geom_abline(intercept = 0,slope=1)


qSE_summary <- dataset_summary %>% filter(iter==10)%>% group_by(ng,sigma_e,sigma_g,model) %>% 
  summarise(mean_qSE=mean(quasiSE,na.rm=T),SD_mean_qSE=sd(quasiSE,na.rm=T),
            mean_SE=mean(SE,na.rm=T),SD_mean_se=sd(SE,na.rm=T),
            mean_cont_SE=mean(SE_cont[genotype!='G1'],na.rm=T),SD_mean_cont_SE=sd(SE[genotype!='G1'],na.rm=T),
            median_cont_SE=median(SE_cont[genotype!='G1'],na.rm=T),
            n=n())


qSE_summary2 <- qSE_summary %>% filter(sigma_e==0.9, sigma_g %in% c(0.8,2))

qSE_summary2 %>% ggplot(aes(x=ng,y=mean_qSE,colour = model))+geom_point(aes(shape=factor(sigma_g)))+ylim(0,1.5)+
  geom_errorbar(aes(ymax=mean_qSE+1.96*SD_mean_qSE/sqrt(n),
                    ymin=mean_qSE-1.96*SD_mean_qSE/sqrt(n)))+theme_bw()

as <- dataset_summary %>% filter(ng==100,iter==11,sigma_g %in% c(0.8),sigma_e==0.9)
as
mmqse <- qSE_summary %>% select(ng:mean_qSE,mean_cont_SE,median_cont_SE) %>% filter(model=='MM') %>% select(-'model') %>% rename(metric_mean_qse='mean_qSE',
                                                                                                                                 metric_mean_sed='mean_cont_SE',
                                                                                                                                 metric_median_cont_SED="median_cont_SE")

rnqse <- qSE_summary %>% select(ng:mean_qSE,mean_cont_SE,median_cont_SE) %>% filter(model!='MM')

qSE_summary %>% filter(ng==100,sigma_g %in% c(0.8),sigma_e==0.6)

fullqSE <- left_join(rnqse,mmqse)

fullqSE <- fullqSE %>% mutate(nfactor=mean_qSE*sigma_e/(metric_mean_qse),
                              nfactor2=mean_cont_SE*sigma_e/(metric_mean_sed),
                              nfactor3=median_cont_SE*sigma_e/(metric_median_cont_SED))

fullqSE %>% filter(ng>=10) %>%ggplot(aes(x=ng,y=nfactor,color=model)) +
  geom_point()+geom_line()+facet_grid(sigma_e~sigma_g,scales= 'free') + ylim(0,3)


fullqSE %>% filter(ng==100) %>% group_by(model,sigma_e,sigma_g) %>% 
  summarise(mean_n_factor=mean(nfactor)) %>% mutate(ratio=sigma_g/sigma_e) %>%
  ggplot(aes(y=mean_n_factor,x=ratio,color=model))+geom_point()+geom_line()+xlim(0,4)+ylim(0,2.5)#500


fullqSE %>% filter(ng==100) %>% group_by(model,sigma_e,sigma_g) %>% 
  summarise(mean_n_factor=mean(nfactor2)) %>% mutate(ratio=sigma_g/sigma_e) %>%
  ggplot(aes(y=mean_n_factor,x=ratio,color=model))+geom_point()+geom_line()+xlim(0,4)+ylim(0,2.5)#500


fullqSE %>% filter(ng==100) %>% group_by(model,sigma_e,sigma_g) %>% 
  summarise(mean_n_factor=mean(nfactor3)) %>% mutate(ratio=sigma_g/sigma_e) %>%
  ggplot(aes(y=mean_n_factor,x=ratio,color=model))+geom_point()+geom_line()+xlim(0,4)+ylim(0,2.5)#500


SED_analysis <- fullqSE %>% group_by(ng,model,sigma_e,sigma_g) %>% 
  summarise(mean_n_factor=mean(nfactor,na.rm=T)) %>% mutate(ratio=sigma_g/sigma_e) %>% filter(round(ratio,6) %in% c(0.02500000, 0.8, 3.75000000)) 


SED_analysis %>% 
  ggplot(aes(y=mean_n_factor,x=ng ,color=model))+geom_point()+geom_line()+
  facet_grid(ratio~. ) +theme_bw()+ylim(1,3)


SED_analysis <- fullqSE %>% group_by(ng,model,sigma_e,sigma_g) %>% 
  summarise(mean_n_factor=mean(nfactor2,na.rm=T)) %>% mutate(ratio=sigma_g/sigma_e) %>% filter(round(ratio,6) %in% c(0.02500000, 0.8, 3.75000000)) 

SED_analysis %>% 
  ggplot(aes(y=mean_n_factor,x=ng ,color=model))+geom_point()+geom_line()+
  facet_grid(ratio~. ) +theme_bw()+ylim(1,2)

SED_analysis <- fullqSE %>% group_by(ng,model,sigma_e,sigma_g) %>% 
  summarise(mean_n_factor=mean(nfactor3,na.rm=T)) %>% mutate(ratio=sigma_g/sigma_e) %>% filter(round(ratio,6) %in% c(0.02500000, 0.8, 3.75000000)) 

SED_analysis %>% 
  ggplot(aes(y=mean_n_factor,x=ng ,color=model))+geom_point()+geom_line()+
  facet_grid(ratio~. ) +theme_bw()+ylim(1,5)


# Check number of parameters

parms_consol <- dataset_summary[,2:4]

table(parms_consol$ng) #5

table(parms_consol$sigma_e) #14

table(parms_consol$sigma_g) # 11

5*14*11*30*3*20 #ng*se_e*sg*iter*mod*meth*gen

# add label text

#--------------------------------RUN HERE--------------------------------------

dataset_summary <- dataset_summary %>% mutate(ps = paste(ng, sigma_g, sigma_e, iter, model,genotype,sep='-'))

# Reference case - average parameters 

dataset_summary %>% filter(ng==100,sigma_e==0.9,sigma_g==0.8,model!="MM") %>% 
  group_by(iter) %>% ggplot(aes(x=scl_means,y=value.genotype,color=model)) + 
  geom_point() +geom_abline(intercept = 0,slope = 1)+ 
  facet_wrap(~iter,scales = 'free')+ geom_smooth(method = 'lm',se=F) + 
  theme_bw()


reference_cases_recovering_scale <- dataset_summary %>% 
  filter(ng==100,sigma_e==0.9,sigma_g==0.8,iter==22,model!="MM")

write.csv(reference_cases_recovering_scale,
          'Results/Plot_inputs/recovering_scale.csv',row.names = F)


#------------------------ General pre-processing -------------------------------

# Adding a column with reference genotypic value


# Creating a subset with only ranking outputs

dataset_summary1 <- dataset_summary %>% filter(model != 'MM') # Remove metric data results

dataset_summary1 <- droplevels(dataset_summary1)

#----------------------- Select  valid estimates  ------------------------------

remove_outlyrs_singular <- dataset_summary1 %>% ungroup() %>% filter(Singular==T) %>% 
  dplyr::select(Experiment0) %>% unique() %>% as.data.frame()


remove_outlyrs_negative_SE <- dataset_summary %>% ungroup() %>% filter(is.na(SE)) %>% 
  dplyr::select(Experiment0) %>% unique() %>% as.data.frame()


# Clean dataset of experiment with a least one singular case

dataset_summary_noS0 <- dataset_summary %>% 
  mutate( Exp_valid = case_when(Experiment0 %in% remove_outlyrs_singular$Experiment0 ~ "Singular",
                                !(Experiment0 %in% remove_outlyrs_singular$Experiment0) & 
                                  Experiment0 %in% remove_outlyrs_negative_SE$Experiment0 ~ "Negative SE" ,.default="Valid") # There was NA but not negaviate SE
  ) 

table(dataset_summary_noS0$Exp_valid)

ds_sm_noS <- unique(dataset_summary_noS0[c('Experiment0','model','Exp_valid')])

ds_sm_noS <- ds_sm_noS[ds_sm_noS$Exp_valid!='Valid',]

table(ds_sm_noS$Experiment0,ds_sm_noS$model,ds_sm_noS$Exp_valid)


length(unique(ds_sm_noS$Experiment0))/23100*100 # percentage of missing data

write.csv(dataset_summary_noS0,'Results/full_experiment_valid_novalidcases.csv')

# Remove experiment from the original dataset

dataset_summary_noS <- dataset_summary_noS0 %>% filter(Exp_valid=="Valid")

summary(dataset_summary_noS)



# Variation ranks

dataset_summary_noS %>% filter(model!="MM") %>% select(c(iter:sigma_g,genotype,est_means,model)) %>% group_by(iter,ng, sigma_e, sigma_g,model) %>% 
  summarise(mean_observ=mean(est_means),mean_Se=1.96*sd(est_means)/sqrt(n()),n=n())

write.csv(dataset_summary_noS,'Results/dataset_summary_noS.csv',row.names = F)

dataset_summary_noS <- read.csv('Results/dataset_summary_noS.csv')


##------------------------  Correlation - analysis -----------------------------

# index_summary_0 <- dataset_summary %>% group_by(iter, ng, sigma_e, sigma_g , model) %>% 
#   summarise(cor = cor(value.genotype, est_means),
#             R2= summary(lm(value.genotype~est_means))$r.squared,
#             #cor_sp = cor(value.genotype, scl_means,method = 'spearman'),
#             cor2 = cor(value.genotype, est_means,method = 'kendall'),
#             RMSE=sqrt(mean((value.genotype-scl_means)^2)),
#             Bias= mean(sum(scl_means-scl_est_means)),
#   ) %>% mutate(Ratio=sigma_g/sigma_e)


# Correct the scale of plackett-luce estimates

#--------------------------------RUN HERE--------------------------------------

dataset_summary <- dataset_summary %>% ungroup() %>% mutate(id=1:nrow(dataset_summary)) %>% relocate(id)



re_scale_PL <- dataset_summary %>% filter(model!="MM") %>%
  select(iter,ng, sigma_e, sigma_g,model,est_means,genotype,quasiSE) %>%
  pivot_wider(values_from =  c(est_means,quasiSE) ,names_from = c(model))

re_scale_PL <- re_scale_PL %>% group_by(iter,  ng, sigma_e, sigma_g) %>% 
  mutate(slope_th_pl=sum(est_means_Th*est_means_PL)/sum(est_means_PL^2)) 

re_scale_PL <- re_scale_PL %>% mutate(full_scl_est_means=est_means_PL*slope_th_pl,
                                      full_quasiSE_scl=quasiSE_PL*slope_th_pl)

data.frame(head(re_scale_PL))

# robustness of scale factor


# Integration of the scaled variable

re_scaled_vars_PL <- re_scale_PL %>% select(iter, ng, sigma_e, sigma_g ,  
                                            genotype,full_scl_est_means ,
                                            full_quasiSE_scl,slope_th_pl) %>% mutate(model='PL') 

#--------------------------------RUN HERE--------------------------------------

dataset_summary2 <- dataset_summary %>% left_join(re_scaled_vars_PL,by=c('iter','ng','sigma_e','sigma_g','genotype','model'))

dataset_summary2 <- dataset_summary2 %>% mutate(full_scl_est_means=ifelse(is.na(full_scl_est_means),scl_est_means,full_scl_est_means),
                                                full_quasiSE_scl=ifelse(is.na(full_quasiSE_scl),quasiSE,full_quasiSE_scl))

data.frame(head(dataset_summary2,20))
quasiSE_scl

table(table(paste0(dataset_summary$iter,dataset_summary$ng,dataset_summary$sigma_e,dataset_summary$sigma_g,dataset_summary$genotype)))


# RMSE and bias analysis 
index_summary_0 <- dataset_summary2 %>% group_by(iter, ng, sigma_e, sigma_g , model) %>% 
  filter(genotype!="G1") %>%
  summarise(cor = cor(scl_true_means, full_scl_est_means,method = 'kendall'),
            #R2= summary(lm(value.genotype~est_means))$r.squared,
            #cor_sp = cor(value.genotype, scl_means,method = 'spearman'),
            #cor2 = cor(value.genotype, est_means,method = 'kendall'),
            RMSE=sqrt(mean((scl_true_means-full_scl_est_means)^2)),
            Bias= mean((full_scl_est_means-scl_true_means)),
  ) %>% mutate(Ratio=sigma_g/sigma_e)

write.csv(index_summary_0,"Results/Plot_inputs/index_summary_RMSE_CCOR_R2_V3_fulldataset.csv",row.names = F)

#index_summary_0 <- read.csv("Results/Plot_inputs/index_summary_RMSE_CCOR_R2_V2_fulldataset.csv")

performance_Values_0 <- index_summary_0 %>% pivot_longer(-c(iter,ng,sigma_e,sigma_g,model,Ratio),
                                                         names_to = 'variable',values_to = 'value') %>%
  group_by( ng, sigma_e, sigma_g, model, Ratio, variable) %>% 
  summarise(mu_value=mean(value),SE=1.96*sd(value)/sqrt(n()),n=n())

performance_Values_0 %>% filter(ng %in% c(3,30,100) & !(variable %in% c('Bias','RMSE'))) %>%
  ggplot(aes(x=Ratio,y=mu_value,group=interaction(model,ng))) + 
  geom_point(aes(colour = model,shape = factor(ng)),alpha=0.6) +
  geom_errorbar(aes(ymax = mu_value+SE,ymin=mu_value-SE,colour = model))+
  stat_smooth(geom="line",aes(colour = model),se=F,alpha = 0.5,show.legend = F)+
  facet_grid(variable~.,scales = "free")+theme_bw()


#------------ Two step-approach to estimate genotypic varaince  new version---------

dataset_summary2 <- dataset_summary2 %>% group_by(iter,ng, sigma_e, sigma_g,model) %>%
  mutate(SE_median = median(SE_cont[genotype!='G1']),h2_SEd2 = sd_g_w.qV^2/(sd_g_w.qV^2+SE_median^2/2))

dataset_summary3_agg <- dataset_summary2 %>% group_by( iter, ng, sigma_e, sigma_g,model) %>%
  mutate(error_w = ifelse(model=="MM",SE,quasiSE), invalid_SE = sum(is.na(SE))>0, invalid_qSE=sum(is.na(quasiSE))>0 )



dataset_summary3_agg <- dataset_summary3_agg %>% select(iter,ng,sigma_e,sigma_g,model,
                                                Singular,sd_g_w.nW,sd_g_w.qV,
                                                sd_g_w.GM,sd_g_sc.scl,sd_g_w.Cotr,qh2.2st.w.qV,
                                                qh2.2st.w.qV,qh2.2st.w.GM,
                                                h2_broad_metric1, h2_broad_metric2,
                                                h2_metric_cullis,invalid_SE,invalid_qSE,h2_SEd2,SE_median) %>% unique()


# bring slope

slope_th_pl <- re_scale_PL %>% select(iter,  ng, sigma_e, sigma_g,slope_th_pl) %>% unique()




dataset_summary3_agg <- dataset_summary3_agg %>% left_join(slope_th_pl,by=c('iter','ng','sigma_e','sigma_g' ))


dataset_summary3_agg <- dataset_summary3_agg %>% group_by(iter,ng,sigma_e,sigma_g,model,Singular,h2_metric_cullis) %>%
  mutate(sg_g_estimate_scl=ifelse(model=='MM',sd_g_w.GM/sigma_e,sd_g_w.qV),
         sg_g_estimate_scl2=ifelse(model=='MM',sd_g_w.qV/sigma_e,sd_g_w.qV),
         sg_gse_estimate_scl2=ifelse(model=='MM',sd_g_w.Cotr/sigma_e,sd_g_w.Cotr),
         heritability_fnl=ifelse(model=='MM',qh2.2st.w.GM,qh2.2st.w.qV),heritability_cullisEstQV=qh2.2st.w.qV) %>%
  mutate(sg_g_estimate_scl =ifelse(model=='PL',sg_g_estimate_scl*slope_th_pl,sg_g_estimate_scl),
         sg_g_estimate_scl2=ifelse(model=='PL',sg_g_estimate_scl2*slope_th_pl,sg_g_estimate_scl2),
         sg_gse_estimate_scl2=ifelse(model=='PL',sd_g_w.Cotr*slope_th_pl,sd_g_w.Cotr)) 


write.csv(dataset_summary3_agg,'Results/Plot_inputs/estimation_genotypic_variance.csv')


genotypic_se_g <- dataset_summary3_agg %>% mutate(true_sigma_sc=sigma_g/sigma_e) %>% 
  group_by(ng,sigma_g, sigma_e, true_sigma_sc, model) %>% filter(!is.na(sg_g_estimate_scl)) %>% 
  summarise(mu_sd_g_wQV = mean(sg_g_estimate_scl,na.rm=T),
            SE_mu_sd_g_wQV = sd(sg_g_estimate_scl,na.rm=T)/sqrt(n()) ,n=n())

genotypic_se_g %>% ggplot(aes(x=true_sigma_sc,y=mu_sd_g_wQV,color=model)) +
  geom_point()+ geom_line() + facet_wrap(~ng)

# Heritability estimation




#------------ Two step-approach to estimate genotypic varaince old version----------



dataset_summary_genetic_pars <- dataset_summary_noS %>% dplyr::select(c(iter:h2_metric_cullis,model,sd_g_wQV,Experiment)) %>%
  unique()

write.csv(dataset_summary_genetic_pars,'Results/dataset_summary_genetic_pars.csv',row.names = F)

dataset_summary_genetic_pars <- read.csv('Results/dataset_summary_genetic_pars.csv')

outalayers_QSE <- dataset_summary_noS %>% mutate(outlayer1_QSE = quasiSE > 1.191)

exp_outalayers_QSE <- outalayers_QSE %>% dplyr::filter(outlayer1_QSE) %>% select(Experiment) %>% unique()

dataset_summary_genetic_pars_outRemoved <- dataset_summary_genetic_pars %>% dplyr::filter(!(Experiment %in%  exp_outalayers_QSE[,1]))

# --------------------  Including outlayers ------------------------------------

genotypic_se_g <- dataset_summary_genetic_pars %>% 
  mutate(sd_g_wQV_scaled = ifelse(model=='MM',sd_g_wQV,sd_g_wQV*sigma_e)) %>% 
  group_by(ng,sigma_g, sigma_e,  model) %>% 
  summarise(mu_sd_g_wQV = mean(sd_g_wQV,na.rm=T),
            SE_mu_sd_g_wQV = sd(sd_g_wQV,na.rm=T)/sqrt(n()) ,
            mu_sd_g_wQV_scaled= mean(sd_g_wQV_scaled,na.rm=T),
            SE_mu_sd_g_wQV_scaled=sd(sd_g_wQV_scaled,na.rm=T)/sqrt(n()),n=n())

write.csv(genotypic_se_g,'Results/Plot_inputs/estimation_genotypic_variance.csv',row.names = F)                          

genotypic_se_g %>% filter(model!='MM' & sigma_e == 0.8) %>% ggplot(aes(x=mu_sd_g_wQV,y=sigma_g))+
  geom_smooth(method = 'lm',se=F,aes(shape=model),color='black',linewidth=0.7) + 
  geom_point(aes(shape=model)) + 
  geom_errorbar(aes(xmax=mu_sd_g_wQV+SE_mu_sd_g_wQV,xmin=mu_sd_g_wQV-SE_mu_sd_g_wQV),width=0.02) + 
  facet_wrap(~ng) +theme_bw() + geom_abline(intercept = 0, slope = 0.8,linetype= 'dashed')


genotypic_se_g %>%ggplot(aes(x=mu_sd_g_wQV_scaled,y=sigma_g,group = interaction(model,factor(sigma_e),factor(ng))))+
  geom_point(aes(color=model,shape=factor(sigma_e)),alpha=0.7)+ 
  geom_errorbar(aes(xmin=mu_sd_g_wQV_scaled-SE_mu_sd_g_wQV_scaled,xmax=mu_sd_g_wQV_scaled+SE_mu_sd_g_wQV_scaled,colour = model),width=0.01, alpha=0.4, linewidth=0.01) +
  geom_line(aes(color=model),alpha=0.7)+
  geom_abline(intercept =  0,slope=1,linetype=2)+
  facet_wrap(~ng)+theme_bw()+ggtitle(label = 'include out')


# differences analysis (scatterplot)

lapply(1:30,function(w){  as <- dataset_summary_genetic_pars %>% filter(iter==w) %>% 
  ggplot(aes(y=sigma_g,x=sd_g_wQV)) + geom_point(aes(color=model)) + facet_grid(ng~sigma_e)+
  stat_poly_eq(aes(label =  paste(..eq.label.., ..rr.label.., 'rho','`=`',round(sqrt(r.squared),2),sep = "~~~~"),colour = model),formula = y ~ x,parse = TRUE,size=2.7)+
  geom_smooth(method = 'lm',se=F,aes(color=model))+theme_bw()
ggsave(paste0('Results/Plot_genotype_variance/iter_',w,'.png'),as,height = 14,width =16 )
}
)

# pearson correlation 

cor_analysis_se_g <- dataset_summary_genetic_pars %>% group_by(iter, ng,sigma_e,  model) %>% 
  summarise(cor_est_sg_truesg=cor(sigma_g,sd_g_wQV,use="complete.obs"),
            R2_est_sg_truesg=cor_est_sg_truesg^2,n=n()) 

sd_g_wQV_summ <- dataset_summary_genetic_pars %>% group_by(ng,sigma_e,sigma_g,model) %>% 
  summarise(mu_sd_g_wQV = mean(sd_g_wQV,na.rm = T),l1_sd_g_wQV=mu_sd_g_wQV - 1.96*sd(sd_g_wQV,na.rm = T)/n(),
            l2_sd_g_wQV=mu_sd_g_wQV + 1.96*sd(sd_g_wQV,na.rm = T)/n(),n=n())


# sd_g_wQV_summ %>% ggplot(aes(mu_sd_g_wQV,sigma_g))+geom_point() + facet_grid(ng~sigma_e)
# 
# 
rho_genotic_manes <- sd_g_wQV_summ %>% group_by(ng,sigma_e,model) %>%
  summarise(rho_mean = cor(mu_sd_g_wQV,sigma_g),rho_l1=cor(l1_sd_g_wQV,sigma_g),
            rho_2=cor(l2_sd_g_wQV,sigma_g))

# Scarrterplot for correlation btwn means of person correlation between genotypic variance estimates and true variance

ge_sigm <- rho_genotic_manes %>% ggplot(aes(x=sigma_e,y=rho_mean,group=interaction(model,factor(ng))))+geom_point(aes(color=factor(ng),shape = model))+
  geom_line(aes(color=factor(ng)))+
  ylab("Pearson corrletion coef. between actual genotypic variance \nvs mean of estimates genotypic variance")+xlab("Sigma E")+theme_bw()#+ylim(0.99,1)

ge_sigm

ggsave('Results/ge_sigm.png',ge_sigm,height = 4.5,width = 6)

save_cor <- dataset_summary_genetic_pars %>% 
  left_join(cor_analysis_se_g,by = c('iter', 'ng','sigma_e',  'model'))

write.csv(save_cor,'Results/final_correlation.csv')

# Scarrterplot for correlation btwn genotypic variance estimates and true variance

cor_analysis_se_g %>% ggplot(aes(x=factor(sigma_e),y=cor_est_sg_truesg)) +
  facet_wrap(~ng)  +geom_jitter(aes(color=model),width = 0.1, height = 0,alpha=0.6)+
  theme_bw()+
  ylim(0.4,1)

data.frame(dataset_summary_noS %>% filter(iter==13,ng== 30, sigma_e==0.6))

# Fail trys to show that sigma_e affects the correlation between genotypic estimates variance and true estiamate variance

cor_analysis_se_g2 <-  cor_analysis_se_g %>% group_by(ng,model,sigma_e) %>% 
  summarise(mu_cor_est_sg_truesg=mean(cor_est_sg_truesg),SE_mu_cor_est_sg_truesg=1.96*sd(cor_est_sg_truesg),
            mu_R2_est_sg_truesg=mean(R2_est_sg_truesg),SE_R2_est_sg_truesg=1.96*sd(R2_est_sg_truesg))


cor_analysis_se_g2 %>% ggplot(aes(x=sigma_e,y=mu_cor_est_sg_truesg,group = interaction(model,factor(ng)))) + 
  geom_point(aes(color=factor(ng),shape=model)) + 
  geom_errorbar(aes(ymax = mu_cor_est_sg_truesg+SE_mu_cor_est_sg_truesg,
                    ymin = mu_cor_est_sg_truesg-SE_mu_cor_est_sg_truesg,
                    color=factor(ng)))+geom_line(aes(color=factor(ng)))+
  facet_wrap(~model)+
  ylim(0.75,1)+theme_bw()

genotypic_se_g %>%  ggplot(aes(x=mu_sd_g_wQV,y=sigma_g,group = interaction(model)))+
  geom_point(aes(color=model,shape=factor(sigma_e)),alpha=0.7)+ 
  geom_errorbar(aes(xmin=mu_sd_g_wQV-SE_mu_sd_g_wQV,xmax=mu_sd_g_wQV+SE_mu_sd_g_wQV,colour = model),width=0.01, alpha=0.4, linewidth=0.01) +
  geom_line(aes(color=model),alpha=0.7)+
  geom_abline(intercept =  0,slope=1,linetype=2)+
  facet_wrap(~ng+sigma_e)+theme_bw()+ggtitle(label = 'include out')


genotypic_se_g <- read.csv('Results/Plot_inputs/estimatimation_genotypic_variance.csv')

genotypic_se_g %>%  ggplot(aes(x=mu_sd_g_wQV_scaled,y=sigma_g,group = interaction(model,factor(sigma_e),factor(ng))))+
  geom_point(aes(color=model,shape=factor(sigma_e)),alpha=0.7)+ 
  geom_errorbar(aes(xmin=mu_sd_g_wQV_scaled-SE_mu_sd_g_wQV_scaled,xmax=mu_sd_g_wQV_scaled+SE_mu_sd_g_wQV_scaled,colour = model),width=0.01, alpha=0.4, linewidth=0.01) +
  geom_line(aes(color=model),alpha=0.7)+
  geom_abline(intercept =  0,slope=1,linetype=2)+
  facet_wrap(~ng)+theme_bw()+ggtitle(label = 'include out')

# ------- correlation of genotypic variance--------

genotypic_se_g %>% aes(sigma_g,mu_sd_g_wQV_scaled)

# Removing outlayers

genotypic_se_g_outRemoved <- dataset_summary_genetic_pars_outRemoved %>%
  group_by(ng,sigma_g, sigma_e,  model) %>%
  summarise(mu_sd_g_wQV= mean(sd_g_wQV),SE_sd_g=sd(sd_g_wQV))


genotypic_se_g_outRemoved %>% 
  mutate(mu_sd_g_wQV_scaled = ifelse(model=='MM',mu_sd_g_wQV,mu_sd_g_wQV*sigma_e),
         SE_sd_g_wQV_scaled = ifelse(model=='MM',SE_sd_g,SE_sd_g*sigma_e))  %>% ggplot(aes(x=mu_sd_g_wQV_scaled,y=sigma_g,group = interaction(model,factor(sigma_e))))+
  geom_point(aes(color=model,shape=factor(sigma_e)),alpha=0.7)+ geom_errorbar(aes(ymin=mu_sd_g_wQV_scaled-SE_sd_g,ymax=mu_sd_g_wQV+SE_sd_g,colour = model)) +
  geom_line(aes(color=model),alpha=0.7)+
  geom_abline(intercept =  0,slope=1,linetype=2)+
  facet_wrap(~ng)+theme_bw()+ggtitle(label = 'include out')


# Estimation of heritability

dataset_summary_noS <- read.csv('Results/dataset_summary_noS.csv')

# Heritability based in de difference

dataset_summary_noS2 <- dataset_summary_noS %>% group_by(iter,ng,model,sigma_e,sigma_g) %>% 
  mutate(mean_var_cont = mean(SE_cont[genotype!="G1"]^2,na.rm=T),
         median_var_cont = median(SE_cont[genotype!="G1"]^2,na.rm=T),
         h2_diffbased1=sd_g_wQV^2/(sd_g_wQV^2+mean_var_cont/2),
         h2_diffbased2=sd_g_wQV^2/(sd_g_wQV^2+median_var_cont/2))

dataset_summary_genetic_pars_h2 <- dataset_summary_noS2 %>% dplyr::select(c(iter:h2_metric_cullis,model,qh2.2st.2.0:qh2.2st.scl,h2_diffbased1,h2_diffbased2,Experiment))

dataset_summary_genetic_pars_h2a <- dataset_summary_genetic_pars_h2  %>%   dplyr::filter(h2_metric_cullis>0 & qh2.2st.wQV>0)

dataset_summary_genetic_pars_h2 %>% ggplot(aes(x=h2_metric_cullis,y=qh2.2st.wQV)) +
  geom_point(aes(color=model)) + geom_abline(intercept =  0,slope=1,linetype=2) + facet_wrap(~ng)

dataset_summary_genetic_pars_h2a %>% ggplot(aes(x=qh2.2st.wQV,y=h2_metric_cullis)) +
  geom_smooth(aes(color = model),method = 'lm',se=F)+
  geom_point(aes(color=model),alpha=0.4) + geom_abline(intercept =  0,slope=1,linetype=2) + 
  facet_wrap(~ng) + theme_bw()

write.csv(dataset_summary_genetic_pars_h2a,'Results/Plot_inputs/cullis_heritability_analisis.csv',row.names = F)

ggsave('Results/Paper_figures_V3/Figure_sg_sd_g.png',seg,width =8 ,height = 3)


# Observed dataset

observed_dataset <- read.csv("Results/Plot_inputs/genotypic_means_observed_dataset.csv")

observed_dataset$crop <- factor(observed_dataset$crop,levels = c("Cassava","Groundnut",
                                                                 "Early-maturing maize",
                                                                 
                                                                 "Intermerdia-maturing maize"))



obs_dataset_results <- observed_dataset %>% ggplot(aes(x=estimate_blues,y=estimate_yieldrank)) +  
  geom_smooth(method = 'lm',se=F) + facet_wrap(~crop,scales = 'free',nrow=1,ncol=4)+
  xlab(expression(paste('Genotypic values metric-based (t  ha'^"-1",")"))) +ylab('Genotypic values ranked-based')+
  geom_point()  + geom_text_repel(  #+ geom_text(aes(label = genotype))
    aes(label=genotype),
    box.padding = unit(0.2, "lines"),size=2,max.overlaps = 25
    
  )+ 
  stat_cor(
    #aes(label=paste( signif(..r.., 3),  signif(..p.., 3))),
    cor.coef.name = 'rho',
    method = "pearson",
    label.x.npc = 0.001, # Adjust placement to the right
    label.y.npc = 0.91 ,# Adjust placement to the top
    p.accuracy = 0.001, r.accuracy = 0.01
  )+
  stat_cor(
    #aes(label = paste("Coef.Pearson:", ..r.label..,"p-val",..p..)),
    cor.coef.name = 'tau',
    method = "kendall",
    label.x.npc = 0.001, # Adjust placement to the right
    label.y.npc = 0.99 ,# Adjust placement to the top#0.99
    p.accuracy = 0.001, r.accuracy = 0.01
  )+
  # annotate('text',label = "r= 0.85",x = Inf, y = Inf,
  #          hjust = 1, vjust = 1  
  #          ) +
  theme_bw()



ggsave('Results/Paper_figures_V5/observed_dataset_results.png',obs_dataset_results,height = 4 ,width =  13)






