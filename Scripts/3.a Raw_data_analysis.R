
#

# 10,30,(100),500



summary_analysis <- read.csv('Results/summary_analysis_iter_summarised.csv')

# Ratio Genotypic variance / Residual variance

summary_analysis %>% mutate(ratio=sigma_g/sigma_e) %>% ggplot(aes(x=ratio,y=muSlope)) +
  geom_point(aes(colour = model,shape=factor(ng))) + geom_line(aes(colour = model,shape=factor(ng))) +theme_bw()

# ------------------------- RMSE analysis --------------------------------------


Rmse_analysis2 <- summary_analysis %>% group_by(ng,sigma_g,sigma_e,model) %>% 
  summarise(mu_RMSE=mean(muRMSE),L1 = mean(muRMSE)-se_RMSE,
            L2 =  mean(muRMSE)+se_RMSE)


Rmse_analysis2$sg_se <- paste('sg_g=',Rmse_analysis2$sigma_g,',sg_e=', Rmse_analysis2$sigma_e,sep='')

#Rmse_analysis2 <- Rmse_analysis2[Rmse_analysis2$sg_se %in% c("sg_g=0.2,sg_e=0.6","sg_g=0.8,sg_e=0.9","sg_g=1.2,sg_e=1.6"),]

Rmse_analysis2 <- Rmse_analysis2[Rmse_analysis2$sg_se %in% c("sg_g=0.2,sg_e=0.9","sg_g=0.8,sg_e=0.9","sg_g=1.2,sg_e=0.9"),]

custom_colors <- c("Th" = "#FEA800",  # Deep Red
                   "PL" = "#00CCBE",
                   "MM" = "#00A9FF")  # Crimson

Rmse_analysis2 <- Rmse_analysis2 %>% rename(Model=model)


shapes <- expression('  '*sigma[g]==0.2*","~sigma[epsilon]==0.9,'  '*sigma[g]==0.8*","~sigma[epsilon]==0.9,' '*sigma[g]==1.2*","~sigma[epsilon]==0.9)

ghap <- ggplot(Rmse_analysis2,aes(x=ng,y=mu_RMSE))+
  geom_point(aes(colour=Model,shape = sg_se ),alpha=0.7)+geom_errorbar(aes(ymin = L1, ymax = L2,colour=Model,shape=sg_se ), width = 3)+
  geom_line(aes(colour=Model,shape=sg_se ),alpha=0.7)+scale_color_manual(values = custom_colors)+
  theme_bw()+xlab('Number of replicates') + ylab(expression(paste("Means of RMSE t ha"^"-1")))+ scale_x_continuous(breaks=c(10,30,50,seq(100, 500, 100)))+
  scale_shape_manual(values=c(15,16,17),labels=shapes,name = 'Variance components')

ghap

ggsave('Results/Replicates_se_sg2.png',ghap,height = 4  ,width = 8)

windows()

Rmse_analysis3 <- Rmse_analysis2 %>% mutate(ratio_sg_se=sigma_g/sigma_e)%>% 
  group_by(ng,model,ratio_sg_se) %>% summarise(mu_RMSE=mean(RMSE),L1 = mean(RMSE)-1.96*sd(RMSE)/30,
                                               L2 =  mean(RMSE)+1.96*sd(RMSE)/30) %>% filter(ng%in% c(30,100,500))

Rmse_analysis3 <- Rmse_analysis3 %>% rename(Model=model)

custom_colors <- c("Th" = "#FEA800",  # Deep Red
                   "PL" = "#00CCBE")  # Crimson

ghap <- ggplot(Rmse_analysis3,aes(x=ratio_sg_se,y=mu_RMSE))+
  geom_point(aes(colour=Model,shape = factor(ng) ),alpha=0.7)+geom_errorbar(aes(ymin = L1, ymax = L2,colour=Model,shape=factor(ng) ), width = 0.05)+
  geom_line(aes(colour=Model,shape=factor(ng )),alpha=0.7)+scale_color_manual(values = custom_colors)+
  theme_bw()+xlab('Ratio sigma g/sigma e') + ylim(c(0,5))
  ylab(expression(paste("Means of RMSE t ha"^"-1")))+ scale_x_continuous(breaks=unique(Rmse_analysis3$ratio_sg_se),labels  = round(unique(Rmse_analysis3$ratio_sg_se),2))#+
#scale_shape_manual(values=c(15,16,17),labels=shapes,name = 'Sample size')


ggsave('Results/ratio_sigma_e.png',ghap,height = 4  ,width = 12)

# adding graphic correlation person vs correlation sperman (maybe sperman isn't necessary)








