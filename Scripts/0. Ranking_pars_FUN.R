
# Functions use ranking data analysis for plant breading
# Hugo Dorado
# creation: 10-06/2024

# Create a contrast matrix to be use within the Thurstonian model

contrast.diff.mat <- function(len){ # of items to compare
  
  if(len==2){return(t(matrix(c(1,-1))))}
  
  c0 <- matrix(c(1,rep(0,len-2)),ncol=1)
  
  dM <- cbind(c0,matrix(0,ncol = len-1,nrow = len-1,byrow=T))
  
  for(i in 2:(len-1)){
    cs <- dM[,i-1]
    cs[cs==-1] <-0 
    cs <- cs*-1
    cs[i] <-1
    dM[,i] <- cs
    
  }
  
  cf <- matrix(c(rep(0,len-2),-1),ncol=1)
  
  dM[,ncol(dM)] <- cf
  
  dM
}

# contrast.diff.mat(2);contrast.diff.mat(3)

# Loglikehood function of the thurstonian model

thurstonian_loglikehood <- function(alpha_bet){
  require(mvtnorm)
  ln <- length(alpha_bet)
  C <- contrast.diff.mat(ln)
  mu <- as.matrix(alpha_bet)
  
  Sm <- C%*%diag(length(alpha_bet))%*%t(C) # here the sigma is included, we are constraint to the identify matrix
  D <- if(nrow(Sm)==1){Sm^(-1/2)}else{diag(diag(Sm)^(-1/2))}
  mean <-  C %*% mu
  Smean <- D%*% C %*% mu
  Smf <- D%*%Sm%*%D
  pmv<-pmvnorm(lower=rep(-Inf,ln-1),upper=Smean[,1],mean=rep(0,ln-1),sigma= Smf) # First derivatives regard x and then y multivariate normal mean = c(0,0) Böckenholt 
  # pmvnorm(lower=rep(0,ln-1),upper=rep(Inf,ln-1),mean=mean[,1],sigma=Sm) # non-standarized form (this is equivalent)
  #1 comparison pnorm(-2/sqrt(2),mean = 0,sd = 1)
  log(pmv[1])
}

# thurstonian_loglikehood(alpha_bet=c(1,5,3))


# Placket Luce Loglikehood - only available for three items comparison

pl_likehood <- function(alpha_vec){
  log(alpha_vec[1])-log(alpha_vec[1]+alpha_vec[2]+alpha_vec[3])+
    log(alpha_vec[2])-log(alpha_vec[2]+alpha_vec[3])
}



# Loglikehood function for the sample (restriction mu = 0)

#par = matMu[-1]
#data= rankings

tricotNLLMU0 <- function(data,par,model= "t"){
  
  library(gtools)
  
  par <- c(mu1=0,par)
  
  # Process as a data.frame
  
  matrixTricot <- as.data.frame(as.matrix(data))
  
  nam.par <- names(par)
  
  colnames(matrixTricot) <- nam.par
  
  # Initial values
  
  muMat <- data.frame(nam.mu=nam.par,par=par)
  
  # Create worth values matrix
  
  matrixTricot <- cbind(Ind=row.names(matrixTricot),matrixTricot)
  
  matrixTricot <- matrixTricot %>% pivot_longer(cols = !Ind, names_to = "Item" ,values_to = "Ranking") 
  
  matrixTricot <- matrixTricot[mixedorder( matrixTricot$Ind),]
  
  matrixTricot <- matrixTricot %>% dplyr::filter(Ranking!=0)
  
  matrixTricot <- left_join(matrixTricot,muMat,by=c('Item'='nam.mu'))
  
  matrixTricot <- matrixTricot %>% dplyr::select(!Item) %>% pivot_wider(names_from =  Ranking,values_from = par)
  
  matrixTricot <- matrixTricot[,-1]
  
  matrixTricot <- matrixTricot[,mixedorder(colnames(matrixTricot))]
  
  matrixTricot <- as.matrix(matrixTricot)
  
  loglk <- if(model=='t'){apply(matrixTricot, 1, thurstonian_loglikehood)}else{logworth <- exp(matrixTricot); apply(logworth, 1, pl_likehood)}
  
  -1*sum(loglk,na.rm = T)
}

# likehood mean

#par = matMu2[-2]

tricotNLLMU0_ref <- function(data,par,model= "t"){
  
  library(gtools)
  
  par <- c(par['mu'],muG1=0,par[-which(names(par)=='mu')])
  
  # Process as a data.frame
  
  matrixTricot <- as.data.frame(as.matrix(data))
  
  nam.par <- names(par)
  
  colnames(matrixTricot) <- nam.par[-1]
  
  # Initial values
  
  muMat <- data.frame(nam.mu=nam.par,par=par)
  
  # Create worth values matrix
  
  matrixTricot <- cbind(Ind=row.names(matrixTricot),matrixTricot)
  
  matrixTricot <- matrixTricot %>% pivot_longer(cols = !Ind, names_to = "Item" ,values_to = "Ranking") 
  
  matrixTricot <- matrixTricot[mixedorder( matrixTricot$Ind),]
  
  matrixTricot <- matrixTricot %>% dplyr::filter(Ranking!=0)
  
  matrixTricot <- left_join(matrixTricot,muMat,by=c('Item'='nam.mu'))
  
  matrixTricot <- matrixTricot %>% dplyr::select(!Item) %>% pivot_wider(names_from =  Ranking,values_from = par)
  
  matrixTricot <- matrixTricot[,-1]
  
  matrixTricot <- matrixTricot[,mixedorder(colnames(matrixTricot))]
  
  matrixTricot <- as.matrix(matrixTricot) + par['mu']
  
  loglk <- if(model=='t'){apply(matrixTricot, 1, thurstonian_loglikehood)}else{logworth <- exp(matrixTricot); apply(logworth, 1, pl_likehood)}
  
  -1*sum(loglk,na.rm = T)
}

# the intecept is not succesfull because any value will canceled in difference of the parameters





# convert to ranking observed metric values from a tricot design.

convert_rankings <- function(data_simulated,names.genotype,env='env.vec',farm='farm.vec',genotype='genotype', trait.value.rank='trait.value.rank'){
  
  require(PlackettLuce)
  
  data_simulated <- do.call(rbind,split(data_simulated[c('env.vec','farm.vec','genotype','trait.value.rank')],data_simulated$farm.vec))
  
  
  data_simulated_wide <- stats::reshape(data_simulated,direction = 'wide',idvar=c("env.vec",'farm.vec'),timevar='genotype' ) 
  
  
  names(data_simulated_wide) <- gsub(paste0(trait.value.rank,"."),'',names(data_simulated_wide))
  
  
  data_simulated_wide <- data_simulated_wide[,c('env.vec','farm.vec',names.genotype)]
  
  
  data_simulated_wide[is.na(data_simulated_wide)] <- 0
  
  
  ranks <- as.rankings(data_simulated_wide[,-c(1:2)])
  
  rownames(ranks) <- paste0(data_simulated_wide$env.vec,".",data_simulated_wide$farm.vec)
  
  ranks
}


# Multi enviroment simulation

simulate_multi_env_tricot <- function(n.genotypes=27,mu.genotype =13.03,n.envs =4,
                                      sigma.genotype =0.91,sigma.env.farm=4.39,sigma.env=7.01,
                                      sigma.gen.env=1.11,sigma.plot=5.36,n.rep.genotype.env=16){
  
  require(PlackettLuce)
  require(ClimMobTools)
  
  n.packages <- round((n.genotypes*n.rep.genotype.env/3))
  
  names.genotype <- paste0('G',1:n.genotypes)
  
  # first random process
  
  geno.comb.matrix <- lapply(seq(n.envs),function(env){
    randomise(npackages=n.packages, itemnames= names.genotype)
  })
  
  sapply(geno.comb.matrix,nrow)
  geno.comb.matrix <- do.call(rbind,geno.comb.matrix)
  
  n.farms  <- nrow(geno.comb.matrix)
  n.plots  <- ncol(geno.comb.matrix)
  farm     <- paste0('F',1: n.farms ) # farms will be generated across the enviroments with not repetions
  plot     <- paste0('P',1:n.plots)
  env      <- paste0('Env',1:n.envs) # Note the optimal, enviroment must be stratified
  gen.env  <- paste(names.genotype,env,sep='.')
  
  ##set up data frame with assignment of varieties to farms, plots
  farm.vec <- rep(farm,each= n.plots)
  plot.vec <- rep(plot,length= length(farm.vec))
  env.vec  <- rep(env,array(n.packages*3,n.envs))
  
  
  # set up (small) deviations due to farm and plot (residual in this case)
  
  long.design <- data.frame(env.vec=env.vec,farm.vec = farm.vec, plot.vec = plot.vec , env.farm.vec = paste0(env.vec,'.',farm.vec), farm.plot.vec= paste0(env.vec,'.',farm.vec,'.',plot.vec), genotype = c(t(geno.comb.matrix)))
  
  long.design$gen.env.vec <- factor(paste0(long.design$genotype,'.',long.design$env.vec))
  
  long.design$gen.env.vec <- factor(long.design$gen.env.vec,levels = unique(long.design$gen.env.vec))
  
  long.design$env.vec <- factor(long.design$env.vec,levels = unique(env.vec))
  
  long.design$farm.vec <- factor(long.design$farm.vec,levels = unique(farm.vec)) 
  
  long.design$plot.vec <- factor(long.design$plot.vec,levels = unique(plot.vec))
  
  long.design$env.farm.vec <- factor(long.design$env.farm.vec,levels = unique(long.design$env.farm.vec))
  
  long.design$farm.plot.vec <- factor(long.design$farm.plot.vec,levels = unique(long.design$farm.plot.vec))
  
  long.design$genotype <- factor(long.design$genotype,levels = names.genotype)
  
  # Design matrix
  env.m.design <- if(n.envs!=1){model.matrix(~-1+env.vec,data=long.design)}else{matrix(array(0,nrow(long.design)),ncol=1)}
  env.farm.m.design      <- model.matrix(~-1+env.farm.vec,data=long.design)
  farm.plot.vec.m.design <- model.matrix(~-1+farm.plot.vec,data=long.design)
  genotype.m.design      <- model.matrix(~-1+genotype,data=long.design)
  gen.env.m.design       <- model.matrix(~-1+gen.env.vec,data=long.design)
  
  # Coefficients matrix (second random process)
  
  coef.genotypes <- matrix(rnorm(n.genotypes,mean=mu.genotype,sd=sigma.genotype),ncol = 1 , nrow = n.genotypes)
  coef.envs      <- matrix(rnorm(n.envs,mean=0,sd=sigma.env) ,ncol = 1 , nrow =  n.envs)
  coef.gen.env   <- matrix(rnorm(n.genotypes*n.envs,mean=0,sd=sigma.gen.env),ncol = 1,nrow = n.genotypes*n.envs)
  coef.env.farms <- matrix(rnorm(n.farms,mean=0,sd=sigma.env.farm) ,ncol = 1 , nrow =  n.farms) # farms only can be dependend to enviroments
  coef.plots     <- matrix(rnorm(n.plots*n.farms,mean=0,sd=sigma.plot),ncol = 1, nrow = n.plots*n.farms)
  
  
  # Sum of effects
  
  trait.value.obs <- genotype.m.design %*% coef.genotypes +  env.m.design %*% coef.envs+ gen.env.m.design %*% coef.gen.env + env.farm.m.design %*% coef.env.farms + farm.plot.vec.m.design %*% coef.plots
  
  # complete effect
  
  data_simulated <-
    data.frame(long.design,value.genotype = genotype.m.design %*% coef.genotypes,
               value.env = env.m.design %*% coef.envs,
               value.gen.env = gen.env.m.design %*% coef.gen.env,
               value.env.farm = env.farm.m.design %*% coef.env.farms,
               value.env.farm.plot = farm.plot.vec.m.design %*% coef.plots,trait.value.obs = trait.value.obs)
  
  
  # This line of code is valid as long the line long.design$farm.vec <- factor(long.design$farm.vec,levels = unique(farm.vec)) is writen, other its going to be a mess
  
  data_simulated$trait.value.rank <- do.call(c,with(data_simulated,
                                                    tapply(-1*trait.value.obs, farm.vec,rank)))
  
  
  #temds <- data_simulated %>% group_by(farm.vec) %>% mutate(rank2=rank(-trait.value.obs))
  
  #temds$trait.value.rank == temds$rank2
  
  
  # Output
  
  pars = data.frame(n.rep.genotype.env=n.rep.genotype.env,n.packages=n.packages,n.genotypes=n.genotypes,mu.genotype=mu.genotype,n.envs =n.envs,
                    sigma.genotype=sigma.genotype,sigma.env.farm=sigma.env.farm,sigma.env=sigma.env,sigma.plot=sigma.plot,sigma.gen.env =sigma.gen.env)
  
  list(dataset=data_simulated,rankings= convert_rankings(data_simulated = data_simulated,names.genotype=names.genotype),pars = pars)
}


# Example
# 
# a <- 
# simulate_multi_env_tricot(n.genotypes=10,mu.genotype =13.03,n.envs =4,
#                           sigma.genotype =0.91,sigma.env.farm=4.39,
#                           sigma.env=7.01,sigma.gen.env=1.11,sigma.plot=5.36,
#                           n.rep.genotype.env=10)
# 
# a$
# 
# a$rankings
# 
# a$dataset



tricotNLLMU0_summZ <- function(data,par,model= "t"){
  
  library(gtools)
  
  par <- c(par['mu'],muG1=0,par[-which(names(par)=='mu')])
  
  # Process as a data.frame
  
  matrixTricot <- as.data.frame(as.matrix(data))
  
  nam.par <- names(par)
  
  colnames(matrixTricot) <- nam.par[-1]
  
  # Initial values
  
  muMat <- data.frame(nam.mu=nam.par,par=par)
  
  # Create worth values matrix
  
  matrixTricot <- cbind(Ind=row.names(matrixTricot),matrixTricot)
  
  matrixTricot <- matrixTricot %>% pivot_longer(cols = !Ind, names_to = "Item" ,values_to = "Ranking") 
  
  matrixTricot <- matrixTricot[mixedorder( matrixTricot$Ind),]
  
  matrixTricot <- matrixTricot %>% dplyr::filter(Ranking!=0)
  
  matrixTricot <- left_join(matrixTricot,muMat,by=c('Item'='nam.mu'))
  
  matrixTricot <- matrixTricot %>% dplyr::select(!Item) %>% pivot_wider(names_from =  Ranking,values_from = par)
  
  matrixTricot <- matrixTricot[,-1]
  
  matrixTricot <- matrixTricot[,mixedorder(colnames(matrixTricot))]
  
  matrixTricot <- as.matrix(matrixTricot) + par['mu']
  
  loglk <- if(model=='t'){apply(matrixTricot, 1, thurstonian_loglikehood)}else{logworth <- exp(matrixTricot); apply(logworth, 1, pl_likehood)}
  
  -1*sum(loglk,na.rm = T)
}

