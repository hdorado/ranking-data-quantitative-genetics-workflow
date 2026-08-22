
n.genotypes=10
mu.genotype =3.63
n.envs =1
sigma.genotype =0.8
sigma.env.farm=0
sigma.env=0
sigma.gen.env=0
sigma.plot=0.9
n.rep.genotype.env=100

require(PlackettLuce)
require(ClimMobTools)
require(lme4)
require(nlme)

source('Scripts/0. Ranking_pars_FUN.R')

setwd("Paper_1_Estimation_of_genetic_parameters")

set.seed(123)

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



data_simulated$trait.value.rank <- do.call(c,with(data_simulated,
                                                  tapply(-1*trait.value.obs, farm.vec,rank)))

# Output

pars = data.frame(n.rep.genotype.env=n.rep.genotype.env,n.packages=n.packages,n.genotypes=n.genotypes,mu.genotype=mu.genotype,n.envs =n.envs,
                  sigma.genotype=sigma.genotype,sigma.env.farm=sigma.env.farm,sigma.env=sigma.env,sigma.plot=sigma.plot,sigma.gen.env =sigma.gen.env)

ls_out <- list(dataset=data_simulated,rankings= convert_rankings(data_simulated = data_simulated,names.genotype=names.genotype),pars = pars)

ls_out$dataset$env.farm.vec <- factor(ls_out$dataset$env.farm.vec)





mod <- lmer(trait.value.obs~genotype+(1|env.farm.vec),data=ls_out$dataset)

summary(mod)

summary(mod)$coefficients

0.9097/sqrt(100)

emmeans::emmeans(mod,'genotype')

qvcalc(PlackettLuce(ls_out$rankings))




#-------------------------------------------------------------------------------

sigma.env.farm=1

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

# Using exactly the same variation of the first simulation

#coef.genotypes <- matrix(rnorm(n.genotypes,mean=mu.genotype,sd=sigma.genotype),ncol = 1 , nrow = n.genotypes)
#coef.envs      <- matrix(rnorm(n.envs,mean=0,sd=sigma.env) ,ncol = 1 , nrow =  n.envs)
#coef.gen.env   <- matrix(rnorm(n.genotypes*n.envs,mean=0,sd=sigma.gen.env),ncol = 1,nrow = n.genotypes*n.envs)
coef.env.farms <- matrix(rnorm(n.farms,mean=0,sd=sigma.env.farm) ,ncol = 1 , nrow =  n.farms) # farms only can be dependend to enviroments
#coef.plots     <- matrix(rnorm(n.plots*n.farms,mean=0,sd=sigma.plot),ncol = 1, nrow = n.plots*n.farms)


# Sum of effects

trait.value.obs <- genotype.m.design %*% coef.genotypes +  env.m.design %*% coef.envs+ gen.env.m.design %*% coef.gen.env + env.farm.m.design %*% coef.env.farms + farm.plot.vec.m.design %*% coef.plots

# complete effect

data_simulated <-
  data.frame(long.design,value.genotype = genotype.m.design %*% coef.genotypes,
             value.env = env.m.design %*% coef.envs,
             value.gen.env = gen.env.m.design %*% coef.gen.env,
             value.env.farm = env.farm.m.design %*% coef.env.farms,
             value.env.farm.plot = farm.plot.vec.m.design %*% coef.plots,trait.value.obs = trait.value.obs)



data_simulated$trait.value.rank <- do.call(c,with(data_simulated,
                                                  tapply(-1*trait.value.obs, farm.vec,rank)))

# Output

pars = data.frame(n.rep.genotype.env=n.rep.genotype.env,n.packages=n.packages,n.genotypes=n.genotypes,mu.genotype=mu.genotype,n.envs =n.envs,
                  sigma.genotype=sigma.genotype,sigma.env.farm=sigma.env.farm,sigma.env=sigma.env,sigma.plot=sigma.plot,sigma.gen.env =sigma.gen.env)

ls_out <- list(dataset=data_simulated,rankings= convert_rankings(data_simulated = data_simulated,names.genotype=names.genotype),pars = pars)

ls_out$dataset$env.farm.vec <- factor(ls_out$dataset$env.farm.vec)


mod <- lmer(trait.value.obs~genotype+(1|env.farm.vec),data=ls_out$dataset)

summary(mod)

summary(mod)$coefficients

emmeans::emmeans(mod,'genotype')

qvcalc(PlackettLuce(ls_out$rankings))

