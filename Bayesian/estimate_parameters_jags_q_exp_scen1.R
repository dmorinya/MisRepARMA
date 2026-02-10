library(readxl)
library(dplyr)
library(MMWRweek)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(R2jags)
library(doParallel)

ncores <- detectCores() - 1  
registerDoParallel(cores=ncores)  
cl <- makeCluster(ncores, type="FORK")  

load("/home/dmorina/Insync/2102177@uab.cat/OneDrive Biz/Projectes/2022/0022022. DaSciVio/UnderreportedGBV/R/prop_survey.RData")
data_ts_Week <- read_xlsx("/home/dmorina/Insync/2102177@uab.cat/OneDrive Biz/Projectes/2022/0022022. DaSciVio/Data/GBV_TS.xlsx")

# JAGS to estimate parameters
jags_GBV <- function() {
  for (i in 1:N) {
  #  lambda0[i] <- ifelse(i < cp, lambda1, lambda2)
    thin[i] <- q0*step(cp-i)+(1-(1-q0)*exp(-alpha*(i-cp)))*step(i-cp) #exponential
    #thin[i] <- q0*step(cp-i)+(q0+(i-cp)/(1/(1-q0)*(alpha-cp)))*step(i-cp) #linear
    #z[i] ~ dpois(lambda0[i])
    #x[i] ~ dpois(lambda + Ind[i]*beta)
    y[i] ~ dpois(thin[i]*(lambda + Ind[i]*beta))
  }  
  # priors
  #lambda1 ~ dgamma(1, 1)
  #lambda2 ~ dgamma(3, 3)
  q0 ~ dbeta(1, 1)
  lambda ~ dnorm(lambda_mean, 100)
  beta ~ dgamma(2, 1)
  cp ~ dunif(0, N)
  #alpha ~ dunif(N, 5000) #linear
  alpha ~ dgamma(1, 3) #exp
}

#### BASE 25
lambda_mean <- prop_survey[names(prop_survey)=="BASE 25"]
datos <- list(y=data_ts_Week[data_ts_Week$BASE==25, ]$Cases, z=data_ts_Week[data_ts_Week$BASE==25, ]$Cases, 
              N=length(data_ts_Week[data_ts_Week$BASE==25, ]$Cases), Ind=data_ts_Week$I, lambda_mean=lambda_mean)
init_values <- function(){
  list(q0=runif(1), lambda=rnorm(1, lambda_mean, 0.1), beta=runif(1, 2, 1), 
       cp=runif(1, 0, length(data_ts_Week[data_ts_Week$BASE==25, ]$Cases)),
       #alpha=runif(1, length(data_ts_Week[data_ts_Week$BASE==25, ]$Cases), 5000), #linear
       alpha=rgamma(1, 1, 3), #exp
       .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
}
params <- c("q0", "lambda", "beta", "cp", "alpha")
clusterExport(cl, list("datos", "params", "init_values"))
fit_GBV_25_exp <- jags.parallel(data = datos, inits = NULL, parameters.to.save = params, model.file = jags_GBV,
                            n.chains = 5, n.iter = 50000, n.burnin = 5000, n.thin = 10, DIC = T)

#### BASE 26
lambda_mean <- prop_survey[names(prop_survey)=="BASE 26"]
datos <- list(y=data_ts_Week[data_ts_Week$BASE==26, ]$Cases, z=data_ts_Week[data_ts_Week$BASE==26, ]$Cases, 
              N=length(data_ts_Week[data_ts_Week$BASE==26, ]$Cases), Ind=data_ts_Week$I, lambda_mean=lambda_mean)
init_values <- function(){
  list(q0=runif(1), lambda=rnorm(1, lambda_mean, 0.1), beta=runif(1, 2, 1), 
       cp=runif(1, 0, length(data_ts_Week[data_ts_Week$BASE==26, ]$Cases)),
       #alpha=runif(1, length(data_ts_Week[data_ts_Week$BASE==26, ]$Cases), 5000), #linear
       alpha=rgamma(1, 1, 3), #exp
       .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
}
params <- c("q0", "lambda", "beta", "cp", "alpha")
clusterExport(cl, list("datos", "params", "init_values"))
fit_GBV_26_exp <- jags.parallel(data = datos, inits = NULL, parameters.to.save = params, model.file = jags_GBV,
                            n.chains = 5, n.iter = 50000, n.burnin = 5000, n.thin = 10, DIC = T)

#### BASE 27
lambda_mean <- prop_survey[names(prop_survey)=="BASE 27"]
datos <- list(y=data_ts_Week[data_ts_Week$BASE==27, ]$Cases, z=data_ts_Week[data_ts_Week$BASE==27, ]$Cases, 
              N=length(data_ts_Week[data_ts_Week$BASE==27, ]$Cases), Ind=data_ts_Week$I, lambda_mean=lambda_mean)
init_values <- function(){
  list(q0=runif(1), lambda=rnorm(1, lambda_mean, 0.1), beta=runif(1, 2, 1), 
       cp=runif(1, 0, length(data_ts_Week[data_ts_Week$BASE==27, ]$Cases)),
       #alpha=runif(1, length(data_ts_Week[data_ts_Week$BASE==27, ]$Cases), 5000), #linear
       alpha=rgamma(1, 1, 3), #exp
       .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
}
params <- c("q0", "lambda", "beta", "cp", "alpha")
clusterExport(cl, list("datos", "params", "init_values"))
fit_GBV_27_exp <- jags.parallel(data = datos, inits = NULL, parameters.to.save = params, model.file = jags_GBV,
                            n.chains = 5, n.iter = 50000, n.burnin = 5000, n.thin = 10, DIC = T)

#### BASE 31
lambda_mean <- prop_survey[names(prop_survey)=="BASE 31"]
datos <- list(y=data_ts_Week[data_ts_Week$BASE==31, ]$Cases, z=data_ts_Week[data_ts_Week$BASE==31, ]$Cases, 
              N=length(data_ts_Week[data_ts_Week$BASE==31, ]$Cases), Ind=data_ts_Week$I, lambda_mean=lambda_mean)
init_values <- function(){
  list(q0=runif(1), lambda=rnorm(1, lambda_mean, 0.1), beta=runif(1, 2, 1), 
       cp=runif(1, 0, length(data_ts_Week[data_ts_Week$BASE==31, ]$Cases)),
       #alpha=runif(1, length(data_ts_Week[data_ts_Week$BASE==31, ]$Cases), 5000), #linear
       alpha=rgamma(1, 1, 3), #exp
       .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
}
params <- c("q0", "lambda", "beta", "cp", "alpha")
clusterExport(cl, list("datos", "params", "init_values"))
fit_GBV_31_exp <- jags.parallel(data = datos, inits = NULL, parameters.to.save = params, model.file = jags_GBV,
                            n.chains = 5, n.iter = 50000, n.burnin = 5000, n.thin = 10, DIC = T)

#### BASE 34
lambda_mean <- 1 #### According to the macrosurvey this value is 0, which is not possible
datos <- list(y=data_ts_Week[data_ts_Week$BASE==34, ]$Cases, z=data_ts_Week[data_ts_Week$BASE==34, ]$Cases, 
              N=length(data_ts_Week[data_ts_Week$BASE==34, ]$Cases), Ind=data_ts_Week$I, lambda_mean=lambda_mean)
init_values <- function(){
  list(q0=runif(1), lambda=rnorm(1, lambda_mean, 0.1), beta=runif(1, 2, 1), 
       cp=runif(1, 0, length(data_ts_Week[data_ts_Week$BASE==34, ]$Cases)),
       #alpha=runif(1, length(data_ts_Week[data_ts_Week$BASE==34, ]$Cases), 5000), #linear
       alpha=rgamma(1, 1, 3), #exp
       .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
}
params <- c("q0", "lambda", "beta", "cp", "alpha")
clusterExport(cl, list("datos", "params", "init_values"))
fit_GBV_34_exp <- jags.parallel(data = datos, inits = NULL, parameters.to.save = params, model.file = jags_GBV,
                            n.chains = 5, n.iter = 50000, n.burnin = 5000, n.thin = 10, DIC = T)

#### BASE 35
lambda_mean <- prop_survey[names(prop_survey)=="BASE 35"]
datos <- list(y=data_ts_Week[data_ts_Week$BASE==35, ]$Cases, z=data_ts_Week[data_ts_Week$BASE==35, ]$Cases, 
              N=length(data_ts_Week[data_ts_Week$BASE==35, ]$Cases), Ind=data_ts_Week$I, lambda_mean=lambda_mean)
init_values <- function(){
  list(q0=runif(1), lambda=rnorm(1, lambda_mean, 0.1), beta=runif(1, 2, 1), 
       cp=runif(1, 0, length(data_ts_Week[data_ts_Week$BASE==35, ]$Cases)),
       #alpha=runif(1, length(data_ts_Week[data_ts_Week$BASE==35, ]$Cases), 5000), #linear
       alpha=rgamma(1, 1, 3), #exp
       .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
}
params <- c("q0", "lambda", "beta", "cp", "alpha")
clusterExport(cl, list("datos", "params", "init_values"))
fit_GBV_35_exp <- jags.parallel(data = datos, inits = NULL, parameters.to.save = params, model.file = jags_GBV,
                            n.chains = 5, n.iter = 50000, n.burnin = 5000, n.thin = 10, DIC = T)

#### GLOBAL NORTHERN METROPOLITAN BASIC HEALTH AREA
y <- data_ts_Week %>% group_by(Year, Week) %>% summarise(Cases=sum(Cases))
lambda_mean <- sum(prop_survey)
datos <- list(y=y$Cases, z=y$Cases, 
              N=length(y$Cases), Ind=data_ts_Week$I, lambda_mean=lambda_mean)
init_values <- function(){
  list(q0=runif(1), lambda=rnorm(1, lambda_mean, 0.1), beta=runif(1, 2, 1), 
       cp=runif(1, 0, length(y$Cases)),
       #alpha=runif(1, length(y$Cases), 5000), #linear
       alpha=rgamma(1, 1, 3), #exp
       .RNG.name = "base::Wichmann-Hill", .RNG.seed = 132023)
}
params <- c("q0", "lambda", "beta", "cp", "alpha")
clusterExport(cl, list("datos", "params", "init_values"))
fit_GBV_global_exp <- jags.parallel(data = datos, inits = NULL, parameters.to.save = params, model.file = jags_GBV,
                                n.chains = 5, n.iter = 200000, n.burnin = 20000, n.thin = 10, DIC = T)

save(list=c("fit_GBV_25_exp", "fit_GBV_26_exp", "fit_GBV_27_exp", "fit_GBV_31_exp", 
            "fit_GBV_34_exp", "fit_GBV_35_exp", "fit_GBV_global_exp"),
     file="Underreporting/R/parameters_exp.RData")
