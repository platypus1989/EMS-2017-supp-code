#######################################
### This script simulates some ST data
### then fits time series models.
### The entire script, as written, 
### takes XXX seconds to run on a  
### modest desktop
###
### This code replicates that used for time series
### models in the EMS paper, "Forecasting Urban 
### Water Demand with Statistical and
### Machine Learning Methods Using Large 
### Space-Time Data" 
### (Duerr & Merrill et. al. 2017+)
###
### Isaac Duerr
### 2017.10.05
#######################################

rm(list=ls()) #clear data
sapply(dev.list(), dev.off); cat('\014') #clear figures and console
set.seed(1) #set random seed
library(randomForest)
library(quantregForest)
library(caret)
library(gbm)
library(BayesTree)


#Load functions for Gini and NOIS calculations
source("/Users/ufgi/Dropbox (UFL)/uf_FFL_irrigation/gini.r")
#source("gini.r")
source("20160720_NOIS.R")

### Simulate some data, we will use the same data as the other examples ###


##Spatial/temporal data

ns = 25; nt = 100; nn = ns*nt #25 locations, 100 time points

X1 = rnorm(nn); X2 = rnorm(nn) #two exogenous predictors
locs = cbind(runif(ns), runif(ns)) #locations on [0,1]^2
locID <- 1:25
times = 1:nt

phi = c(1,.25) #one corr parameter for space, one for time
sigma2 = .25 #nugget
tau2 = 1 #ST variation scale parameter

dist_S = as.matrix(dist(locs))
dist_T = as.matrix(dist(times))

my.dat = data.frame('X1' = X1,
                    'X2' = X2,
                    'times' = rep(times,each=ns),
                    'locs' = locs,
                    "ID" = locID)

t_s = 10 + 1*my.dat$X1 + .5*my.dat$X2 + #linear predictor part...
  t(chol(tau2*exp(-phi[1]*dist_S) %x% exp(-phi[2]*dist_T)))%*%rnorm(nn) + #...ST part...
  rnorm(nn, sd = sqrt(sigma2)) #...and white noise.

which.test = which(my.dat$times >= 80) #we will forecast time points 80 - 100.
which.train = which(my.dat$times < 80)

dat.train = t_s[which.train,]
dat.test = t_s[which.test,]

num_loc = length(unique(my.dat$ID))

 

indices <- unlist(lapply(1:num_loc, function(x) (x-1)*80+1))
all.loc <- insert(dat.train, at=indices, values=list(rep(NA,100)))

fit.ar1.all <- arima(all.loc, order=c(1,0,0), method="ML")
fixed.phi <- as.numeric(fit.ar1.all$coef[1])
fixed.intercept <- as.numeric(fit.ar1.all$coef[2])

ar1.fixed.all.predictions <- rep(NA, num_loc*20)
ar1.random.intercept.predictions <- rep(NA, num_loc*20)
ar1.random.intercept.upper <- rep(NA, num_loc*20)
ar1.global.intercept.upper <- rep(NA, num_loc*20)

for(i in 1:num_loc){
  train.start <- 79 * (i-1) + 1
  train.end <- 79 * i
  test.start <- 21 * (i-1) + 1
  test.end <- 21 * i
  
  # Predict for each house using the global model with fixed intercept
  fit.fixed.new <- arima(dat.train[train.start:train.end], 
                         order=c(1,0,0), fixed=c(fixed.phi, fixed.intercept), transform.pars=FALSE)
  predict.ar1.fixed <- predict(fit.fixed.new, n.ahead=21, se.fit=T)
  predict.ar1.fixed.upper <- predict.ar1.fixed$pred + 1.645 * predict.ar1.fixed$se
  ar1.fixed.all.predictions[test.start:test.end] <- as.numeric(predict.ar1.fixed$pred)
  ar1.global.intercept.upper[test.start:test.end] <- as.numeric(predict.ar1.fixed.upper)
  # Predict for each house using the model with fixed phi and random intercept
  fit.random.intercept.new <- arima(dat.train[train.start:train.end], 
                                    order=c(1,0,0), fixed=c(fixed.phi, NA), transform.pars=FALSE)
  predict.ar1.random.intercept <- predict(fit.random.intercept.new, n.ahead=21, se.fit=T)
  predict.ar1.random.intercept.upper <- predict.ar1.random.intercept$pred + 1.645 * predict.ar1.random.intercept$se
  ar1.random.intercept.predictions[test.start:test.end] <- as.numeric(predict.ar1.random.intercept$pred)
  ar1.random.intercept.upper[test.start:test.end] <- as.numeric(predict.ar1.random.intercept.upper)
  
}

# Find overall MSEs
ar1.mses.global <- sqrt( mean( (ar1.fixed.all.predictions - dat.test)^2) ) 
ar1.gini.global <- NormalizedGini(dat.test,ar1.fixed.all.predictions)
ar1.ECPI.global <- mean(dat.test< ar1.global.intercept.upper)
ar1.AWPI.global <- mean(ar1.global.intercept.upper)
ar1.NOIS.global <- NOIS(dat.test,cbind(NA,ar1.global.intercept.upper))  

ar1.mses.random <- sqrt( mean( (ar1.random.intercept.predictions 
                                    - dat.test)^2) ) 
ar1.gini.random <- NormalizedGini(dat.test,ar1.random.intercept.predictions)
ar1.ECPI.random <- mean(dat.test< ar1.random.intercept.upper)
ar1.AWPI.random <- mean(ar1.random.intercept.upper)
ar1.NOIS.random <- NOIS(dat.test,cbind(NA,ar1.random.intercept.upper))  


