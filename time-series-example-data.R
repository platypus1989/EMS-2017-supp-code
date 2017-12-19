#######################################
### This script simulates some ST data
### then fits time series models.
### The entire script, as written, 
### takes 5 seconds to run on a  
### modest laptop
###
### This code replicates that used for time series
### models in the EMS paper, "Forecasting Urban 
### Water Demand with Statistical and
### Machine Learning Methods Using Large 
### Space-Time Data" 
### (Duerr & Merrill et. al. 2017+)
###
### Chuan Wang
### 2017.12.18
#######################################

rm(list=ls()) #clear data
sapply(dev.list(), dev.off); cat('\014') #clear figures and console
tic <- proc.time()
set.seed(1) #set random seed
library(forecast)
library(R.utils)
#Load functions for Gini and NOIS calculations
source("gini.R")
source("NOIS.R")

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

my.dat$Y = as.vector(10 + 1*sin(my.dat$X1) + .5*my.dat$X2 + #linear predictor part...
                       t(chol(tau2*exp(-phi[1]*dist_S) %x% exp(-phi[2]*dist_T)))%*%rnorm(nn) + #...ST part...
                       rnorm(nn, sd = sqrt(sigma2))) #...and white noise.

which.test = which(my.dat$times >= 80) #we will forecast time points 80 - 100.
which.train = which(my.dat$times < 80)

ts.train = matrix(my.dat[which.train,'Y'], ncol = ns, byrow=TRUE)
ts.test = matrix(my.dat[which.test,'Y'], ncol = ns, byrow=TRUE)

X1.train = matrix(my.dat[which.train,'X1'], ncol = ns, byrow=TRUE)
X1.test = matrix(my.dat[which.test,'X1'], ncol = ns, byrow=TRUE)

X2.train = matrix(my.dat[which.train,'X2'], ncol = ns, byrow=TRUE)
X2.test = matrix(my.dat[which.test,'X2'], ncol = ns, byrow=TRUE)


### Evaluation Metrics ###
Eval_metric <- function(x, y, int){
  return(c(
    RMSE = sqrt(mean((x-y)^2)),
    gini = NormalizedGini(as.vector(t(x)),as.vector(t(y))),
    NOIS = NOIS(as.vector(t(x)),cbind(as.vector(t(int[,,1])),as.vector(t(int[,,2])))),
    AWPI = mean(int[,,2] - int[,,1]),
    ECPI = mean((x < int[,,2]) & x > int[,,1])
  ))
}

### AR1 Models ###
t_s = 10 + 1*my.dat$X1 + .5*my.dat$X2 + #linear predictor part...
  t(chol(tau2*exp(-phi[1]*dist_S) %x% exp(-phi[2]*dist_T)))%*%rnorm(nn) + #...ST part...
  rnorm(nn, sd = sqrt(sigma2)) #...and white noise.

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
ar1_global_stat <- c(RMSE = sqrt( mean( (ar1.fixed.all.predictions - dat.test)^2) ),
                     gini = NormalizedGini(dat.test,ar1.fixed.all.predictions),
                     NOIS = NOIS(dat.test,cbind(NA,ar1.global.intercept.upper)),
                     AWPI = mean(ar1.global.intercept.upper),
                     ECPI = mean(dat.test< ar1.global.intercept.upper))

ar1_random_stat <- c(RMSE = sqrt( mean( (ar1.random.intercept.predictions - dat.test)^2) ),
                     gini = NormalizedGini(dat.test,ar1.random.intercept.predictions),
                     NOIS = NOIS(dat.test,cbind(NA,ar1.random.intercept.upper)),
                     AWPI = mean(ar1.random.intercept.upper),
                     ECPI = mean(dat.test< ar1.random.intercept.upper))


### ARIMA Models ###

# ARIMA with fixed order (1,1,1)
pred <- matrix(NA, nrow = nrow(ts.test), ncol = ncol(ts.test))
pred_int <- array(NA, dim=c(nrow(ts.test), ncol(ts.test),2))
for (i in 1:ns){
  arima_model <- Arima(ts.train[,i], order = c(1,1,1), method="CSS")
  temp_forecast <- forecast(arima_model, h=21)
  pred[,i] <- temp_forecast$mean
  pred_int[,i,1] <- temp_forecast$lower[,2]
  pred_int[,i,2] <- temp_forecast$upper[,2]
}
arima_stat <- Eval_metric(ts.test, pred, pred_int)

# ARIMA with exogenous variables
pred <- matrix(NA, nrow = nrow(ts.test), ncol = ncol(ts.test))
pred_int <- array(NA, dim=c(nrow(ts.test), ncol(ts.test),2))
for (i in 1:ns){
  train_xreg <- cbind(X1.train[,i], X2.train[,i])
  arima_model <- Arima(ts.train[,i], order = c(1,1,1),
                       xreg=train_xreg, method="CSS")
  test_xreg <- cbind(X1.test[,i], X2.test[,i])
  temp_forecast <- forecast(arima_model, h=21,
                            xreg=test_xreg)
  pred[,i] <- temp_forecast$mean
  pred_int[,i,1] <- temp_forecast$lower[,2]
  pred_int[,i,2] <- temp_forecast$upper[,2] 
}
exo_arima_stat <- Eval_metric(ts.test, pred, pred_int)

# auto ARIMA
pred <- matrix(NA, nrow = nrow(ts.test), ncol = ncol(ts.test))
pred_int <- array(NA, dim=c(nrow(ts.test), ncol(ts.test),2))
for (i in 1:ns){
  arima_model <- auto.arima(ts.train[,i])
  temp_forecast <- forecast(arima_model, h=21)
  pred[,i] <- temp_forecast$mean
  pred_int[,i,1] <- temp_forecast$lower[,2]
  pred_int[,i,2] <- temp_forecast$upper[,2] 
}
auto_arima_stat <- Eval_metric(ts.test, pred, pred_int)


### Summary ###
summary_table <- rbind(ar1_global_stat,
                       ar1_random_stat,
                       arima_stat,
                       exo_arima_stat,
                       auto_arima_stat)
rownames(summary_table) <- c('AR1 with fixed intercept',
                             'AR1 with random intercept',
                             'fixed order ARIMA', 
                             'ARIMA with exogenous variables', 
                             'auto ARIMA')
knitr::kable(summary_table)

toc <- proc.time()

cat(paste('it takes', round((toc-tic)[3]), 'seconds to finish on a 2015 Mac Book Pro with Intel Core i5 CPU.'))
