#######################################
### This script simulates some ST data
### then fits time series models.
### The entire script, as written, 
### takes 4 seconds to run on a  
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

### train and test model ###

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

### 21 step ahead forecast ###

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

one_step_table <- rbind(arima_stat,
                        exo_arima_stat,
                        auto_arima_stat)
rownames(one_step_table) <- c('fixed order ARIMA', 'ARIMA with exogenous variables', 'auto ARIMA')
knitr::kable(one_step_table)

toc <- proc.time()

cat(paste('it takes', round((toc-tic)[3]), 'seconds to finish on a 2015 Mac Book Pro with Intel Core i5 CPU.'))
