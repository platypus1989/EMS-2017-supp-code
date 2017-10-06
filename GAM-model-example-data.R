#######################################
### This script simulates some ST data
### then fits a GAM and forecasts.
###
### This code was used for GAM modeling
### in the EMS paper, "Forecasting Urban 
### Water Demand with Statistical and
### Machine Learning Methods Using Large 
### Space-Time Data" 
### (Duerr & Merrill et. al. 2017+)
###
### Hunter R. Merrill
### 2017.09.07
#######################################

rm(list=ls()) #clear data
graphics.off(); cat('\014') #clear figures and console
set.seed(1) #set random seed
library(fields) #for distance matrix
library(Rcpp); library(RcppArmadillo) #for C++ stuff

### Simulate some data ###

ns = 25; nt = 100; nn = ns*nt #25 locations, 100 time points

X1 = rnorm(nn); X2 = rnorm(nn); X3 = rnorm(nn); X4 = rnorm(nn) #four exogenous predictors
locs = cbind(runif(ns), runif(ns)) #locations on [0,1]^2
times = 1:nt #100 time points, evenly spaced (this code/model works even if not evenly spaced)

phi = c(1,.25) #one corr parameter for space, one for time
sigma2 = .25 #nugget
tau2 = 1 #ST variation scale parameter

dist_S = as.matrix(dist(locs))
dist_T = as.matrix(dist(times)) #distance matrices

my.dat = data.frame('X1' = X1,
                    'X2' = X2,
                    'X3' = X3,
                    'X4' = X4,
                    'times' = rep(times,each=ns),
                    'locs' = locs) #create data frame

my.dat = my.dat[order(my.dat$locs.1, my.dat$times),] #order for tensor product representation

my.dat$Y = 10 + 1*sin(my.dat$X1) + .5*my.dat$X2 + 1*my.dat$X3*my.dat$X4 + #linear predictor part...
  t(chol(tau2*exp(-phi[1]*dist_S) %x% exp(-phi[2]*dist_T)))%*%rnorm(nn) + #...ST part...
  rnorm(nn, sd = sqrt(sigma2)) #...and white noise.

### train and test model ###

which.test = which(my.dat$times >= 80) #we will forecast time points 80 - 100.
which.train = which(my.dat$times < 80)

dat.train = my.dat[which.train,]
dat.test = my.dat[which.test,]

#fit GAM, get predictions
library(mgcv)
fit.gam = gam(Y ~ s(X1) + s(X2) + te(X3, X4) + 
                te(locs.1, locs.2, times, d = c(2,1), k = c(10,10), bs = c('tp','tp')), 
              data = dat.train)
plot(fit.gam, pages=1)
pred.gam = predict(fit.gam, newdata = dat.test)

par(mfrow = c(1,1), mar=c(3,3,1,1), mgp = c(2,1,0))
plot(dat.test$Y, pred.gam, xlab = 'Observed', ylab = 'Forecasted - GAM'); abline(0,1)
legend('bottomright', pch = NA, legend = paste('R2 =', round(cor(dat.test$Y, pred.gam),2)))
