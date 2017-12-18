#######################################
### This script simulates some ST data
### then fits a model and forecasts.
### The entire script, as written, 
### takes 140 seconds to run on a  
### modest laptop.
###
### This code was used for ST modeling
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

ns = 50; nt = 50; nn = ns*nt #50 locations, 50 time points

X1 = rnorm(nn); X2 = rnorm(nn); X3 = rnorm(nn); X4 = rnorm(nn) #four exogenous predictors
locs = cbind(runif(ns), runif(ns)) #locations on [0,1]^2
times = 1:nt #100 time points, evenly spaced (this code/model works even if not evenly spaced)

phi = c(.25,.15) #one corr parameter for space, one for time
sigma2 = .25 #nugget
tau2 = 2.5 #ST variation scale parameter

dist_S = as.matrix(dist(locs))
dist_T = as.matrix(dist(times)) #distance matrices

my.dat = data.frame('X1' = X1,
                    'X2' = X2,
                    'X3' = X3,
                    'X4' = X4,
                    'times' = rep(times,each=ns),
                    'locs' = locs) #create data frame

my.dat = my.dat[order(my.dat$locs.1, my.dat$times),] #order for tensor product representation

my.dat$Y = 10 + 1*my.dat$X1 + .5*my.dat$X2 + 1*my.dat$X3*my.dat$X4 + #linear predictor part...
  t(chol(tau2*exp(-phi[1]*dist_S) %x% exp(-phi[2]*dist_T)))%*%rnorm(nn) + #...ST part...
  rnorm(nn, sd = sqrt(sigma2)) #...and white noise.

### train and test model ###

which.test = which(my.dat$times >= 40) #we will forecast time points 40 - 50.
which.train = which(my.dat$times < 40)

dat.train = my.dat[which.train,]
dat.test = my.dat[which.test,]

dist_T.new = as.matrix(dist(unique(dat.test$times)))
dist_T.on = as.matrix(rdist(unique(dat.train$times),unique(dat.test$times)))
dist_T.old = as.matrix(dist(unique(dat.train$times))) #correlation matrices, for prediction

sourceCpp('ST-model.cpp') #load C++ code

Xlin = model.matrix(as.formula(paste("~ (", 
                                     paste(names(dat.train)[grep('X',names(dat.train))], 
                                           collapse = " + "), 
                                     ")^2")), data = dat.train) #make design matrix

n.samples = 2500 #number of MCMC samples
print_step = 100 #print progress every 100 samples
tune = TRUE #should MH proposals be tuned?
tune_after = 200 #when should they start? (let samples get somewhere first)
tune_every = 100 #how often should they be tuned? (how often to calculate previous acceptance rate)
tune_for = 800 #how long should they be tuned? (stop tuning after tune_for + tune_after)
burn.in = tune_after + tune_for #throw these out

ST_samples = ST_MCMC(X_lin = Xlin, #linear model matrix
                    Y = dat.train$Y, #response data
                    dist_S = dist_S, #spatial distance matrix
                    dist_T = dist_T.old, #temporal distance matrix
                    n_samples = n.samples, #number of samples
                    print_step = print_step, #print every 100 samples
                    tune = tune, #should MH proposals be tuned?
                    tune_after = tune_after, #when should they start? (let samples get somewhere first)
                    tune_every = tune_every, #how often should they be tuned? (how often to calculate previous acceptance rate)
                    tune_for = tune_for) #how long should they be tuned? (stop tuning after tune_for + tune_after)

#check estimates:
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(2,1,0))
betas = c(10,1,.5,0,0,0,0,0,0,0,0) #true values
for(i in 1:11) {
  plot(ST_samples$b.lin.samples[i,burn.in:n.samples], type = 'l',
       ylab = bquote(beta[.(i)]),
       ylim = range(c(ST_samples$b.lin.samples[i,burn.in:n.samples], betas[i])))
  abline(h = betas[i], col = 'blue',lwd=2)
}

par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(2,1,0))
plot(ST_samples$tau2[burn.in:n.samples],  type = 'l',
     ylab = bquote(tau^2),
     ylim = range(c(ST_samples$tau2[burn.in:n.samples], tau2)))
abline(h = tau2, col = 'blue',lwd=2)
plot(ST_samples$sig2[burn.in:n.samples],  type = 'l',
     ylab = bquote(sigma^2),
     ylim = range(c(ST_samples$sig2[burn.in:n.samples], sigma2)))
abline(h = sigma2, col = 'blue',lwd=2)
plot(ST_samples$phi[1,burn.in:n.samples],  type = 'l',
     ylab = bquote(phi['S']),
     ylim = range(c(ST_samples$phi[1,burn.in:n.samples], phi[1])))
abline(h = phi[1], col = 'blue',lwd=2)
plot(ST_samples$phi[2,burn.in:n.samples],  type = 'l',
     ylab = bquote(phi['T']),
     ylim = range(c(ST_samples$phi[2,burn.in:n.samples], phi[2])))
abline(h = phi[2], col = 'blue',lwd=2)

#make forecasts with Bayesian BLUP:
thin = 10
Xnew = model.matrix(as.formula(paste("~ (", 
                                     paste(names(dat.train)[grep('X',names(dat.train))], 
                                           collapse = " + "), 
                                     ")^2")), data = dat.test)

Y.pred = matrix(0, nrow = nrow(dat.test), ncol = (n.samples-burn.in)/thin)
count = 0
for(i in (burn.in + 1):n.samples){
  if(i%%thin == 0){
    count = count + 1
    Y.pred[,count] = predict_new(mn_new = Xnew%*%ST_samples$b.lin.samples[,i], 
                                 mn_old = Xlin%*%ST_samples$b.lin.samples[,i],
                                 Y_old = dat.train$Y,
                                 Sigma_S_old = exp(-ST_samples$phi.samples[1,i]*dist_S),
                                 Sigma_T_old = exp(-ST_samples$phi.samples[2,i]*dist_T.old),
                                 Sigma_S_new = exp(-ST_samples$phi.samples[1,i]*dist_S),
                                 Sigma_T_new = exp(-ST_samples$phi.samples[2,i]*dist_T.new),
                                 Sigma_S_on = exp(-ST_samples$phi.samples[1,i]*dist_S),
                                 Sigma_T_on = exp(-ST_samples$phi.samples[2,i]*dist_T.on),
                                 sigma2 = ST_samples$sig2.samples[i],
                                 tau2 = ST_samples$tau2.samples[i])
  }
}

Y.pred.mns = matrix(nrow=ns, ncol=length(unique(dat.test$times)))
Y.obs.vals = matrix(nrow=ns, ncol=length(unique(dat.test$times)))
for(i in 1:length(unique(dat.test$times))){
  Y.pred.mns[,i] = apply(Y.pred[1:ns*length(unique(dat.test$times)) - 
                                  length(unique(dat.test$times)) + i,], 1, mean)
  Y.obs.vals[,i] = dat.test$Y[dat.test$times == i + 39]
}

# better than lm?
pred.lm = predict(lm(Y ~ (X1 + X2 + X3 + X4)^2, data = dat.train), newdata = dat.test)

par(mfrow = c(2,1), mar=c(3,3,1,1), mgp = c(2,1,0))
plot(c(Y.obs.vals), c(Y.pred.mns), xlab = 'Observed', ylab = 'Forecasted - ST'); abline(0,1)
legend('bottomright', pch = NA, legend = paste('R2 =', round(cor(c(Y.obs.vals), c(Y.pred.mns)),2)))
plot(dat.test$Y, pred.lm, xlab = 'Observed', ylab = 'Forecasted - LM'); abline(0,1)
legend('bottomright', pch = NA, legend = paste('R2 =', round(cor(dat.test$Y, pred.lm),2)))
# yes.
