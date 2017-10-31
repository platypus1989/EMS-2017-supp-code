#######################################
### This script simulates some data
### then fits machine learning models
### The entire script, as written, 
### takes 120 seconds to run on a  
### modest desktop
###
### This code replicates that used for machine learning
### in the EMS paper, "Forecasting Urban 
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
library(truncnorm)


#Load functions for Gini and NOIS calculations
source("gini.R")
#source("gini.r")
source("NOIS.R")

### Simulate some data ###


##Spatial/temporal data

ns = 25; nt = 100; nn = ns*nt #25 locations, 100 time points

X1 = rnorm(nn); X2 = rnorm(nn) #two exogenous predictors
locs = cbind(runif(ns), runif(ns)) #locations on [0,1]^2
times = 1:nt

phi = c(1,.25) #one corr parameter for space, one for time
sigma2 = .25 #nugget
tau2 = 1 #ST variation scale parameter

dist_S = as.matrix(dist(locs))
dist_T = as.matrix(dist(times))

tmp.dat = data.frame('X1' = X1,
                     'X2' = X2,
                     'times' = rep(times,each=ns),
                     'locs' = locs)

tmp.dat$Y = 10 + 1*tmp.dat$X1 + .5*tmp.dat$X2 + #linear predictor part...
  t(chol(tau2*exp(-phi[1]*dist_S) %x% exp(-phi[2]*dist_T)))%*%rnorm(nn) + #...ST part...
  rnorm(nn, sd = sqrt(sigma2)) #...and white noise.


sim_dat <- data.frame(matrix(0,nt*ns,11))
colnames(sim_dat) <- c("RAW1_mn","P_mn","ET","HeatArea","LndVal","Days","AW_percent_month_mean","FPE","GreenSpace","Bldval","Act_max")

#Covariates
load("params.RData")

for(var in colnames(sim_dat)){
  sim_dat[,var] <- rtruncnorm(ns*nt, a=0, mean = (params_lm[params_lm$vars == var, "int"] + params_lm[params_lm$vars == var, "mean"]*10*tmp.dat$Y), sd = params_lm[params_lm$vars == var,"sd"]/2)
}


#simulated data data
my.dat = data.frame("X1" = tmp.dat$X1,
                         "X2" = tmp.dat$X2,
                         "times" = tmp.dat$times,
                         "RAW1_mn" = sim_dat$RAW1_mn,
                         "P_mn" = sim_dat$P_mn,
                         "ET" = sim_dat$ET,
                         "HeatArea" = sim_dat$HeatArea,
                         "LndVal" = sim_dat$LndVal,
                         "AW_percent_month_mean" = sim_dat$AW_percent_month_mean,
                         "FPE" = sim_dat$FPE,
                         "GreenSpace" = sim_dat$GreenSpace,
                         "Bldval" = sim_dat$Bldval,
                         "Act_max" = sim_dat$Act_max,
                         "TotWater" = tmp.dat$Y)



### train and test model ###

which.test = which(my.dat$times >= 80) #we will forecast time points 80 - 100.
which.train = which(my.dat$times < 80)

dat.train = my.dat[which.train,]
dat.test = my.dat[which.test,]



##Crossvalidation set up for a selected model (random forest)##
control <- trainControl(method="repeatedcv", number=2, repeats=1)

rf_tuning <- data.frame(mtry = seq(2,13, by=2))
#Train and Tune the RF
rf.mod <- train(x=dat.train[,-which(colnames(my.dat)=="TotWater")],
                  y= dat.train$TotWater,
                  method = "rf",   # Randomforest
                  metric="RMSE",
                  verbose = FALSE,
                  trControl=control,
                  importance = T)


#Tuning parameters from best performing model
best_mtry <- rf.mod$finalModel$mtry
best_ntrees <- rf.mod$finalModel$ntree

#Refit for confidence intervals
best_model <- quantregForest(dat.train[,-which(colnames(my.dat)=="TotWater")], dat.train$TotWater,mtry = best_mtry, ntree=best_ntrees, importance=TRUE)

#variable importance measure from best performing model, in this set up the spatial components are most important
VI_RF <- varImp(rf.mod, scale=T)
plot(VI_RF)

##predict future timepoints
predicts<-predict(best_model, dat.test, what = c(0.05,.5,.95))

#Calculate performance metrics
rMSE<-sqrt(sum((dat.test$TotWater - predicts[,2])^2)/dim(dat.test)[1])
gini<-NormalizedGini(dat.test$TotWater,predicts[,2])
PI_width <- abs(predicts[,1] - predicts[,3])
nois <- NOIS(dat.test$TotWater,cbind(predicts[,1],predicts[,3]))
AWPI <- mean(PI_width)
ECPI <- sum(predicts[,1] <= dat.test$TotWater & dat.test$TotWater < predicts[,3])/length(dat.test$TotWater)
results <- data.frame(rMSE=rMSE, gini=gini, NOIS = nois, AWPI = AWPI, ECPI = ECPI)

results



##Fitting BART and GBM, for more details on tuning and calculating confidence intervals for each model please refer to An Introduction to Statistical Learning
#with Applications in R by Gareth James, Daniela Witten, Trevor Hastie and Robert Tibshirani. Coding tutorials are available at http://www-bcf.usc.edu/~gareth/ISL/


#BART
bart.mod <- bart(dat.train[,-which(colnames(my.dat)=="TotWater")],dat.train$TotWater, x.test = dat.test[,-which(colnames(my.dat)=="TotWater")], ntree=100)
predicts <- bart.mod$yhat.test.mean
(RMSE<-sqrt(sum((dat.test$TotWater-predicts)^2)/dim(dat.test)[1]))
(gini<-NormalizedGini(dat.test$TotWater,predicts))



#GBM
GBM_form <- paste("TotWater ~ ",names(dat.train[,-which(colnames(my.dat)=="TotWater")])[1])
for(name in names(dat.train[,-which(colnames(my.dat)=="TotWater")])[-1]){
  GBM_form <- paste(GBM_form,"+",name)
}
GBM_form <- formula(GBM_form)

GBM.fit <- gbm( form = GBM_form,
                data = dat.train,
                distribution = "gaussian",  
                var.monotone = NULL,
                n.trees = 200,
                interaction.depth = 2,  
                n.minobsinnode = 100,
                shrinkage = 0.01,
                bag.fraction = 0.5,
                train.fraction = 0.8)

GBM.predict <- predict(GBM.fit,newdata=dat.test,n.trees=GBM.fit$n.trees)

(RMSE <-  sqrt(mean((GBM.predict-dat.test$TotWater)^2)))     
(gini <- NormalizedGini(dat.test$TotWater,GBM.predict))

