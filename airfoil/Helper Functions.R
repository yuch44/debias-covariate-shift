library(glmnet)
library(caret)
library(ranger)
library(mgcv)
library(nnet)
library(xgboost)
library(brulee)
library(recipes)
library(kernlab)
library(SuperLearner)

# Estimates marginal density ratio given probability
# @param eta estimated probability of being in class 1
# @param n0 sample size of class 0
# @param n1 sample size of class 1
marg <- function(eta,n0,n1){
  return((n0/n1)*(eta/(1-eta)))
}

# Returns a function that computes estimated marginal density ratio at a point
# @param data data to estimate density ratio
# @param n0 number of points from class 0
# @param n1 number of points from class 1
# @param type "ld" for low dimensional "hd" for high dimensional
estimate_marginal_ratio <- function(data,n0,n1,type){
  if(type == "linear") {
    model <- glm(class ~. -y, data=data, family=binomial())
    marg_ratio <- function(x) {
      eta <- predict(model, newdata = x, type = "response")
      eta[eta<0.01] <- 0.01
      eta[eta>0.99] <- 0.99
      return(sapply(eta,function(x){return(marg(x,n0,n1))}))
    }
  }
  if(type == "superlearner") {
    xgb_grid <- create.SL.xgboost(tune = list(ntrees = c(100,200,500,1000), max_depth = c(2,4,6), shrinkage =
                                                c(0.3), minobspernode = c(1)), 
                                  detailed_names = F, env = .GlobalEnv,
                                  name_prefix = "SL.xgb")
    model <- SuperLearner(Y=data$class,X=data[,!names(data) %in% 
                                                c("class","y")],
                          SL.library = c(c("SL.ranger","SL.lm","SL.ksvm"),
                                         xgb_grid$names),family = binomial())
    marg_ratio <- function(x) {
      eta <- predict(model,x[,!names(x) %in% 
                               c("class","y")],onlySL = T)$pred
      eta[eta<0.01] <- 0.01
      eta[eta>0.99] <- 0.99
      return(sapply(eta,function(x){return(marg(x,n0,n1))}))
    }
  }
  if(type == "nn"){
     rec <- recipe(class ~., data = data[,-5]) %>%
       step_normalize(all_numeric_predictors())
     model <- brulee_mlp(rec,data=data,epochs=15000,hidden_units=c(10,10),learn_rate=0.001,verbose=T)
     marg_ratio <- function(x) {
       eta <- predict(model,x[,1:4],type="prob")$.pred_1
       eta[eta<0.01] <- 0.01
       eta[eta>0.99] <- 0.99
       return(sapply(eta,function(x){return(marg(x,n0,n1))}))
     }
  }
  return(marg_ratio)
}

#Compute joint ratio given probability
joint <- function(eta){
  return((1-eta)/eta)
}

# Returns a function which that computes estimated joint density ratio at given points
# @param data data to estimate density ratio
# @param type "ld" for low-dimensional setting, "hd" for high dimensional
# penalty.factor = c(rep(1,ncol(data)-2),0)
estimate_joint_ratio <- function(data,type){
  if(type == "linear") {
    model <- glm(class ~. , data=data, family=binomial())
    joint_ratio <- function(point) {
      eta <- predict(model, newdata = point, type = "response")
      eta[eta<0.01] <- 0.01
      eta[eta>0.99] <- 0.99
      return(sapply(eta,joint))
    }
  }
  
  if(type == "superlearner") {
    SL.nnet.mod = function(...){
      SL.nnet(...,size=10,maxit=1000)
    }
    xgb_grid <- create.SL.xgboost(tune = list(ntrees = c(100,200,500), max_depth = c(2,6), shrinkage =
                                                c(0.3), minobspernode = c(1)), 
                                  detailed_names = F, env = .GlobalEnv,
                                  name_prefix = "SL.xgb")
    model <- SuperLearner(Y=data$class,X=data[,!names(data) %in% 
                                                c("class")],
                          SL.library = c(c("SL.ranger","SL.lm"),
                                         xgb_grid$names),family = binomial())
    marg_ratio <- function(x) {
      eta <- predict(model,x[,!names(x) %in% 
                               c("class")],onlySL = T)$pred
      eta[eta<0.01] <- 0.01
      eta[eta>0.99] <- 0.99
      return(sapply(eta,joint))
    }
  }
  
  return(joint_ratio)
}

