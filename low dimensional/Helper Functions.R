library(MASS)
library(glmnet)
library(caret)
library(doRNG)
library(mvtnorm)

generate_x <- function(n, mu, Sigma,s,mu1){
  x <- mvrnorm(2*n, mu=mu, Sigma = Sigma)
  ind <- which(g.orac(x, mu1,Sigma)<50&g.orac(x, mu1,Sigma)>1/50)
  x <- x[ind,]
  n_res <- n - length(ind)
  while(n_res>0){
    xx <- mvrnorm(2*n_res, mu=mu, Sigma = Sigma)
    ind <- which(g.orac(xx,mu1,Sigma)<50&g.orac(xx,mu1,Sigma)>1/50)
    x <- rbind(x, xx[ind,])
    n_res <- n_res - length(ind)
  }
  return(x[1:n,])
}

generate_y <- function(x, beta, c){
  n = dim(x)[1]
  res = rnorm(n)
  y = tcrossprod(x, t(beta)) + res + c
  return(y)
}

# Generates a multivariate normal dataset with the following parameters
# @param p dimension
# @param s sparse dimension of y
# @param n number of datapoints for each sample
# @param mu mean
# @param beta coefficients of linear model
# @param class whether this is sample 0 or 1

generate_data <- function(s,n,mu,beta,class,sigma,alt,mu1) {
  Sigma <- sigma * diag(length(mu))
  data <- generate_x(n,mu,Sigma,s,mu1)
  ys <- generate_y(data,beta,alt)
  data <- as.data.frame(data)
  data$y <- ys
  data$class <- rep(class,n)
  return(data)
}

#oracle gamma function
g.orac <- function(x,mu,Sigma){
  x<-as.matrix(x)
  return(dmvnorm(x,mean=mu,sigma = Sigma)/dmvnorm(x,mean=rep(0,length(mu)),sigma = Sigma))
}


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
  if(type == "hd"){
    X <- model.matrix(~.  - y - class-1, data = data)
    model <- cv.glmnet(X,data$class, family="binomial",alpha = 1)
    marg_ratio <- function(x) {
      eta <- predict(model, newx = x, type = "response",s="lambda.min")
      return(sapply(eta,function(x){return(marg(x,n0,n1))}))
    }
  }
  if(type == "ld") {
    model <- glm(class ~. -y, data=data, family=binomial())
    marg_ratio <- function(x) {
      eta <- predict(model, newdata = x, type = "response")
      return(sapply(eta,function(x){return(marg(x,n0,n1))}))
    }
  }
  return(marg_ratio)
}


#Compute joint ratio given proabability
joint <- function(eta){
  return((1-eta)/eta)
}

# Returns a function which that computes estimated joint density ratio at given points
# @param data data to estimate density ratio
# @param type "ld" for low-dimensional setting, "hd" for high dimensional
# penalty.factor = c(rep(1,ncol(data)-2),0) force y
estimate_joint_ratio <- function(data,type){
  if(type == "hd") {
    X <- model.matrix(~.  - class-1, data = data)
    model <- cv.glmnet(X,data$class, family="binomial",alpha = 1,penalty.factor = c(rep(1,ncol(data)-2),0))
    joint_ratio <- function(point) {
      eta <- predict(model, newx = point, type = "response",s="lambda.min")
      return(sapply(eta,joint))
    }
  }
  if(type == "ld") {
    model <- glm(class ~. , data=data, family=binomial())
    joint_ratio <- function(point) {
      eta <- predict(model, newdata = point, type = "response")
      return(sapply(eta,joint))
    }
  }
  return(joint_ratio)
}

# The oracle alpha nuisance function
oracle_alpha <- function(mu,eta1,eta2,beta) {
  estimate <- function(x) {
    x <- as.numeric(x)
    mean <- eta1 + (eta2*beta)
    mean <- mean %*% mu
    var <- eta1 + (eta2*beta)
    var <- var^2
    var <- sum(var) + (2*(eta2^2))
    num <- eta1 %*% x + (eta2 * (beta %*% x))
    num <- num - mean
    return(1-pnorm(num/sqrt(var)))
  }
  if(sum(eta1^2)==0 & eta2 == 0){
    estimate <- function(x){
      return(.5)
    }
  }
  return(estimate)
}

