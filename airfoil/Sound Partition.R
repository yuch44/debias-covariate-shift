library(readr)
library(SuperLearner)
source("Helper Functions.R")

set.seed(85768)

# Sound partition

dat <- read.table("airfoil_self_noise.dat")
colnames(dat)<- c("Frequency","Angle","Chord","Velocity",
                  "Suction","Sound")
dat.x <- as.matrix(dat[,1:5])
dat.y <- as.numeric(dat[,6])

dat.x[,1] = log(dat.x[,1]) # Log transform
dat.x[,5] = log(dat.x[,5]) # Log transform
N <- nrow(dat.x); p <- ncol(dat.x)

dat.x <- dat.x[order(dat.y),]
dat.y <- sort(dat.y)

n1 <- round(N*0.5); n2 <- N-n1
x1 <- dat.x[1:n1,]; y1 <- dat.y[1:n1]
x2 <- dat.x[-(1:n1),]; y2 <- dat.y[-(1:n1)]

###group flip

flipsize <- round(n1*0.05)
ind1 <- sample(1:n1, flipsize, replace = F)
ind2 <- sample(1:n2, flipsize, replace = F)
x2temp <- x2[ind2,]; y2temp <- y2[ind2]
x1temp <- x1[ind1,]; y1temp <- y1[ind1]
x1[ind1, ]<- x2temp; y1[ind1] <- y2temp
x2[ind2, ] <- x1temp; y2[ind2] <- y1temp

## form data sets
data0 <- as.data.frame(x1)
data0$y <- y1
data0$class <- rep(0,nrow(data0))
data1 <- as.data.frame(x2)
data1$y <- y2
data1$class <- rep(1,nrow(data1))

M <- 500
numFolds <- 2

########################## linear gamma, linear alpha ##########################
ps_linear <- rep(0,M)
values_linear <- rep(0,M)
sd_linear <- rep(0,M)
pb = txtProgressBar(min = 0, max = M, initial = 0) 

for(i in 1:M) {
  setTxtProgressBar(pb,i)
  # sample split
  split0_ind <- sample(seq_len(nrow(data0)),size=floor(0.5 * nrow(data0)))
  split1_ind <- sample(seq_len(nrow(data1)),size=floor(0.5 * nrow(data1)))
  split0 <- data0[split0_ind,]
  split1 <- data1[split1_ind,]
  split <- rbind(split0,split1)
  
  test0 <- data0[-split0_ind,]
  test1 <- data1[-split1_ind,]
  
  # construct a function
  marg_ratio <- estimate_marginal_ratio(split,nrow(split0),nrow(split1),"linear")
  joint_ratio <- estimate_joint_ratio(split,"linear")
  
  a_calc <- function(r0,r1) {
    if(r0<r1){
      return(1)
    } else {
      return(0)
    }
  }
  
  a_outer <- function(point0,point1) {
    cr0 <- joint_ratio(point0) * marg_ratio(point0)
    cr1 <- joint_ratio(point1) * marg_ratio(point1)
    return(outer(cr0,cr1,Vectorize(a_calc)))
  }
  
  # form cross-fit folds
  ind0 <- sample(seq(nrow(test0)))
  ind1 <- sample(seq(nrow(test1)))
  fold0 <- split(ind0, ceiling(seq_along(ind0)/(nrow(test0)/numFolds)))
  fold1 <- split(ind1, ceiling(seq_along(ind1)/(nrow(test1)/numFolds)))
  
  store_values <- matrix(nrow=nrow(test0),ncol = nrow(test1))
  
  for(j in 1:numFolds){
    for(k in 1:numFolds){
      # split data into estimate and nuisance
      est_data0 <- test0[fold0[[j]],]
      est_data1 <- test1[fold1[[k]],]
      est_data <- rbind(est_data0,est_data1)
      nuisance0 <- test0[-fold0[[j]],]
      nuisance1 <- test1[-fold1[[k]],]
      nuisance <- rbind(nuisance0,nuisance1)
      
      # estimate nuisance parameters
      alpha_data <- nuisance0[,1:5]
      as <- a_outer(nuisance0,nuisance1)
      alpha_data$successes <- rowSums(as)
      alpha_data$failures <- nrow(nuisance1) - rowSums(as)
      alpha_model <- glm(cbind(successes,failures)~., data=alpha_data,family=binomial())
      gamma <- estimate_marginal_ratio(nuisance,nrow(nuisance0),nrow(nuisance1),"linear")
      
      # estimate 
      as <- a_outer(est_data0,est_data1)
      gamma0 <- gamma(est_data0)
      alpha0 <- predict(alpha_model,est_data0,type="response")
      alpha1 <- predict(alpha_model,est_data1,type="response")
      for(u in 1:(length(fold0[[j]]))){
        for(v in 1:(length(fold1[[k]]))){
          store_values[fold0[[j]][u],fold1[[k]][v]] <- gamma0[u]*as[u,v] + alpha1[v]-gamma0[u]*alpha0[u]
        }
      }
    }
  }
  
  value <- sum(store_values)/(nrow(test0)*nrow(test1))
  values_linear[i] <- value
  
  #variance estimate
  variance0 <- rep(0,nrow(test0))
  variance1 <- rep(0,nrow(test1))
  
  for(u in 1:(nrow(test0))){
    variance0[u] <- mean(store_values[u,]) - .5
  }
  
  for(u in 1:(nrow(test1))){
    variance1[u] <- mean(store_values[,u]) - .5
  }
  
  variance0 <- sum(variance0^2) / (nrow(test0)-1)
  
  variance1 <- sum(variance1^2) / (nrow(test1)-1)
  
  variance <- (2 *variance0 + 2*variance1) / (nrow(test0)+nrow(test1))
  
  t <- (value-.5) / (sqrt(variance))
  sd_linear[i] <- sqrt(variance)
  ps_linear[i] <- pnorm(t)
}
close(pb)

p_linear <- 2 * median(ps_linear)

####################### linear gamma, superlearner alpha #######################
ps_linear_sl <- rep(0,M)
sd_linear_sl <- rep(0,M)
values_linear_sl <- rep(0,M)

pb = txtProgressBar(min = 0, max = M, initial = 0) 
for(i in 1:M) {
  setTxtProgressBar(pb,i)
  # sample split
  split0_ind <- sample(seq_len(nrow(data0)),size=floor(0.5 * nrow(data0)))
  split1_ind <- sample(seq_len(nrow(data1)),size=floor(0.5 * nrow(data1)))
  split0 <- data0[split0_ind,]
  split1 <- data1[split1_ind,]
  split <- rbind(split0,split1)
  
  test0 <- data0[-split0_ind,]
  test1 <- data1[-split1_ind,]
  
  # construct a function
  marg_ratio <- estimate_marginal_ratio(split,nrow(split0),nrow(split1),"linear")
  joint_ratio <- estimate_joint_ratio(split,"linear")
  
  a_calc <- function(r0,r1) {
    if(r0<r1){
      return(1)
    } else {
      return(0)
    }
  }
  
  a_outer <- function(point0,point1) {
    cr0 <- joint_ratio(point0) * marg_ratio(point0)
    cr1 <- joint_ratio(point1) * marg_ratio(point1)
    return(outer(cr0,cr1,Vectorize(a_calc)))
  }
  
  # form cross-fit folds
  ind0 <- sample(seq(nrow(test0)))
  ind1 <- sample(seq(nrow(test1)))
  fold0 <- split(ind0, ceiling(seq_along(ind0)/(nrow(test0)/numFolds)))
  fold1 <- split(ind1, ceiling(seq_along(ind1)/(nrow(test1)/numFolds)))
  
  store_values <- matrix(nrow=nrow(test0),ncol = nrow(test1))
  
  for(j in 1:numFolds){
    for(k in 1:numFolds){
      # split data into estimate and nuisance
      est_data0 <- test0[fold0[[j]],]
      est_data1 <- test1[fold1[[k]],]
      est_data <- rbind(est_data0,est_data1)
      nuisance0 <- test0[-fold0[[j]],]
      nuisance1 <- test1[-fold1[[k]],]
      nuisance <- rbind(nuisance0,nuisance1)
      
      # estimate nuisance parameters
      alpha_data <- nuisance0[,1:5]
      as <- a_outer(nuisance0,nuisance1)
      as <- rowMeans(as)
      xgb_grid <- create.SL.xgboost(tune = list(ntrees = c(100,200,500,1000), max_depth = c(2,4,6), shrinkage =
                                                  c(0.3), minobspernode = c(1)), detailed_names = F, env = .GlobalEnv,
                                    name_prefix = "SL.xgb")
      
      alpha_model <- SuperLearner(Y=as,X=alpha_data,SL.library = c(c("SL.ranger","SL.lm","SL.ksvm"),xgb_grid$names))
      gamma <- estimate_marginal_ratio(nuisance,nrow(nuisance0),nrow(nuisance1),"linear")
      
      # estimate 
      as <- a_outer(est_data0,est_data1)
      gamma0 <- gamma(est_data0)
      alpha0 <- predict(alpha_model,est_data0[,1:5],onlySL = T)$pred
      alpha1 <- predict(alpha_model,est_data1[,1:5],onlySL = T)$pred
      for(u in 1:(length(fold0[[j]]))){
        for(v in 1:(length(fold1[[k]]))){
          store_values[fold0[[j]][u],fold1[[k]][v]] <- gamma0[u]*as[u,v] + alpha1[v]-gamma0[u]*alpha0[u]
        }
      }
    }
  }
  
  value <- sum(store_values)/(nrow(test0)*nrow(test1))
  values_linear_sl[i] <- value
  
  #variance estimate
  variance0 <- rep(0,nrow(test0))
  variance1 <- rep(0,nrow(test1))
  
  for(u in 1:(nrow(test0))){
    variance0[u] <- mean(store_values[u,]) - .5
  }
  
  for(u in 1:(nrow(test1))){
    variance1[u] <- mean(store_values[,u]) - .5
  }
  
  variance0 <- sum(variance0^2) / (nrow(test0)-1)
  
  variance1 <- sum(variance1^2) / (nrow(test1)-1)
  
  variance <- (2 *variance0 + 2*variance1) / (nrow(test0)+nrow(test1))
  
  t <- (value-.5) / (sqrt(variance))
  sd_linear_sl[i] <- sqrt(variance)
  ps_linear_sl[i] <- pnorm(t)
}
close(pb)

p_linear_sl <- 2 * median(ps_linear_sl)

####################### superlearner gamma, superlearner alpha #######################
ps_sl_sl <- rep(0,M)
sd_sl_sl <- rep(0,M)
values_sl <- rep(0,M)

pb = txtProgressBar(min = 0, max = M, initial = 0) 
for(i in 1:M) {
  setTxtProgressBar(pb,i)
  # sample split
  split0_ind <- sample(seq_len(nrow(data0)),size=floor(0.5 * nrow(data0)))
  split1_ind <- sample(seq_len(nrow(data1)),size=floor(0.5 * nrow(data1)))
  split0 <- data0[split0_ind,]
  split1 <- data1[split1_ind,]
  split <- rbind(split0,split1)
  
  test0 <- data0[-split0_ind,]
  test1 <- data1[-split1_ind,]
  
  # construct a function
  marg_ratio <- estimate_marginal_ratio(split,nrow(split0),nrow(split1),"linear")
  joint_ratio <- estimate_joint_ratio(split,"linear")
  
  a_calc <- function(r0,r1) {
    if(r0<r1){
      return(1)
    } else {
      return(0)
    }
  }
  
  a_outer <- function(point0,point1) {
    cr0 <- joint_ratio(point0) * marg_ratio(point0)
    cr1 <- joint_ratio(point1) * marg_ratio(point1)
    return(outer(cr0,cr1,Vectorize(a_calc)))
  }
  
  # form cross-fit folds
  ind0 <- sample(seq(nrow(test0)))
  ind1 <- sample(seq(nrow(test1)))
  fold0 <- split(ind0, ceiling(seq_along(ind0)/(nrow(test0)/numFolds)))
  fold1 <- split(ind1, ceiling(seq_along(ind1)/(nrow(test1)/numFolds)))
  
  store_values <- matrix(nrow=nrow(test0),ncol = nrow(test1))
  
  for(j in 1:numFolds){
    for(k in 1:numFolds){
      # split data into estimate and nuisance
      est_data0 <- test0[fold0[[j]],]
      est_data1 <- test1[fold1[[k]],]
      est_data <- rbind(est_data0,est_data1)
      nuisance0 <- test0[-fold0[[j]],]
      nuisance1 <- test1[-fold1[[k]],]
      nuisance <- rbind(nuisance0,nuisance1)
      
      # estimate nuisance parameters
      alpha_data <- nuisance0[,1:5]
      as <- a_outer(nuisance0,nuisance1)
      as <- rowMeans(as)
      xgb_grid <- create.SL.xgboost(tune = list(ntrees = c(100,200,500,1000), max_depth = c(2,4,6), shrinkage =
                                                  c(0.3), minobspernode = c(1)), detailed_names = F, env = .GlobalEnv,
                                    name_prefix = "SL.xgb")
      
      alpha_model <- SuperLearner(Y=as,X=alpha_data,SL.library = c(c("SL.ranger","SL.lm","SL.ksvm"),xgb_grid$names))
      gamma <- estimate_marginal_ratio(nuisance,nrow(nuisance0),nrow(nuisance1),"superlearner")
      
      # estimate 
      as <- a_outer(est_data0,est_data1)
      gamma0 <- gamma(est_data0[,1:5])
      alpha0 <- predict(alpha_model,est_data0[,1:5],onlySL = T)$pred
      alpha1 <- predict(alpha_model,est_data1[,1:5],onlySL = T)$pred
      for(u in 1:(length(fold0[[j]]))){
        for(v in 1:(length(fold1[[k]]))){
          store_values[fold0[[j]][u],fold1[[k]][v]] <- gamma0[u]*as[u,v] + alpha1[v]-gamma0[u]*alpha0[u]
        }
      }
    }
  }
  
  value <- sum(store_values)/(nrow(test0)*nrow(test1))
  values_sl[i] <- value
  
  #variance estimate
  variance0 <- rep(0,nrow(test0))
  variance1 <- rep(0,nrow(test1))
  
  for(u in 1:(nrow(test0))){
    variance0[u] <- mean(store_values[u,]) - .5
  }
  
  for(u in 1:(nrow(test1))){
    variance1[u] <- mean(store_values[,u]) - .5
  }
  
  variance0 <- sum(variance0^2) / (nrow(test0)-1)
  
  variance1 <- sum(variance1^2) / (nrow(test1)-1)
  
  variance <- (2 *variance0 + 2*variance1) / (nrow(test0)+nrow(test1))
  
  t <- (value-.5) / (sqrt(variance))
  sd_sl_sl[i] <- sqrt(variance)
  ps_sl_sl[i] <- pnorm(t)
}

close(pb)

p_sl_sl <- 2 * median(ps_sl_sl)

