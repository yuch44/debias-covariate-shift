library(readr)
source("Helper Functions.R")
set.seed(598058)

airfoil <- read.table("airfoil_self_noise.dat")
colnames(airfoil)<- c("Frequency","Angle","Chord","Velocity",
                      "Suction","y")

airfoil$Frequency <- log(airfoil$Frequency)
airfoil$Suction <- log(airfoil$Suction)

M<-500
numFolds<-2

#######################Gamma linear, alpha linear###############################

values_linear <- rep(0,M)
rejections_linear <- rep(0,M)

for(i in 1:M){
  ind <- sample(seq_len(nrow(airfoil)), size = 751)
  sample0 <- airfoil[ind,]
  sample1 <- airfoil[-ind,]
  
  ind0 <- sample(seq_len(nrow(sample0)), size = 375)
  ind1 <- sample(seq_len(nrow(sample1)), size = 375)
  data0 <- sample0[-ind0,]
  split0 <- sample0[ind0,]
  data1 <- sample1[-ind1,]
  split1 <- sample1[ind1,]
  split0$class <- rep(0,nrow(split0))
  split1$class <- rep(1,nrow(split1))
  data0$class <- rep(0,nrow(data0))
  data1$class <- rep(1,nrow(data1))
  split <-rbind(split0,split1)
  
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
  ind0 <- sample(seq(nrow(data0)))
  ind1 <- sample(seq(nrow(data1)))
  fold0 <- split(ind0, ceiling(seq_along(ind0)/(nrow(data0)/numFolds)))
  fold1 <- split(ind1, ceiling(seq_along(ind1)/(nrow(data1)/numFolds)))
  
  store_values <- matrix(nrow=nrow(data0),ncol = nrow(data1))
  
  
  for(j in 1:numFolds){
    for(k in 1:numFolds){
      # split data into estimate and nuisance
      est_data0 <- data0[fold0[[j]],]
      est_data1 <- data1[fold1[[k]],]
      est_data <- rbind(est_data0,est_data1)
      nuisance0 <- data0[-fold0[[j]],]
      nuisance1 <- data1[-fold1[[k]],]
      nuisance <- rbind(nuisance0,nuisance1)
      
      #estimate nuisance parameters with nuisance split
      alpha_data <- nuisance0[,1:5]
      as <- a_outer(nuisance0,nuisance1)
      alpha_data$successes <- rowSums(as)
      alpha_data$failures <- nrow(nuisance1) - rowSums(as)
      alpha_model <- glm(cbind(successes,failures)~., data=alpha_data,family=binomial())
      gamma <- estimate_marginal_ratio(nuisance,nrow(nuisance0),nrow(nuisance1),"linear")
      
      #Estimate parameter on other split
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
  
  value <- sum(store_values)/(nrow(data0)*nrow(data1))
  values_linear[i] <- value
  
  #variance estimate
  variance0 <- rep(0,nrow(data0))
  variance1 <- rep(0,nrow(data1))
  
  for(u in 1:(nrow(data0))){
    variance0[u] <- mean(store_values[u,]) - .5
  }
  
  for(u in 1:(nrow(data1))){
    variance1[u] <- mean(store_values[,u]) - .5
  }
  
  variance0 <- sum(variance0^2) / (nrow(data0)-1)
  
  variance1 <- sum(variance1^2) / (nrow(data1)-1)
  
  variance <- (2 *variance0 + 2*variance1) / (nrow(data0)+nrow(data1))
  
  t <- (.5-value) / (sqrt(variance))
  rejections_linear[i] <- ifelse(t >= qnorm(1-.05), 1, 0)
}

#############################Gamma sl, alpha sl#################################

values_sl <- rep(0,M)
rejections_sl<- rep(0,M)

for(i in 1:M){
  ind <- sample(seq_len(nrow(airfoil)), size = 751)
  sample0 <- airfoil[ind,]
  sample1 <- airfoil[-ind,]

  ind0 <- sample(seq_len(nrow(sample0)), size = 375)
  ind1 <- sample(seq_len(nrow(sample1)), size = 375)
  data0 <- sample0[-ind0,]
  split0 <- sample0[ind0,]
  data1 <- sample1[-ind1,]
  split1 <- sample1[ind1,]
  split0$class <- rep(0,nrow(split0))
  split1$class <- rep(1,nrow(split1))
  data0$class <- rep(0,nrow(data0))
  data1$class <- rep(1,nrow(data1))
  split <-rbind(split0,split1)

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
  ind0 <- sample(seq(nrow(data0)))
  ind1 <- sample(seq(nrow(data1)))
  fold0 <- split(ind0, ceiling(seq_along(ind0)/(nrow(data0)/numFolds)))
  fold1 <- split(ind1, ceiling(seq_along(ind1)/(nrow(data1)/numFolds)))

  store_values <- matrix(nrow=nrow(data0),ncol = nrow(data1))


  for(j in 1:numFolds){
    for(k in 1:numFolds){
      # split data into estimate and nuisance
      est_data0 <- data0[fold0[[j]],]
      est_data1 <- data1[fold1[[k]],]
      est_data <- rbind(est_data0,est_data1)
      nuisance0 <- data0[-fold0[[j]],]
      nuisance1 <- data1[-fold1[[k]],]
      nuisance <- rbind(nuisance0,nuisance1)

      #estimate nuisance parameters with nuisance split
      alpha_data <- nuisance0[,1:5]
      as <- a_outer(nuisance0,nuisance1)
      as <- rowMeans(as)
      xgb_grid <- create.SL.xgboost(tune = list(ntrees = c(100,200,500), max_depth = c(2,4,6), shrinkage =
                                                  c(0.3), minobspernode = c(1)), detailed_names = F, env = .GlobalEnv,
                                    name_prefix = "SL.xgb")

      alpha_model <- SuperLearner(Y=as,X=alpha_data,SL.library = c(c("SL.ranger","SL.lm"),xgb_grid$names))
      gamma <- estimate_marginal_ratio(nuisance,nrow(nuisance0),nrow(nuisance1),"superlearner")

      #Estimate parameter on other split
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

  value <- sum(store_values)/(nrow(data0)*nrow(data1))
  values_sl[i] <- value

  #variance estimate
  variance0 <- rep(0,nrow(data0))
  variance1 <- rep(0,nrow(data1))

  for(u in 1:(nrow(data0))){
    variance0[u] <- mean(store_values[u,]) - .5
  }

  for(u in 1:(nrow(data1))){
    variance1[u] <- mean(store_values[,u]) - .5
  }

  variance0 <- sum(variance0^2) / (nrow(data0)-1)

  variance1 <- sum(variance1^2) / (nrow(data1)-1)

  variance <- (2 *variance0 + 2*variance1) / (nrow(data0)+nrow(data1))

  t <- (.5-value) / (sqrt(variance))
  rejections_sl[i] <- ifelse(t >= qnorm(1-.05), 1, 0)
}

rej_per_linear <- mean(rejections_linear)
rej_per_sl <- mean(rejections_sl)
