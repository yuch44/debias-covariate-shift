library(doParallel)
library(foreach)
library(doSNOW)
library(doRNG)

num_cores <- detectCores()

low_dim_lin <- function(M,alpha_set,mu,beta,n,d,s,numFolds,sigma,alt){
  pb <- txtProgressBar(max = M, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  result <- foreach(i =1:M, .options.snow = opts
  ) %dorng% {
    source("Helper Functions.R",local = TRUE)
    data <- rbind(generate_data(s,n,rep(0,d),beta,0,sigma,0,mu),generate_data(s,n,mu,beta,1,sigma,alt,mu))

    split <- rbind(generate_data(s,2*n,rep(0,d),beta,0,sigma,0,mu),generate_data(s,2*n,mu,beta,1,sigma,alt,mu))
    marg_ratio <- estimate_marginal_ratio(split,n,n,"ld")
    joint_ratio <- estimate_joint_ratio(split,"ld")

    a_calc <- function(r0,r1) {
      if(r0 < r1) {
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

    #cross-fit folds
    data0 <- data[1:n,]
    data1 <- data[(n+1):nrow(data),]
    store_values <- matrix(nrow=n,ncol=n)
    ind0 <- sample(seq(n))
    ind1 <- sample(seq(n))
    fold0 <- split(ind0, ceiling(seq_along(ind0)/(n/numFolds)))
    fold1 <- split(ind1, ceiling(seq_along(ind1)/(n/numFolds)))

    for(j in 1:numFolds){
      for(k in 1:numFolds){
        #form crossfit split
        est_data0 <- data0[fold0[[j]],]
        est_data1 <- data1[fold1[[k]],]
        est_data <- rbind(est_data0,est_data1)
        nuisance0 <- data0[-fold0[[j]],]
        nuisance1 <- data1[-fold1[[k]],]
        nuisance <- rbind(nuisance0,nuisance1)

        alpha_data <- nuisance0[,1:d]
        as <- a_outer(nuisance0,nuisance1)
        alpha_data$successes <- rowSums(as)
        alpha_data$failures <- nrow(nuisance1) - rowSums(as)
        alpha_model <- glm(cbind(successes,failures)~.,family="binomial",data=alpha_data)
        gamma <- estimate_marginal_ratio(nuisance,nrow(nuisance0),nrow(nuisance1),"ld")

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

    value <- sum(store_values)/(n*n)


    #Variance estimate

    variance0 <- rep(0,n)
    variance1 <- rep(0,n)

    for(u in 1:(n)){
      variance0[u] <- mean(store_values[u,]) - .5
    }

    for(u in 1:(n)){
      variance1[u] <- mean(store_values[,u]) -.5
    }

    variance0 <- sum(variance0^2) / (n-1)

    variance1 <- sum(variance1^2) / (n-1)


    #compute variance and normalized test-statistic
    variance <- (2 *variance0 + 2*variance1)/(n+n)
    t <- (.5-value) / (sqrt(variance))
    num_rejected <- ifelse(t >= qnorm(1-alpha_set), 1, 0)
    return(c(value,sqrt(variance),num_rejected,t))
  }
  return(result)
}


low_dim_plugin <- function(M,alpha_set,mu,beta,n,d,s,numFolds,sigma,alt,seed){
  pb <- txtProgressBar(max = M, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  result <- foreach(i =1:M, .options.snow = opts) %dorng% {
    source("Helper Functions.R",local = TRUE)
    data <- rbind(generate_data(s,n,rep(0,d),beta,0,sigma,0,mu),generate_data(s,n,mu,beta,1,sigma,alt,mu))
    store_values <- matrix(nrow=n,ncol=n)
    
    split <- rbind(generate_data(s,n,rep(0,d),beta,0,sigma,0,mu),generate_data(s,n,mu,beta,1,sigma,alt,mu))
    marg_ratio <- estimate_marginal_ratio(split,n,n,"ld")
    joint_ratio <- estimate_joint_ratio(split,"ld")
    
    a_calc <- function(r0,r1) {
      if(abs(r0-r1) == 0) {
        return(runif(1))
      } else if(r0 < r1) {
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
    
    #cross-fit folds
    data0 <- data[1:n,]
    data1 <- data[(n+1):nrow(data),]
    ind0 <- sample(seq(n))
    ind1 <- sample(seq(n))
    fold0 <- split(ind0, ceiling(seq_along(ind0)/(n/numFolds)))
    fold1 <- split(ind1, ceiling(seq_along(ind1)/(n/numFolds)))
    
    for(j in 1:numFolds){
      for(k in 1:numFolds){
        #form crossfit split
        est_data0 <- data0[fold0[[j]],]
        est_data1 <- data1[fold1[[k]],]
        est_data <- rbind(est_data0,est_data1)
        nuisance0 <- data0[-fold0[[j]],]
        nuisance1 <- data1[-fold1[[k]],]
        nuisance <- rbind(nuisance0,nuisance1)

        gamma <- estimate_marginal_ratio(nuisance,nrow(nuisance0),nrow(nuisance1),"ld")
        
        #Estimate parameter on other split
        as <- a_outer(est_data0,est_data1)
        gamma0 <- gamma(est_data0)
        
        for(u in 1:(length(fold0[[j]]))){
          for(v in 1:(length(fold1[[k]]))){
            store_values[fold0[[j]][u],fold1[[k]][v]] <- gamma0[u]*as[u,v]
          }
        }
      }
    }
    
    value <- sum(store_values)/(n*n)
    
    
    #Variance estimate
    
    variance0 <- rep(0,n)
    variance1 <- rep(0,n)
    
    for(u in 1:(n)){
      variance0[u] <- mean(store_values[u,]) - .5
    }
    
    for(u in 1:(n)){
      variance1[u] <- mean(store_values[,u]) -.5
    }
    
    variance0 <- sum(variance0^2) / (n-1)
    
    variance1 <- sum(variance1^2) / (n-1)
    
    #compute variance and normalized test-statistic
    variance <- (2 *variance0 + 2*variance1) / (n+n)
    t <- (.5-value) / (sqrt(variance))
    num_rejected <- ifelse(t > qnorm(1-alpha_set), 1, 0)
    return(c(value,sqrt(variance),num_rejected,t))
  }
  return(result)
}

