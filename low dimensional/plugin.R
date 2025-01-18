source("Helper Functions.R")
source("Tests.R")

num_cores <- detectCores()
cl <- makeCluster(num_cores-2)
registerDoSNOW(cl)

# null

# 250
M <- 500
alpha_set <- 0.05
mu <- c(1,-1,1,-1,0)
beta <- c(1,1,1,1,1)
n <- 250
d <- 5
numFolds <- 2
s <- 5
sigma <- 1
alt <- 0


results_250 <- low_dim_plugin(M,alpha_set,mu,beta,n,d,s,numFolds,sigma,alt)
results_250 <- as.data.frame(do.call(rbind, results_250))

# 500
M <- 500
alpha_set <- 0.05
mu <- c(1,-1,1,-1,0)
beta <- c(1,1,1,1,1)
n <- 500
d <- 5
numFolds <- 2
s <- 5
sigma <- 1
alt <- 0


results_500 <- low_dim_plugin(M,alpha_set,mu,beta,n,d,s,numFolds,sigma,alt)
results_500 <- as.data.frame(do.call(rbind, results_500))

#1000
M <- 500
alpha_set <- 0.05
mu <- c(1,-1,1,-1,0)
beta <- c(1,1,1,1,1)
n <- 1000
d <- 5
numFolds <- 2
s <- 5
sigma <- 1
alt <- 0

results_1000 <- low_dim_plugin(M,alpha_set,mu,beta,n,d,s,numFolds,sigma,alt)
results_1000 <- as.data.frame(do.call(rbind, results_1000))

#2000
M <- 500
alpha_set <- 0.05
mu <- c(1,-1,1,-1,0)
beta <- c(1,1,1,1,1)
n <- 2000
d <- 5
numFolds <- 2
s <- 5
sigma <- 1
alt <- 0

results_2000 <- low_dim_plugin(M,alpha_set,mu,beta,n,d,s,numFolds,sigma,alt)
results_2000 <- as.data.frame(do.call(rbind, results_2000))


#Alternative

#250
M <- 500
alpha_set <- 0.05
mu <- c(1,-1,1,-1,0)
beta <- c(1,1,1,1,1)
n <- 250
d <- 5
numFolds <- 2
s <- 5
sigma <- 1
alt <- 0.5


results_250_a <- low_dim_plugin(M,alpha_set,mu,beta,n,d,s,numFolds,sigma,alt)
results_250_a <- as.data.frame(do.call(rbind, results_250_a))

# 500
M <- 500
alpha_set <- 0.05
mu <- c(1,-1,1,-1,0)
beta <- c(1,1,1,1,1)
n <- 500
d <- 5
s <- 5
sigma <- 1
alt <- .5
numFolds <- 2

results_500_a <- low_dim_plugin(M,alpha_set,mu,beta,n,d,s,numFolds,sigma,alt)
results_500_a <- as.data.frame(do.call(rbind, results_500_a))

#1000
M <- 500
alpha_set <- 0.05
mu <- c(1,-1,1,-1,0)
beta <- c(1,1,1,1,1)
n <- 1000
d <- 5
s <- 5
sigma <- 1
alt <- .5
numFolds <- 2

results_1000_a <- low_dim_plugin(M,alpha_set,mu,beta,n,d,s,numFolds,sigma,alt)
results_1000_a <- as.data.frame(do.call(rbind, results_1000_a))

#2000
M <- 500
alpha_set <- 0.05
mu <- c(1,-1,1,-1,0)
beta <- c(1,1,1,1,1)
n <- 2000
d <- 5
s <- 5
sigma <- 1
alt <- .5
numFolds <- 2

results_2000_a <- low_dim_plugin(M,alpha_set,mu,beta,n,d,s,numFolds,sigma,alt)
results_2000_a <- as.data.frame(do.call(rbind, results_2000_a))

stopCluster(cl)
