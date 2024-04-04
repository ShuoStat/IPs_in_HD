
source("./fun.R")

#- loading packages
load_packages(c("glmnet", "doParallel", "dplyr", "survival"))

#-------------------------------------------------------------------------------

glm.loo <- function(X, y, folds, loo = 0, lambda){
  
  n <- nrow(X)
  f <- unique(folds) 
  
  mse.loo <- c()
  mse.ful <- c()
  
  for(i in f){
    
    s <- folds != i
    l <- rep(T, n)
    if (loo <= n)
      l[loo] <- F 
    
    mod    <- glmnet(X[s & l,], y[s & l], lambda = lambda)
    p.loo  <- predict(mod, X[!s &  l,])
    p.ful  <- predict(mod, X[!s,])
    
    y.loo <- y[!s & l]
    y.ful <- y[!s]
    
    re.loo <- colSums((p.loo - y.loo)^2)
    re.ful <- colSums((p.ful - y.ful)^2)
    
    mse.loo <- rbind(mse.loo, re.loo)
    mse.ful <- rbind(mse.ful, re.ful)
    
  }
  
  return(list(n = loo, mse.loo = colSums(mse.loo), mse.ful = colSums(mse.ful))) 
}



#- use
set.seed(1)
dat <- read.csv("../../../0- data/edu_bodyfat_std.csv")
dat <- as.matrix(dat)
n <- nrow(dat)
folds <- sample(rep(1:10, length = n), n)

X <- dat[,-1]
y <- dat[,1]


loo <- c()
ful <- c()

lambda <- glmnet(X, y)$lambda

for(i in 1:(n + 1)){
  
  re  <- glm.loo(X, y, folds, loo = i, lambda = lambda)
  print(re$n)
  loo <- rbind(loo, re$mse.loo)
  ful <- rbind(ful, re$mse.ful)

}

sloo <- apply(loo, 1, which.min)

sful <- c()

for(i in 1:(n + 1)){
  sful <- c(sful, ful[i, sloo[i]])
}
  
length(sful)
sful - sful[252]










