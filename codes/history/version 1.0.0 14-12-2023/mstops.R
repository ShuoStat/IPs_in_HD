
source("./fun.R")

#- loading packages
load_packages(c("mboost", "doParallel", "dplyr", "survival"))

#-------------------------------------------------------------------------------
#- Parameters
#-------------------------------------------------------------------------------

dir <- "../../../0- data/"
d <- c("singh02.RData", "ressom06.RData", "petricoin02.RData", "golub99.RData", 
       "wang05.RData", "veer02.RData")

# link function
link <- c("binomial", "binomial", "binomial", "binomial", "binomial", "cox")

#- family
type <- c("binomial" = Binomial(), "cox" =  CoxPH())[link]
  
#-------------------------------------------------------------------------------
#- get the mstops and corresponding foldid using cv
#-------------------------------------------------------------------------------
#- prior mstop of each data
#- double amount of mstops in original model
#- 50 if optimal mstops is too small

mstops <- c(818, 7904, 3660, 136, 262, 218)
cv.mstops <- list()

for(i in 1:6){
  
  no.mstop = mstops[i]
  print(d[i])
  print(Sys.time())
  
  #- load data
  load(paste0(dir, d[i]))
  
  X <- t(X)
  X <- scale(X)
  n = nrow(X)
  family = type[[i]]
  
  if (link[i] == "binomial") 
    y = as.factor(y)
  if (link[i] == "cox")
    y = surv
  
  set.seed(1)
  nfolds = 10
  foldid <- sample(rep(1:nfolds, length = nrow(X)), nrow(X)) 
  cv <- sapply(1:max(foldid), function(x) as.numeric(foldid != x))
  
  ncores <- detectCores() - 1
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  re.mstops <- foreach(k = 1:(n + 1), 
                       .packages = c("mboost")) %dopar% {
                         
                         time1 <- Sys.time()
                         w    <- rep(1, n)
                         if (k <= n) w[k] <- 0
                         boost <- glmboost(X, y, 
                                           family = family,
                                           weights = w,
                                           control = boost_control(mstop = no.mstop, 
                                                                   nu = 0.1, 
                                                                   risk = "inbag"), 
                                           center = F)
                         
                         cvboost <- cvrisk1(boost, folds = cv)

                         out <- mstop(cvboost)
                         time2 <- Sys.time()
                         out <- list(cvboost = cvboost,
                                     "id" = k,
                                     "mstop" = out, 
                                     "time" = as.numeric(difftime(time2, time1, unit = "secs")))

                         return(out)
                       }
  
  
  cv.mstops[[d[i]]] <- list(mstops = re.mstops, foldid = foldid)
  
  registerDoSEQ()
  stopCluster(cl)
}

#-
# save(list = "cv.mstops", file = "../output/cv.mstops.RData")
#-------------------------------------------------------------------------------
#- prediction
#-------------------------------------------------------------------------------

#- loo, leave one out

cv.pred <- function(X, y, family, folds, loo = 0){
  
  n <- nrow(X)
  res <- c()
  for(i in seq_len(ncol(folds))){
    
    f <- folds[,i]
    if (loo <= n)
      f[loo] <- 0
    
    boost <- glmboost(X[f == 1, ], y[f == 1], 
                      family = family,
                      control = boost_control(mstop = no.mstop, 
                                              nu = 0.1, 
                                              risk = "inbag"), 
                      center = F)
    
    out <-  folds[,i] == 0
    yy <- y[out]
    
    re <- mapply(function(k){
      p <- boost[k]$predict(newdata = X[out,])
      Binomial()@risk(yy, p)}, 
      k = 1:no.mstop)
    
    res <- rbind(res, re)
  }
  return(res)
}

mstops <- c(818, 7904, 3660, 136, 262, 218)
cv.mstops <- list()

for(i in 1:6){
  
  no.mstop = mstops[i]
  print(d[i])
  print(Sys.time())
  
  #- load data
  load(paste0(dir, d[i]))
  
  X <- t(X)
  X <- scale(X)
  n = nrow(X)
  family = type[[i]]
  
  if (link[i] == "binomial") 
    y = as.factor(y)
  if (link[i] == "cox")
    y = surv
  
  set.seed(1)
  nfolds = 10
  foldid <- sample(rep(1:nfolds, length = nrow(X)), nrow(X)) 
  cv <- sapply(1:max(foldid), function(x) as.numeric(foldid != x))
  
  ncores <- detectCores() - 1
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  re.loo <- foreach(k = 1:(n + 1), 
                    .packages = c("mboost")) %dopar% {
                      
                      re <- cv.pred(X, y, family = family, folds = folds, loo = 200)
                      colSums(re)
                    }
  
  load("../output/cv.mstops.RData")
  mstops <- cv.mstops$golub99.RData
  
  re <- c()
  for(i in 1:n){
    
    tmp <- re.loo[[i]][m2[i]] - re.loo[[n + 1]][m2[n + 1]]
    re  <- c(re, tmp)
  }

  #- likelihood displacement

  re <- c()
  for(i in seq_len(n + 1)){
    re <- c(re, colSums(re.mstops[[i]]$cvboost)[m2[i]])
  }
  
  
  re - re[n + 1]
  order(abs(re - re[n + 1]), decreasing = T)
  
  re <- c()
  for(i in seq_len(n + 1)){
    re <- c(re, colSums(re.mstops[[i]]$cvboost)[i + 1])
  }
  
  re - re[n + 1]
  cv.mstops[[d[i]]] <- list(mstops = re.mstops, foldid = foldid)
  
  registerDoSEQ()
  stopCluster(cl)
}

#-------------------------------------------------------------------------------


load("../output/cv.mstops.RData")
jpeg("./output/mstops.jpeg", width = 10, height = 3, units = "in",
     res = 300)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
for(s in d[c(3, 4, 6)]){
  
  obj <- cv.mstops[[s]][["mstops"]]
  drop1 <- obj[-length(obj)]
  drop0 <- obj[length(obj)]
  score <- abs(drop1 - drop0)
  
  x = 1:length(drop1)
  y = drop1
  plot(x, y, xlab = "", ylab = "", cex = 0.8)
  abline(h = drop0)
  
  for(i in 1:length(x)){
    lines(c(x[i], x[i]), c(drop0, y[i]), col = "grey")
  }
  points(x, y, xlab = "", ylab = "", pch = 19, cex = 0.8)
  mtext(gsub(".RData", "", s), side = 3, line = 0.2, adj = 0)
  mtext(side = 2, line = 0.2, at = drop0, text = "orig", cex = 0.8, las = 1)
  ord <- order(score, decreasing = T)[1:5]
  text(ord, drop1[ord], pos = 4, cex = 0.9, labels = ord)
  mtext("mstops", side = 2, line = 2.2)
  mtext("observations", side = 1, line = 2.2)
  
}
dev.off()

#-------------------------------------------------------------------------------

jpeg("./output/mstops.A.jpeg", width = 10, height = 3, units = "in",
     res = 300) #- A, for appendix

par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))

for(s in d[c(1, 2, 5)]){
  
  obj <- cv.mstops[[s]][["mstops"]]
  drop1 <- obj[-length(obj)]
  drop0 <- obj[length(obj)]
  score <- abs(drop1 - drop0)
  
  x = 1:length(drop1)
  y = drop1
  plot(x, y, xlab = "", ylab = "", cex = 0.9)
  abline(h = drop0)
  
  for(i in 1:length(x)){
    lines(c(x[i], x[i]), c(drop0, y[i]), col = "grey")
  }
  points(x, y, xlab = "", ylab = "", pch = 19, cex = 0.8)
  mtext(gsub(".RData", "", s), side = 3, line = 0.2, adj = 0)
  mtext(side = 2, line = 0.2, at = drop0, text = "orig", cex = 0.8, las = 1)
  ord <- order(score, decreasing = T)[1:3]
  text(ord, drop1[ord], pos = 4, cex = 0.8, labels = ord)
  mtext("mstops", side = 2, line = 2.2)
  mtext("observations", side = 1, line = 2.2)
  
}
dev.off()

#-------------------------------------------------------------------------------
