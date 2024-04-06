
source("./fun.R")

#- loading packages
load_packages(c("mboost", "doParallel", "dplyr", "survival", "doRNG", "doSNOW"))

#-------------------------------------------------------------------------------

dir <- "../../../0- data/"
d   <- c("singh02.RData", "ressom06.RData", "petricoin02.RData", 
         "golub99.RData", "wang05.RData", "veer02.RData")

# link function
link <- c("binomial", "binomial", "binomial", 
          "binomial", "binomial", "cox")

#- family
type <- c("binomial" = Binomial(), "cox" =  CoxPH())[link]

#-------------------------------------------------------------------------------
#- Likelihood displacement
#-------------------------------------------------------------------------------

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
  
  re.mstops <- foreach(i = 1:(n + 1), 
                       .packages = c("mboost")) %dopar% {
                         
                         time1 <- Sys.time()
                         w    <- rep(1, n)
                         if (i <= n) w[i] <- 0
                         boost <- glmboost(X, y, 
                                           family = family,
                                           weights = w,
                                           control = boost_control(mstop = no.mstop, 
                                                                   nu = 0.1, 
                                                                   risk = "inbag"), 
                                           center = F)
                         
                         cvboost <- cvrisk(boost, folds = cv)
                         
                         out <- mstop(cvboost)
                         time2 <- Sys.time()
                         out <- list(cvboost = cvboost,
                                     "id" = i,
                                     "mstop" = out, 
                                     "time" = as.numeric(difftime(time2, time1, unit = "secs")))
                         
                         return(out)
                       }
  
  
  cv.mstops[[d[i]]] <- list(mstops = re.mstops, foldid = foldid)
  
  registerDoSEQ()
  stopCluster(cl)
}


#-------------------------------------------------------------------------------
#- Prediction: adaptive mstop
#-------------------------------------------------------------------------------

load("../output/re.cv.pred.RData")

plot.pred <- function(obj, mstops, top = 3, ylim = NULL){
  
  n <- length(obj)
  obj <- lapply(obj, colMeans)
  obj <- Reduce(cbind, obj) 
  
  if (length(mstops) == 1) 
    mstops <- rep(mstops, n - 1)
  
  mstops <- mstops + 1 #- because of the 0th step
  
  # m   <- rowMeans(obj[,-n])
  # obj <- apply(obj[,-n], 2, function(x) abs(x - m))
  
  obj <- apply(obj[,-n], 2, function(x) abs(x - obj[,n]))
  scores <- sapply(1:(n - 1), function(x) obj[mstops[x], x])
  scores <- scores / sd(scores)
  
  x = 1:length(scores)
  y = scores
  
  if(is.null(ylim))
    ylim = c(0, max(y) * 1.2) 
  
  plot(x, y, xlab = "", ylab = "", cex = 0.8, yaxs = "i", 
       ylim = ylim)
  
  abline(h = 0)
  for(i in seq_along(x)){
    lines(c(x[i], x[i]), c(0, y[i]), col = "grey")
  }
  
  points(x, y, xlab = "", ylab = "", pch = 19, cex = 0.8)
  mtext("Df-Cvpath", side = 2, line = 2.2)
  mtext("Observations", side = 1, line = 2.2)
  ord <- order(abs(scores), decreasing = T)[1:5]
  text(ord, y[ord], pos = 3, cex = 0.8, labels = x[ord])
}


jpeg("../output/pred.jpeg", width = 9, height = 3, units = "in", res = 300)
par(mar = c(3, 3, 2, 1), oma = c(1, 1, 1, 1), mfrow = c(1, 3))

s <- c("petricoin02.RData", "golub99.RData", "veer02.RData")
load("../output/cv.mstops.RData")

for(i in s){
  
  obj <- re.cv.pred[[i]]
  #- for adaptive mstops
  plot.pred(obj$cv.pred, mstops = cv.mstops[[i]])
  
  #- for fixed mstops
  plot.pred(obj$cv.pred, mstops = rev(cv.mstops[[i]])[1])
  
  mtext(gsub(".RData", "", i), adj = 0, line = 0.2, side = 3)
  
}

dev.off()

#-------------------------------------------------------------------------------

jpeg("../output/pred.A.jpeg", width = 9, height = 3, units = "in", 
     res = 300)
par(mar = c(3, 3, 2, 1), oma = c(1, 1, 1, 1), mfrow = c(1, 3))
s <- c("singh02.RData", "wang05.RData", "ressom06.RData")
for(i in s){
  obj <- re.cv.pred[[i]]
  plot.pred(obj$cv.pred)
  mtext(gsub(".RData", "", i), adj = 0, line = 0.2, side = 3)
}

dev.off()


#-------------------------------------------------------------------------------
#- Prediction: cv.path
#-------------------------------------------------------------------------------

plot.pred.path <- function(obj, top = 3){
  
  n <- length(obj)
  obj <- lapply(obj, colMeans)
  obj <- Reduce(cbind, obj)
  
  # m   <- rowMeans(obj[,-n])
  # obj <- apply(obj[,-n], 2, function(x) abs(x - m))
  obj <- apply(obj[,-n], 2, function(x) abs(x - obj[,n]))
  
  ylim <- c(0, max(obj))
  mstops <- nrow(obj)
  x = 1:mstops
  
  plot.new()
  plot.window(xlim = c(0, nrow(obj)), ylim = ylim)
  axis(1)
  axis(2)
  box()
  ord <- order(colSums(obj), decreasing = T)[1:top]
  
  for(i in 1:ncol(obj)){
    lines(x, obj[,i], col = "grey")
  }
  
  for(i in ord){
    lines(x, obj[,i], col = "black")
  }
  
  #- add text
  for(i in seq_along(ord)){
    text(mstops, obj[mstops, ord[i]], labels = ord[i], pos = 1, cex = 0.8)  
  }
  mtext("Number of iterations", side = 1, line = 2.2)
  mtext(bquote(Delta ~ "cvm"), side = 2, line = 2.2)
}

plot.pred.score <- function(obj, top = 3, ylim = NULL){
  
  n <- length(obj)
  obj <- lapply(obj, colMeans)
  obj <- Reduce(cbind, obj)
  
  # m   <- rowMeans(obj[,-n])
  # obj <- apply(obj[,-n], 2, function(x) abs(x - m))
  obj <- apply(obj[,-n], 2, function(x) abs(x - obj[,n]))
  
  obj <- colMeans(obj)
  score <- obj / sd(obj)
  
  x = 1:length(score)
  y = score
  
  if(is.null(ylim))
    ylim = c(0, max(y) * 1.2) 
  plot(x, y, xlab = "", ylab = "", cex = 0.8, yaxs = "i", 
       ylim = ylim)
  
  abline(h = 0)
  for(i in seq_along(x)){
    lines(c(x[i], x[i]), c(0, y[i]), col = "grey")
  }
  points(x, y, xlab = "", ylab = "", pch = 19, cex = 0.8)
  mtext("Df-Cvpath", side = 2, line = 2.2)
  mtext("Observations", side = 1, line = 2.2)
  ord <- order(abs(score), decreasing = T)[1:5]
  text(ord, y[ord], pos = 3, cex = 0.8, labels = x[ord])

}

load("../output/re.cv.pred.RData")

jpeg("../output/cvpath.jpeg", width = 9, height = 5, units = "in", res = 300)
layout(matrix(1:6, nrow = 2, byrow = F))
par(mar = c(3, 3, 2, 1), oma = c(1, 1, 1, 1))

s <- c("petricoin02.RData", "golub99.RData", "veer02.RData")

for(i in s){
  
  obj <- re.cv.pred[[i]]
  plot.pred.path(obj$cv.pred)
  mtext(gsub(".RData", "", i), adj = 0, line = 0.2, side = 3)
  plot.pred.score(obj$cv.pred, ylim = c(0, 13))
}

dev.off()

jpeg("../output/cvpath.A.jpeg", width = 9, height = 5, units = "in", 
     res = 300)

mat <- matrix(1:6, nrow = 2, byrow = F)
layout(mat)
par(mar = c(3, 3, 2, 1), oma = c(1, 1, 1, 1))
s <- c("singh02.RData", "wang05.RData", "ressom06.RData")
for(i in s){
  obj <- re.cv.pred[[i]]
  plot.pred.path(obj$cv.pred)
  mtext(gsub(".RData", "", i), adj = 0, line = 0.2, side = 3)
  plot.pred.score(obj$cv.pred, ylim = c(0, 13))
}
dev.off()


#-------------------------------------------------------------------------------





