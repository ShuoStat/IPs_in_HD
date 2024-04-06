library(mboost)
library(doParallel)
source("../pure-transformation/transfo.R")
source("./regressionfun.R")
source("../sources/blandartman&correlation.R")

#- function to get data
load.data <- function(dir, filename){
  dir <- paste0(strsplit(dir, '/')[[1]], collapse = "/")
  dir <- paste0(dir, "/", filename)
  load(dir)
}

dir <- "../data/"
d   <- c("scherzer07.RData", "singh02.RData", "ressom06.RData", 
         "petricoin02.RData", "golub99.RData", "wang05.RData",
          "sotiriou06.RData", "veer02.RData")
link <- c("binomial", "binomial", "binomial", "binomial", 
          "binomial", "binomial", "cox", "cox")

type <- c(Binomial(), Binomial(), Binomial(), Binomial(),
          Binomial(), Binomial(), CoxPH(), CoxPH())

#-------------------------------------------------------------------------------
#- get the mstops and corresponding foldid using cv
#-------------------------------------------------------------------------------
#- prior mstop of each data
#- double amount of mstops in original model
#- 50 if optimal mstops is too small

mstops <- c(50, 818, 7904, 3660, 136, 262, 50, 218)
cv.mstops <- list()
for(i in 1:8){
  
  no.mstop = mstops[i]
  print(i)
  print(Sys.time())
  
  load(paste0(dir, d[i]))
  X <- t(X)
  X <- scale(X)
  n = nrow(X)
  family = type[[i]]
  
  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  set.seed(1)
  nfolds = 10
  foldid <- sample(rep(1:nfolds, length = nrow(X)), nrow(X)) 
  cv <- sapply(1:max(foldid), function(x) as.numeric(foldid != x))
  
  ncores <- detectCores() - 1
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  re.mstops <- foreach(i = 1:(n + 1), .combine = c, 
                       .packages = c("mboost")) %dopar% {
                         
                         boost <- glmboost(X[-i,], y[-i], 
                                           family = family,
                                           control = boost_control(mstop = no.mstop, 
                                                                   nu = 0.1, 
                                                                   risk = "inbag"), 
                                           center = T)
                         
                         cvboost <- cvrisk(boost, folds = cv[-i,])
                         out <- mstop(cvboost)
                         names(out) <- i
                         return(out)
                       }
  
  cv.mstops[[d[i]]] <- list(mstops = re.mstops, foldid= foldid)
  
  registerDoSEQ()
  stopCluster(cl)
}

#- 
# load("./output/cv.mstops.RData")
# save(list = "cv.mstops", file = "./output/cv.mstops.RData")

#-------------------------------------------------------------------------------

load("./output/cv.mstops.RData")

jpeg("./output/tuning.main.points.jpeg", width = 10, height = 3, units = "in",
     res = 300)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
for(s in d[c(4, 5, 8)]){
  
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

jpeg("./output/tuning.append.points.jpeg", width = 10, height = 3, units = "in",
     res = 300)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
for(s in d[c(2, 3, 6)]){
  
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
























