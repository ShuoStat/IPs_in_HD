
source("./fun.R")

#- loading packages
load_packages(c("mboost", "doParallel", "dplyr", "survival"))

#-------------------------------------------------------------------------------
#- Parameters
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
  
  re.mstops <- foreach(i = 1:(n + 1), 
                       .packages = c("mboost")) %dopar% {
                         
                         time1 <- Sys.time()
                         boost <- glmboost(X[-i,], y[-i], 
                                           family = family,
                                           control = boost_control(mstop = no.mstop, 
                                                                   nu = 0.1, 
                                                                   risk = "inbag"), 
                                           center = T)
                         
                         cvboost <- cvrisk(boost, folds = cv[-i,])
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

#-

cvpred <- function(boost, folds, newdata) {
  
  boost <- glmboost(X[-i,], y[-i], 
                    family = family,
                    control = boost_control(mstop = no.mstop, 
                                            nu = 0.1, 
                                            risk = "inbag"), 
                    center = T)
  
  
  
  
  
  
  
  
  
  
}




# save(list = "cv.mstops", file = "./output/cv.mstops.RData")
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


