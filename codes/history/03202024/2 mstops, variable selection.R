
source("./fun.R")

#- loading packages
load_packages(c("mboost", "doParallel", "dplyr", "survival"))

#-------------------------------------------------------------------------------
#- Parameters
#-------------------------------------------------------------------------------

d <- c("singh02.RData", "ressom06.RData", "petricoin02.RData", 
       "golub99.RData", "wang05.RData", "veer02.RData")

# link function
link <- c("binomial", "binomial", "binomial", "binomial", "binomial", "cox")

#- family
type <- c("binomial" = Binomial(), "cox" =  CoxPH())[link]
  
#-------------------------------------------------------------------------------
#- get the mstops and corresponding foldid using cv
#-------------------------------------------------------------------------------
#- prior mstop of each data
#- double amount of mstops iin original model
#- 50 if optimal mstops is too small

mstopSet <- c(818, 3230, 1272, 136, 262, 218)
nu     <- c(0.1, 0.3, 0.3, 0.1, 0.1, 0.1)
Ncores <- c(20, 8, 20, 20, 20, 20)

for(i in 1:6){
  
  print(d[i])
  print(Sys.time())
  
  #- load data
  load(paste0("../data/", d[i]))
  
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
  
  ncores <- pmin(detectCores() - 1, Ncores[i])
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  re.mstops <- foreach(k = 1:(n + 1), 
                       .packages = c("mboost")) %dopar% {
                         
                         time1 <- Sys.time()
                         boost <- glmboost(X[-k,], y[-k], 
                                           family = family,
                                           control = boost_control(mstop = mstopSet[i], 
                                                                   nu = nu[i], 
                                                                   risk = "inbag"), 
                                           center = F)
                         
                         cvboost <- cvrisk(boost, folds = cv[-k,])
                         
                         out <- mstop(cvboost)
                         time2 <- Sys.time()
                         out <- list(cvboost = cvboost,
                                     "id" = k,
                                     "mstop" = out, 
                                     "time" = as.numeric(difftime(time2, time1, unit = "secs")))
                         
                         return(out)
                       }
  
  mstops <- list(mstops = re.mstops, foldid = foldid)
  save(list = "mstops", file = paste0("../output/mstop_", d[i]))
  
  registerDoSEQ()
  stopCluster(cl)
}

#-------------------------------------------------------------------------------
#- IPs on variable selection
#-------------------------------------------------------------------------------

drop1vs <- function(datName, family, nu, method = c("fixed", "adaptive"), 
                    nCores = 1){
  
  # method, "fixed", "adaptive"
  # dropvs, drop1 for variable selection
  load(paste0("../output/mstop_", datName))
  load(paste0("../data/", datName))
  
  if (class(family) == "boost_family_glm")
    y = as.factor(y)
  if (class(family) == "boost_family")
    y = surv

  # get mstops
  mstops <- unlist(lapply(mstops$mstops, `[[`, "mstop"))
  
  # mstops
  origMstops = mstops[length(mstops)]
  n <- ncol(X)
  p <- nrow(X)
  X <- scale(t(X))
  
  if (method[1] == "fixed") 
    mstops = rep(origMstops, (n + 1))
  if (method[1] == "adaptive") 
    mstops = mstops
  
  # 
  cl <- makeCluster(pmin(detectCores() - 1, nCores))
  registerDoParallel(cl)
  
  #- (n + 1): boosting without removing any case
  reboost <- foreach(i = 1:(n + 1), 
                     .packages = c("mboost", "survival")) %dopar%{
                       
                       boost <- glmboost(X[-i,], y[-i], 
                                         family = family,
                                         control = boost_control(mstop = mstops[i], 
                                                                 nu = nu, 
                                                                 risk = "inbag"), 
                                         center = F) 
                       return(names(coef(boost)))
                     }
  registerDoSEQ()
  stopCluster(cl)
  return(reboost)
}

# vs, variable selection for fixed weights
vsFixed <- mapply(drop1vs, 
                  datName = d, 
                  family = type, 
                  nu = nu,
                  MoreArgs = list(method = "fixed",
                                  nCores = 20), 
                  SIMPLIFY = F)

# vs, variable selection for adaptive weights
vsAda <- mapply(drop1vs, 
                datName = d, 
                family = type, 
                nu = nu,
                MoreArgs = list(method = "adaptive",
                                nCores = 20), 
                SIMPLIFY = F)

# save(list = c("vsFixed", "vsAda"), file = "../output/vs.RData")

#

load("../output/vs.RData")

# Figure 1, ada_m1vsm2
jpeg("../results/ada_m1vsm2.jpeg", width = 8, height = 5, units = "in", 
     res = 300)
mat <- matrix(1:6, nrow = 2, byrow = F)
layout(mat)
par(mar = c(3.5, 3.5, 1.5, 1.5), oma = c(1, 1, 1, 1))
for(s in d[c(3, 4, 6)]){
  
  obj <- vsAda[[s]]
  
  plot_vs(obj,  method = "M2",
          name = paste0(gsub(".RData", "", s), "(adaptive)"),
          xlim = c(0, 12), plot = "lol", center = F)
  
  plot_m1vsm2(obj, plot.type = "c", name = "")
}

dev.off()

# ada_m1vsm2 in Appendix
jpeg("../results/ada_m1vsm2_Appendix.jpeg", width = 8, height = 5, 
     units = "in", res = 300)
mat <- matrix(1:6, nrow = 2, byrow = F)
layout(mat)
par(mar = c(3.5, 3.5, 1.5, 1.5), oma = c(1, 1, 1, 1))
for(s in d[-c(3, 4, 6)]){
  
  obj <- vsAda[[s]]
  
  plot_vs(obj,  method = "M2",
          name = paste0(gsub(".RData", "", s), "(adaptive)"),
          xlim = c(0, 12), plot = "lol", center = F)
  
  plot_m1vsm2(obj, plot.type = "c", name = "")
}

dev.off()


# Figure 4

vsPathlist <- list()

for(i in c(3, 4, 6)){
  
  print(d[i])
  load(paste0("../output/mstop_", d[i]))
  load(paste0("../data/", d[i]))
  X <- t(X)
  
  family = type[[i]]
  mstops = unlist(lapply(mstops$mstops, `[[`, "mstop"))
  
  if (class(family) == "boost_family_glm")
    y = as.factor(y)
  if (class(family) == "boost_family")
    y = surv
  
  vsPathlist[[d[i]]] <- vs_path(X, y, family, obs.mstop = mstops, nu[i], 
                                nCores = 4)
}

# save(list = "vsPathlist", file = "../output/vsPathlist.RData")
load("../output/vsPathlist.RData")
load("../output/vs.RData")

jpeg("../results/vsPath.jpeg", width = 8, height = 5, units = "in", res = 300)
par(mfcol = c(2, 3), mar = c(3.5, 3.5, 2, 2))

for(i in names(vsPathlist)){
  
  # pathlist
  plot.vs_path(vsPathlist[[i]])
  mtext(gsub(".RData", "", i), side = 3, line = 0.2, adj = 0)
  
  # corplot to compare fixed and adaptive
  fixObj <- vsFixed[[i]]
  adaObj <- vsAda[[i]]
  
  plot_fixvsada(fixObj, adaObj, method = "M2", name = gsub(".RData", "", i))
}

dev.off()

# #-------------------------------------------------------------------------------
# #- prediction
# #-------------------------------------------------------------------------------
# 
# #- loo, leave one out
# 
# cv.pred <- function(X, y, family, folds, loo = 0){
#   
#   n <- nrow(X)
#   res <- c()
#   for(i in seq_len(ncol(folds))){
#     
#     f <- folds[,i]
#     if (loo <= n)
#       f[loo] <- 0
#     
#     boost <- glmboost(X[f == 1, ], y[f == 1], 
#                       family = family,
#                       control = boost_control(mstop = no.mstop, 
#                                               nu = 0.1, 
#                                               risk = "inbag"), 
#                       center = F)
#     
#     out <-  folds[,i] == 0
#     yy <- y[out]
#     
#     re <- mapply(function(k){
#       p <- boost[k]$predict(newdata = X[out,])
#       Binomial()@risk(yy, p)}, 
#       k = 1:no.mstop)
#     
#     res <- rbind(res, re)
#   }
#   return(res)
# }
# 
# mstopSet <- c(818, 3230, 1272, 136, 262, 218)
# cv.mstops <- list()
# 
# for(i in 1:6){
#   
#   no.mstop = mstops[i]
#   print(d[i])
#   print(Sys.time())
#   
#   #- load data
#   load(paste0(dir, d[i]))
#   
#   X <- t(X)
#   X <- scale(X)
#   n = nrow(X)
#   family = type[[i]]
#   
#   if (link[i] == "binomial") 
#     y = as.factor(y)
#   if (link[i] == "cox")
#     y = surv
#   
#   set.seed(1)
#   nfolds = 10
#   foldid <- sample(rep(1:nfolds, length = nrow(X)), nrow(X)) 
#   cv <- sapply(1:max(foldid), function(x) as.numeric(foldid != x))
#   
#   ncores <- detectCores() - 1
#   cl <- makeCluster(ncores)
#   registerDoParallel(cl)
#   
#   re.loo <- foreach(k = 1:(n + 1), 
#                     .packages = c("mboost")) %dopar% {
#                       
#                       re <- cv.pred(X, y, family = family, folds = folds, loo = 200)
#                       colSums(re)
#                     }
#   
#   load("../output/cv.mstops.RData")
#   mstops <- cv.mstops$golub99.RData
#   
#   re <- c()
#   for(i in 1:n){
#     
#     tmp <- re.loo[[i]][m2[i]] - re.loo[[n + 1]][m2[n + 1]]
#     re  <- c(re, tmp)
#   }
# 
#   #- likelihood displacement
# 
#   re <- c()
#   for(i in seq_len(n + 1)){
#     re <- c(re, colSums(re.mstops[[i]]$cvboost)[m2[i]])
#   }
#   
#   
#   re - re[n + 1]
#   order(abs(re - re[n + 1]), decreasing = T)
#   
#   re <- c()
#   for(i in seq_len(n + 1)){
#     re <- c(re, colSums(re.mstops[[i]]$cvboost)[i + 1])
#   }
#   
#   re - re[n + 1]
#   cv.mstops[[d[i]]] <- list(mstops = re.mstops, foldid = foldid)
#   
#   registerDoSEQ()
#   stopCluster(cl)
# }
# 
# #-------------------------------------------------------------------------------
# 
# 
# load("../output/cv.mstops.RData")
# jpeg("./output/mstops.jpeg", width = 10, height = 3, units = "in",
#      res = 300)
# par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
# for(s in d[c(3, 4, 6)]){
#   
#   obj <- cv.mstops[[s]][["mstops"]]
#   drop1 <- obj[-length(obj)]
#   drop0 <- obj[length(obj)]
#   score <- abs(drop1 - drop0)
#   
#   x = 1:length(drop1)
#   y = drop1
#   plot(x, y, xlab = "", ylab = "", cex = 0.8)
#   abline(h = drop0)
#   
#   for(i in 1:length(x)){
#     lines(c(x[i], x[i]), c(drop0, y[i]), col = "grey")
#   }
#   points(x, y, xlab = "", ylab = "", pch = 19, cex = 0.8)
#   mtext(gsub(".RData", "", s), side = 3, line = 0.2, adj = 0)
#   mtext(side = 2, line = 0.2, at = drop0, text = "orig", cex = 0.8, las = 1)
#   ord <- order(score, decreasing = T)[1:5]
#   text(ord, drop1[ord], pos = 4, cex = 0.9, labels = ord)
#   mtext("mstops", side = 2, line = 2.2)
#   mtext("observations", side = 1, line = 2.2)
#   
# }
# dev.off()
# 
# #-------------------------------------------------------------------------------
# 
# jpeg("./output/mstops.A.jpeg", width = 10, height = 3, units = "in",
#      res = 300) #- A, for appendix
# 
# par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
# 
# for(s in d[c(1, 2, 5)]){
#   
#   obj <- cv.mstops[[s]][["mstops"]]
#   drop1 <- obj[-length(obj)]
#   drop0 <- obj[length(obj)]
#   score <- abs(drop1 - drop0)
#   
#   x = 1:length(drop1)
#   y = drop1
#   plot(x, y, xlab = "", ylab = "", cex = 0.9)
#   abline(h = drop0)
#   
#   for(i in 1:length(x)){
#     lines(c(x[i], x[i]), c(drop0, y[i]), col = "grey")
#   }
#   points(x, y, xlab = "", ylab = "", pch = 19, cex = 0.8)
#   mtext(gsub(".RData", "", s), side = 3, line = 0.2, adj = 0)
#   mtext(side = 2, line = 0.2, at = drop0, text = "orig", cex = 0.8, las = 1)
#   ord <- order(score, decreasing = T)[1:3]
#   text(ord, drop1[ord], pos = 4, cex = 0.8, labels = ord)
#   mtext("mstops", side = 2, line = 2.2)
#   mtext("observations", side = 1, line = 2.2)
#   
# }
# dev.off()

#-------------------------------------------------------------------------------
