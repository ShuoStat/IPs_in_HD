#-------------------------------------------------------------------------------
library(mboost)
library(doParallel)
source("../pure-transformation/transfo.R")
source("./regressionfun.cv.R")
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
#- Prediction
#-------------------------------------------------------------------------------
library(doParallel)
library(doRNG)
library(doSNOW)

mstops <- c(50, 818, 7904, 3660, 136, 262, 50, 218)
#- re.cv.pred <- list()

load("./output/re.cv.pred.RData")
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
  
  ncores <- detectCores() - 2
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(min = 1, max = (n + 1), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cv.pred <- foreach(i = 1:(n + 1), .packages = c("mboost"),
                     .options.snow = opts) %dopar% {
                       
                       boost <- glmboost(X[-i,], y[-i],
                                         family = family,
                                         control = boost_control(mstop = no.mstop,
                                                                 nu = 0.1,
                                                                 risk = "inbag"),
                                         center = T)

                       cvboost <- cvrisk(boost, folds = cv[-i,])
                       return(cvboost)
                     }
  
  re.cv.pred[[d[i]]] <- list(cv.pred = cv.pred, foldid = foldid)
  
  close(pb)
  stopCluster(cl)
}

#- save(list = "re.cv.pred", file = "./output/re.cv.pred.RData")

#-------------------------------------------------------------------------------
#- Vasulization of CV prediction
#-------------------------------------------------------------------------------

plot.cv.pred.path <- function(obj, top = 3){
  
  n <- length(obj)
  obj <- lapply(obj, colMeans)
  obj <- Reduce(cbind, obj)
  obj <- apply(obj, 2, function(x) abs(x - obj[,n]))[,-n]
  
  ylim = c(0, max(obj))
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

plot.cv.pred.score <- function(obj, top = 3, ylim = NULL){
  
  n <- length(obj)
  obj <- lapply(obj, colMeans)
  obj <- Reduce(cbind, obj)
  obj <- apply(obj, 2, function(x) abs(x - obj[,n]))[,-n]
  
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

load("./output/re.cv.pred.RData")

jpeg("./output/cvpath.jpeg", width = 9, height = 5, units = "in", res = 300)
layout(matrix(1:6, nrow = 2, byrow = F))
par(mar = c(3, 3, 2, 1), oma = c(1, 1, 1, 1))

s <- c("petricoin02.RData", "golub99.RData", "veer02.RData")

for(i in s){
  obj <- re.cv.pred[[i]]
  plot.cv.pred.path(obj$cv.pred)
  mtext(gsub(".RData", "", i), adj = 0, line = 0.2, side = 3)
  plot.cv.pred.score(obj$cv.pred, ylim = c(0, 13))
}
dev.off()

jpeg("./output/cvpath.append.jpeg", width = 9, height = 5, units = "in", 
     res = 300)
mat <- matrix(1:6, nrow = 2, byrow = F)
layout(mat)
par(mar = c(3, 3, 2, 1), oma = c(1, 1, 1, 1))
s <- c("singh02.RData", "wang05.RData", "ressom06.RData")
for(i in s){
  obj <- re.cv.pred[[i]]
  plot.cv.pred.path(obj$cv.pred)
  mtext(gsub(".RData", "", i), adj = 0, line = 0.2, side = 3)
  plot.cv.pred.score(obj$cv.pred, ylim = c(0, 13))
}
dev.off()

#-------------------------------------------------------------------------------
#- prediction, loglik
#-------------------------------------------------------------------------------
#- fixed mstops------

load("./output/cv.mstops.RData")

re.boost.loglik.fixed <- list()
for(i in 1:8){
  family = type[[i]]
  load(paste0(dir, d[i]))
  X <- t(X)

  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  cv = cv.mstops[[d[i]]][["mstops"]]
  mstop = cv[length(cv)]
  obs.mstop = cv[-length(cv)]

  re.boost.loglik.fixed[[d[i]]] <- boost.pred.loglik(X, y, 
                                                     mstop = mstop, 
                                                     obs.mstop = obs.mstop,
                                                     family = family,
                                                     method = "fixed",
                                                     ncores = NULL)
}
#- save(list = "re.boost.loglik.fixed", file = "./output/re.boost.loglik.fixed.RData")

#- cross validated mstops -#

re.boost.loglik.cv <- list()
for(i in 1:8){
  family = type[[i]]
  load(paste0(dir, d[i]))
  X <- t(X)

  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  cv = cv.mstops[[d[i]]][["mstops"]]
  mstop = cv[length(cv)]
  obs.mstop = cv[-length(cv)]

  re.boost.loglik.cv[[d[i]]] <- boost.pred.loglik(X, y, 
                                                  mstop = mstop, 
                                                  obs.mstop = obs.mstop,
                                                  family = family,
                                                  method = "cv",
                                                  ncores = NULL)
}

# save(list = "re.boost.loglik.cv", file = "./output/re.boost.loglik.cv.RData")


#-------------------------------------------------------------------------------
#- plot, loglik
#-------------------------------------------------------------------------------

load("./output/re.boost.loglik.cv.RData")
load("./output/re.boost.loglik.fixed.RData")

plot.loglik <- function(obj, xlim){
  drop1 <- obj$drop1
  drop0 <- obj$drop0
  
  score <- abs(drop0 - drop1)
  score <- score / sd(score)
  names(score) <- paste0("obs", 1:length(score))
  hor.barplot(score, cex.names = 0.8, xlim = c(0, xlim))
}

jpeg("./output/loglik.maintxt.jpeg", width = 9, height = 6, units = "in",
     res = 300)
layout(matrix(1:6, nrow = 2, byrow = F))
par(mar = c(0.5, 3.5, 1.5, 1.5))

for(s in d[c(4, 5, 8)]){
  
  plot.loglik(re.boost.loglik.fixed[[s]], xlim = 14)
  mtext(paste0(gsub(".RData","", s), "(fixed)"), 
        line = 0, side = 3, adj = 0)
  
  print(s)
  plot.loglik(re.boost.loglik.cv[[s]], xlim = 14)
  mtext(paste0(gsub(".RData","", s), "(adaptive)"),
        line = 0, side = 3, adj = 0)
}
dev.off()


jpeg("./output/loglik.append.jpeg", width = 15, height = 6, units = "in",
     res = 300)
layout(matrix(1:10, nrow = 2, byrow = F))
par(mar = c(0.5, 3.5, 1.5, 1.5))

for(s in d[-c(4, 5, 8)]){
  
  plot.loglik(re.boost.loglik.fixed[[s]], xlim = 10)
  mtext(paste0(gsub(".RData","", s), "(fixed)"), 
        line = 0, side = 3, adj = 0)
  
  plot.loglik(re.boost.loglik.cv[[s]], xlim = 10)
  mtext(paste0(gsub(".RData","", s), "(adaptive)"),
        line = 0, side = 3, adj = 0)
}
dev.off()

#-------------------------------------------------------------------------------
#- loglik vs boosting path
#-------------------------------------------------------------------------------

#- 1

load("./output/re.boost.loglik.cv.RData")
load("./output/re.boost.beta.cv.RData")

jpeg("./output/loglik.vs.beta.main.jpeg", width = 9, height = 3, units = "in",
     res = 300)
par(mar = c(3.5, 3.5, 2, 2), mfrow = c(1, 3))

for (s in d[c(4, 5, 8)]){

  #- loglik - 
  obj.loglik <- re.boost.loglik.cv[[s]]
  drop1 <- obj.loglik$drop1
  drop0 <- obj.loglik$drop0
  s.pred <- abs(drop0 - drop1)
  s.pred <- s.pred / sd(s.pred)
  names(s.pred) <- paste0("obs", 1:length(s.pred))
  
  #- betapath -
  obj.beta <- re.boost.beta.cv[[s]][["cvbeta"]]
  s.beta <- obj.beta / sd(obj.beta)
  names(s.beta) <- paste0("obs", 1:length(s.beta))
  
  lim = c(min(c(s.beta, s.pred)), max(s.beta, s.pred) * 1.1)
  corrplot(s.beta, s.pred, ci.line = F, xlim = lim,
           ylim = lim)
  o.beta <- order(s.beta, decreasing = T)[1:3]
  o.pred <- order(s.pred, decreasing = T)[1:3]
  text(s.beta[o.beta], s.pred[o.beta], labels = o.beta, cex = 0.6, pos = 4)
  text(s.beta[o.pred], s.pred[o.pred], labels = o.pred, cex = 0.6, pos = 4)
  
  mtext("coefficient changes", side = 1, line = 2.2)
  mtext("prediction changes", side = 2, line = 2.2)
  mtext(gsub(".RData", "", s), side = 3, line = 0.2, adj = 0)
}
dev.off()

#2-

jpeg("./output/loglik.vs.beta.append.jpeg", 
     width = 9, height = 6, units = "in", res = 300)
par(mar = c(3.5, 3.5, 2, 2), mfrow = c(2, 3))

for (s in d[-c(4, 5, 8)]){
  
  #- loglik - 
  obj.loglik <- re.boost.loglik.cv[[s]]
  drop1 <- obj.loglik$drop1
  drop0 <- obj.loglik$drop0
  s.pred <- abs(drop0 - drop1)
  s.pred <- s.pred / sd(s.pred)
  names(s.pred) <- paste0("obs", 1:length(s.pred))
  
  #- betapath -
  obj.beta <- re.boost.beta.cv[[s]][["cvbeta"]]
  s.beta <- obj.beta / sd(obj.beta)
  names(s.beta) <- paste0("obs", 1:length(s.beta))
  
  lim = c(min(c(s.beta, s.pred)), max(s.beta, s.pred) * 1.1)
  corrplot(s.beta, s.pred, ci.line = F, xlim = lim,
           ylim = lim)
  o.beta <- order(s.beta, decreasing = T)[1:3]
  o.pred <- order(s.pred, decreasing = T)[1:3]
  text(s.beta[o.beta], s.pred[o.beta], labels = o.beta, cex = 0.6, pos = 4)
  text(s.beta[o.pred], s.pred[o.pred], labels = o.pred, cex = 0.6, pos = 4)
  
  mtext("coefficient changes", side = 1, line = 2.2)
  mtext("prediction changes", side = 2, line = 2.2)
  mtext(gsub(".RData", "", s), side = 3, line = 0.2, adj = 0)
}
dev.off()

#-------------------------------------------------------------------------------
#- Explore the likelihood of logistic
#-------------------------------------------------------------------------------

data <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
n = nrow(data)
y = data$admit
loglik <- c()
for(i in 1:(n+1)){
  mod <- glm(admit ~ gre + gpa + rank, data = data[-i,], family = "binomial")
  p <- predict(mod, newdata = data, type = "response")
  loglik <- cbind(loglik, y * log(p) + (1 - y)*log(1 - p))
}

loglik.sum <- apply(loglik, 2, sum)
loglik.dis <- loglik.sum[n + 1] - loglik.sum 

#-------------------------------------------------------------------------------
#- cook's D 
#-------------------------------------------------------------------------------

#- fixed mstops------
load("./output/cv.mstops.RData")

re.boost.cook.fixed <- list()
for(i in 1:8){
  family = type[[i]]
  load(paste0(dir, d[i]))
  X <- t(X)

  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  cv = cv.mstops[[d[i]]][["mstops"]]
  mstop = cv[length(cv)]
  obs.mstop = cv[-length(cv)]
  
  re.boost.cook.fixed[[d[i]]] <- boost.pred.loglik(X, y, 
                                                   mstop = mstop, 
                                                   obs.mstop = obs.mstop,
                                                   family = family,
                                                   method = "fixed",
                                                   ncores = NULL,
                                                   measure = "cook")
}
#- save(list = "re.boost.cook.fixed", file = "./output/re.boost.cook.fixed.RData")

#- cross validated mstops

re.boost.cook.cv <- list()
for(i in 1:8){
  family = type[[i]]
  load(paste0(dir, d[i]))
  X <- t(X)

  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  cv = cv.mstops[[d[i]]][["mstops"]]
  mstop = cv[length(cv)]
  obs.mstop = cv[-length(cv)]
  
  re.boost.cook.cv[[d[i]]] <- boost.pred.loglik(X, y, 
                                                mstop = mstop, 
                                                obs.mstop = obs.mstop,
                                                family = family,
                                                method = "cv",
                                                ncores = NULL,
                                                measure = "cook")
}
# save(list = "re.boost.cook.cv", file = "./output/re.boost.cook.cv.RData")

#-------------------------------------------------------------------------------

load("./output/re.boost.cook.cv.RData")
load("./output/re.boost.cook.fixed.RData")

# - cook main text
jpeg("./output/boost.cook.main.jpeg", width = 9, height = 6, units = "in",
     res = 300)
par(mfrow = c(2, 3), mar = c(1, 4, 2, 2))
for (i in d[c(4, 5, 8)]){
  obj <- re.boost.cook.fixed[[i]]
  drop1 <- obj$drop1
  drop0 <- obj$drop0
  
  m <- abs(drop1 - drop0)
  score <- apply(m, 2, sum)
  score <- score / sd(score)
  names(score) <- paste0("obs", 1:length(score))
  
  hor.barplot(score, xlim = c(0, 13))
  mtext(paste0(gsub(".RData","", i), "(fixed)"), 
        line = 0, side = 3, adj = 0)
}

for (i in d[c(4, 5, 8)]){
  obj <- re.boost.cook.cv[[i]]
  drop1 <- obj$drop1
  drop0 <- obj$drop0
  
  m <- abs(drop1 - drop0)
  score <- apply(m, 2, sum)
  score <- score / sd(score)
  names(score) <- paste0("obs", 1:length(score))
  
  hor.barplot(score, xlim = c(0, 13))
  mtext(paste0(gsub(".RData","", i), "(adaptive)"), 
        line = 0, side = 3, adj = 0)
}

dev.off()

#- cook append 

jpeg("./output/boost.cook.append.jpeg", width = 15, height = 6, units = "in",
     res = 300)
par(mfrow = c(2, 5), mar = c(1, 4, 2, 2))
for (i in d[-c(4, 5, 8)]){
  obj <- re.boost.cook.fixed[[i]]
  drop1 <- obj$drop1
  drop0 <- obj$drop0
  
  m <- abs(drop1 - drop0)
  score <- apply(m, 2, sum)
  score <- score / sd(score)
  names(score) <- paste0("obs", 1:length(score))
  
  hor.barplot(score, xlim = c(0, 13))
  mtext(paste0(gsub(".RData","", i), "(fixed)"), 
        line = 0, side = 3, adj = 0)
}

for (i in d[-c(4, 5, 8)]){
  obj <- re.boost.cook.cv[[i]]
  drop1 <- obj$drop1
  drop0 <- obj$drop0
  
  m <- abs(drop1 - drop0)
  score <- apply(m, 2, sum)
  score <- score / sd(score)
  names(score) <- paste0("obs", 1:length(score))
  
  hor.barplot(score, xlim = c(0, 13))
  mtext(paste0(gsub(".RData","", i), "(adaptive)"), 
        line = 0, side = 3, adj = 0)
}
dev.off()

#-------------------------------------------------------------------------------
#- cook, fixed vs cv
#-------------------------------------------------------------------------------

#- 1
load("./output/re.boost.cook.cv.RData")
load("./output/re.boost.cook.fixed.RData")
load("./output/re.boost.beta.cv.RData")

jpeg("./output/cook.fix.vs.cv.jpeg", width = 9, height = 3, units = "in",
     res = 300)
par(mar = c(3.5, 3.5, 2, 2), mfrow = c(1, 3))

for (s in d[c(4, 5, 8)]){
  
  #- fixed 
  obj <- re.boost.cook.fixed[[s]]
  drop1 <- obj$drop1
  drop0 <- obj$drop0
  
  m <- abs(drop1 - drop0)
  score <- apply(m, 2, sum)
  c.fixed <- score / sd(score)

  #- cook'D - 
  obj <- re.boost.cook.cv[[s]]
  drop1 <- obj$drop1
  drop0 <- obj$drop0
  
  m <- abs(drop1 - drop0)
  score <- apply(m, 2, sum)
  c.cv <- score / sd(score)

  lim = c(min(c(c.fixed, c.cv)), max(c.fixed, c.cv) * 1.1)
  corrplot(c.fixed, c.cv, ci.line = F, xlim = lim,
           ylim = lim)
  o.fixed <- order(c.fixed, decreasing = T)[1:3]
  o.cv <- order(c.cv, decreasing = T)[1:3]
  text(c.fixed[o.fixed], c.cv[o.fixed], labels = o.fixed, cex = 0.6, pos = 4)
  text(c.fixed[o.cv], c.cv[o.cv], labels = o.cv, cex = 0.6, pos = 4)
  mtext("fixed", side = 1, line = 2.2)
  mtext("adaptive", side = 2, line = 2.2)
  mtext(gsub(".RData", "", s), side = 3, line = 0.2, adj = 0)
  
  r = format(round(cor(c.fixed, c.cv), 3), nsmall = 3)
  legend("topleft", legend = bquote(rho ~ " = " ~ .(r)), 
         bty = "n", inset = 0.02, cex = 1.5)
  
}
dev.off()


#-------------------------------------------------------------------------------
#- Lasso cv path
#-------------------------------------------------------------------------------

il.scale<-function(x) {
  tmp = scale(x)
  return (replace(tmp, is.nan(tmp), 0))
}

whichfolds <- function(nobs, nfolds) {
  
  howmany <- rep(floor(nobs/nfolds), nfolds)
  temp <- nobs%%nfolds
  if(temp>0) howmany[1:temp] <- howmany[1:temp] + 1
  return(sample(rep(1:nfolds,howmany), size = nobs, replace = FALSE))
}

lambda.int <- function(cvob, mult) {
  ind <- sort.list(cvob$cvm)[1]
  upper <- min(cvob$cvm) + mult * cvob$cvsd[ind]
  temp <- cvob$lambda[which(cvob$cvm <= upper)]
  return(c(min(temp), max(temp)))
}

internal.df.cvpath<-function(fit1, fit2) {
  #Change in cv-path between fit1 and fit2
  #fit1: fitted cv.glmnet
  #fit2: fitted cv.glmnet
  
  if (length(fit1$cvm) == length(fit2$cvm)) {
    x = rev(fit1$lambda)
    y = rev(abs(fit1$cvm - fit2$cvm))
    return(caTools::trapz(x, y))
  }
  else {
    comlam = intersect(fit1$lambda, fit2$lambda)
    idx1 = fit1$lambda%in%comlam
    idx2 = fit2$lambda%in%comlam        
    x = rev(fit1$lambda[idx1])
    y = rev(abs(fit1$cvm[idx1]-abs(fit2$cvm[idx2])))
    out = caTools::trapz(x,y)
    return(out)
  }
}

df.cvpath <- function(X, y, obs = 1:nrow(X), family = "gaussian",
                      mult = 2, assign = NULL, nfolds = NULL, ncors = NULL) {
  X = scale(X)
  n = length(y)
  res <- rep(NA, length(obs))
  
  if (is.null(assign)) {
    if (is.null(nfolds)) {
      assign = whichfolds(n, min(10,n))
    } else {
      assign = whichfolds(n, nfolds)
    }
  }
  
  fitfull <- glmnet::cv.glmnet(X, y, family = family,
                               standardize = TRUE, 
                               foldid = assign)
  
  lamint <- lambda.int(cvob = fitfull, mult = mult)
  lamint <- exp(seq(log(lamint[1]), log(lamint[2]), length.out = 100))
  fitfull <- glmnet::cv.glmnet(X, y, family = family,
                               standardize = TRUE, 
                               foldid = assign, 
                               lambda = lamint)
  
  
  require(doParallel)
  if(is.null(ncors))
    ncors <- detectCores() - 1
  cl <- makeCluster(ncors)
  registerDoParallel(cl)
  
  cvm.mat <- foreach(i = obs, .packages = "glmnet",
                     .combine = "rbind") %dopar% {
                       fitred = glmnet::cv.glmnet(X[-i,], y[-i], family = family,
                                                  standardize = TRUE,
                                                  foldid = assign[-i], 
                                                  lambda = lamint)
                       
                       abs(fitfull$cvm - fitred$cvm)
                     }
  stopCluster(cl)
  
  x = rev(fitfull$lambda)
  for(i in 1:nrow(cvm.mat)){
    
    y = rev(as.vector(cvm.mat[i,]))
    res[i] <- caTools::trapz(x, y)
  }
  
  return(list(score = res,
              cvm = cvm.mat,
              lambda = fitfull$lambda))
}


plot.cvpath <- function(obj, topn = 1){
  
  x <- obj$lambda
  y <- obj$cvm
  n = ncol(y)
  plot.new()
  plot.window(ylim = c(min(y), max(y)), xlim = c(min(x), max(x)))
  axis(1)
  axis(2)
  box()
  apply(y, 1, function(yy){
    lines(x, yy, col = "grey")
  })
  
  score <- obj$score
  ind <- order(score, decreasing = T)[1:topn]
  
  apply(y[ind,,drop = F], 1, function(yy){
    lines(x, yy, col = "black")
  })
  
  text(x = rep(min(x), topn), y = y[ind, n], labels = ind, cex = 0.9)
  
  mtext("Lambda", side = 1, line = 2.2)
  mtext(bquote("|" ~ Delta ~ "cve" ~ "|"), side = 2, line = 2.2)
}

plot.cvpath.score <- function(obj, center = T){
  
  score <- scale(obj$score, center)
  x = 1:length(score)
  y = score
  ylim = c(min(y) * center * 1.2 , max(y) * 1.2) 
  plot(x, y, xlab = "", ylab = "", cex = 0.8, yaxs = "i", 
       ylim = ylim)
  
  abline(h = 0)
  abline(h = 2, lty = "dashed")
  for(i in seq_along(x)){
    lines(c(x[i], x[i]), c(0, y[i]), col = "grey")
  }
  points(x, y, xlab = "", ylab = "", pch = 19, cex = 0.8)
  mtext("Df-Cvpath", side = 2, line = 2.2)
  mtext("observations", side = 1, line = 2.2)
  ord <- order(abs(score), decreasing = T)[1:5]
  text(ord, y[ord], pos = 4, cex = 0.8, labels = x[ord])
}  
  
#- seed 1, 135846
#- seed 2, 956723

set.seed(956723)
assign <- list()

for(i in 1:8){
  
  load(paste0(dir, d[i]))
  n <- ncol(X)
  
  assign[[i]] <- sample(rep(1:10, length = n), n)
}
  
re.lasso.cv.path <- list()
for(i in 1:8){
  
  load(paste0(dir, d[i]))
  X <- t(X)
  
  if (link[i] == "cox")
    y = surv
  
  re.lasso.cv.path[[d[i]]] <- df.cvpath(X, y, 
                                        obs = 1:nrow(X), 
                                        family = link[i],
                                        mult = 2,
                                        assign = assign[[i]],
                                        nfolds = NULL)
}

# save(list = "re.lasso.cv.path", file = "./output/lasso.cv.path.RData")

load("./output/lasso.cv.path.RData")

jpeg("./output/cvpath.jpeg", width = 10, height = 10, units = "in", res = 300)
mat1 <- matrix(1:8, nrow = 2, byrow = F)
mat2 <- matrix(9:16, nrow = 2, byrow = F)

layout(rbind(mat1, mat2))

for(i in d){
  
  par(mar = c(3.5, 3.5, 1.5, 1.5))
  obj <- re.lasso.cv.path[[i]]
  plot.cvpath(obj, 3)
  mtext(gsub(".RData", "", i), side = 3, line = 0.2, cex = 1.1)
  
  lollipop(scale(obj$score), decreasing = T, topn = 5, ylim = c(-2, 10), 
           refline = c(-2, 2), ylab = "Df-Cvpath")
  mtext(gsub(".RData", "", i), side = 3, line = 0.2, cex = 1.1)
}

dev.off()


#-------------------------------------------------------------------------------
#- tmp

plot.cvpath(obj, 3)
axis(1, at = obj$lambda, tick = T, labels = F)

#-

load("./output/lasso.cv.path.RData")

lam <- re.lasso.cv.path$golub99.RData$lambda

i = 5
load(paste0(dir, d[i]))
X <- t(X)

jpeg("./output/tmp.jpeg", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))

n = nrow(X)
assign = sample(rep(1:10, length = n), n)

rm <- c(100, 1, 2, 17)
nam <- c("orig", 1, 2, 17)
for(i in 1:4){
  fitfull <- glmnet::cv.glmnet(X[-rm[i],], y[-rm[i]], 
                               family = "binomial", 
                               standardize = TRUE, 
                               foldid = assign[-rm[i]])
  plot(fitfull)
  legend("topleft", legend = paste0("(-", nam[i], ")"))
  abline(v = log(c(min(lam), max(lam))))
}

dev.off()

#-------------------------------------------------------------------------------
#- tmp- calulation of LD
#-------------------------------------------------------------------------------

boost.pred <- function(X, y, family, 
                       mstop = NULL, 
                       obs.mstop = NULL, 
                       method = c("fixed", "cv"),
                       ncores = NULL,
                       measure = c("loglik", "cook"),
                       output = c("p", "loglik")){
  
  n <- nrow(X)
  p <- ncol(X)
  
  #- (n + 1) is the boosting without removing obs
  if (is.null(ncores)) ncores <- detectCores() - 1
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  if (method == "fixed") mstops = rep(mstop, (n + 1))
  if (method == "cv")    mstops = c(obs.mstop, mstop)
  
  reboost <- foreach(i = 1:(n + 1), j = mstops,
                     .combine = cbind,
                     .packages = c("mboost", "survival"),
                     .export = "loglik.coxph") %dopar% {
                       
                       boost <- glmboost(X[-i,], y[-i], 
                                         family = family,
                                         control = boost_control(mstop = j, 
                                                                 nu = 0.1, 
                                                                 risk = "inbag"), 
                                         center = T)
                       
                       
                       if (isTRUE(all.equal(family, Binomial()))){
                         y0   <- as.numeric(y) - 1
                         pred <- predict(boost, newdata = X, 
                                         type = "response")
                         
                         if (output == "loglik")
                           pred <- y0 * log(pred) + (1 - y0)*log(1 - pred)
                       }
                       
                       if (isTRUE(all.equal(family, CoxPH()))){
                         pred   <- predict(boost, newdata = X)
                         if (output == "loglik")
                           pred <- loglik.coxph(pred, 
                                                surv.time  = y[,1], 
                                                surv.event = y[,2])	
                       }
                       return(pred)
                     }
  
  registerDoSEQ()
  stopCluster(cl)
  
  drop1 = reboost[,-(n + 1)]
  drop0 = reboost[,(n + 1)]
  
  out <- list(drop1 = drop1, drop0 = drop0, y = y, family = family, 
              mstop = mstops, method = method, measure = measure, 
              output = output)
  return(out)
}

#-------------------------------------------------------------------------------
load("./output/cv.mstops.RData")

re.boost.loglik.fixed <- list()
i = 5
family = type[[i]]
load(paste0(dir, d[i]))
X <- t(X)

y <- as.factor(y)

cv = cv.mstops[[d[i]]][["mstops"]]
mstop = cv[length(cv)]
obs.mstop = cv[-length(cv)]

out.p <- boost.pred(X, y, 
                    mstop = mstop, 
                    obs.mstop = obs.mstop,
                    family = family,
                    method = "cv",
                    measure = "loglik",
                    ncores = NULL,
                    output = "loglik")

s1 <- colSums(out.p$drop1)
s0 <- sum(out.p$drop0)
s <- -2 * (s1 - s0)
n <- length(s)

jpeg("./output/LD.jpeg", width = 6, height = 4, units = "in", res = 300)
par(mar = c(3.5, 3.5, 1, 1))
plot.new()
plot.window(xlim = c(1, length(s)), ylim = c(min(s), max(s)))
axis(1)
axis(2)
box()

for(i in 1:n){
  lines(c(i, i), c(0, s[i]), col = "grey")
}

sel <- s > 0
points((1:n)[sel], s[sel], pch = 19)
points((1:n)[!sel], s[!sel], pch = 19, col = "red")

abline(h = 0)
mtext("LD", side = 2, line = 2.2)
mtext("Observations", side = 1, line = 2.2)
dev.off()

#- willi's measure

p.drop1 <- diag(out.p$drop1)
p.orig  <- out.p$drop0
y0 <- as.numeric(y) - 1
dif <- abs(p.drop1 - p.orig)

dif.min <- pmin(abs(p.drop1 - y0), abs(p.orig - y0))
score <- dif / dif.min
boxplot(split(score, y0))


dat <- data.frame(paste0("obs", 1:72), y = y0, p.orig = p.orig,
                  p.drop1 = p.drop1, a = abs(p.drop1 - y0),
                  b = abs(p.orig - y0), min = dif.min, score = score)

write.csv(dat, "./output/dat.csv")

#-------------------------------------------------------------------------------

out.lik <- boost.pred(X, y, 
                      mstop = mstop, 
                      obs.mstop = obs.mstop,
                      family = family,
                      method = "fixed",
                      ncores = NULL,
                      output = "loglik")

dev.drop1 <- diag(out.lik$drop1) * (-2)
dev.orig  <- out.lik$drop0  * (-2)
y0 <- as.numeric(y) - 1
dif <- abs(dev.drop1 - dev.orig)
names(dif) <- paste0("obs", 1:72)

dat <- data.frame(paste0("obs", 1:72), dev.orig = dev.orig,
                  dev.drop1 = dev.drop1, dif = dif)

library(PlotForMFP)

jpeg("./output/tmp.jpeg", width = 10, height = 5, units = "in",
     res = 300)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 2))
corrplot(dev.orig, dev.drop1, xlim = c(0, 10), xlab = "orig",
         ylab = "L-1-O", ci.line = F)

blandartman(dev.orig, dev.drop1, ci.limit = F, ci.spline = F,
            xlab = "(orig +L-1-O) / 2", ylab = "orig - L-1-O", xlim = c(0, 10),
            ylim = c(-10, 10))

dev.off()

#- -obs17
i = 17
out.lik <- boost.pred(X[-i,], y[-i], 
                      mstop = mstop, 
                      obs.mstop = obs.mstop[-i],
                      family = family,
                      method = "fixed",
                      ncores = NULL,
                      output = "loglik")

dev.drop1 <- diag(out.lik$drop1) * (-2)
dev.orig  <- out.lik$drop0  * (-2)
y0 <- as.numeric(y) - 1
dif <- abs(dev.drop1 - dev.orig)
names(dif) <- paste0("obs", (1:72)[-i])

dat <- data.frame(paste0("obs", (1:72)[-i]), dev.orig = dev.orig,
                  dev.drop1 = dev.drop1, dif = dif)

write.csv(dat, "./output/dat.csv")

jpeg("./output/tmp.jpeg", width = 5, height = 4, res = 300, units = "in")
par(mar = c(4, 4, 2, 2))
diff <- split(dif, y0[-i])
boxplot(diff)

mapply(function(i){
  x_i   <- diff[[i]]
  iqr <- quantile(x_i, 0.75) - quantile(x_i, 0.25)
  out <- which(x_i < (quantile(x_i, 0.25) - 1.5 * iqr) | 
                 x_i > (quantile(x_i, 0.75) + 1.5 * iqr))
  if (length(out)  > 0)
    text(rep(i, length(out)), x_i[out], labels = names(x_i)[out], 
         pos = 4, cex = 0.7)
}, i = 1:2)

dev.off()

jpeg("./output/tmp.jpeg", width = 5, height = 4, res = 300, units = "in")
par(mar = c(4, 4, 2, 2))
x = 1:length(dif)
y = dif
center = F
ylim = c(min(y) * center * 1.2 , max(y) * 1.2) 
plot(x, y, xlab = "", ylab = "", cex = 0.8, yaxs = "i", 
     ylim = ylim)

for(i in seq_along(x)){
  lines(c(x[i], x[i]), c(0, y[i]), col = "grey")
}
points(x, y, xlab = "", ylab = "", pch = 19, cex = 0.8)
mtext("DIF in deviance", side = 2, line = 2.2)
mtext("observations", side = 1, line = 2.2)
ord <- order(abs(dif), decreasing = T)[1:5]
text(ord, y[ord], pos = 4, cex = 0.8, labels = x[ord])
dev.off()

#- correlation plot
p.drop1 <- out.p$drop1
p.drop1 <- diag(out)

p.orig <- out.p$drop0

plot(p.drop1, p.orig)
text(p.drop1, p.orig, labels = 1:72, pos = 4, cex = 0.6)

#---
ld <- abs(colSums(tab.lik) - colSums(tab.lik)[1])[-1] * 2
ld <- sort(ld, decreasing = T)
ld.sd <- ld / sd(ld)
ld.sd

re <- round(cbind(ld[1:10], ld.sd[1:10]), 2)
write.csv(as.data.frame(re), "./output/tmp.csv")

#- explanatory plot for the cv-path

lambda.int <- function(cvob, mult) {
  ind <- sort.list(cvob$cvm)[1]
  upper <- min(cvob$cvm) + mult * cvob$cvsd[ind]
  temp <- cvob$lambda[which(cvob$cvm <= upper)]
  return(c(min(temp), max(temp)))
}


library(glmnet)
i = 5
family = "binomial"
load(paste0(dir, d[i]))
X <- t(X)

set.seed(956723)
assign <- list()

for(i in 1:8){
  
  load(paste0(dir, d[i]))
  n <- ncol(X)
  
  assign[[i]] <- sample(rep(1:10, length = n), n)
}

i = 5
family = "binomial"
load(paste0(dir, d[i]))
X <- t(X)

foldid <- assign[[5]]
cv.las <- cv.glmnet(X, y, family = family, foldid = foldid)

lam <- lambda.int(cv.las, 2)
lamint <- exp(seq(log(lam[1]), log(lam[2]), length.out = 100))

orig <- cv.glmnet(X, y, family = family, foldid = foldid, lambda = lamint)

plot.i <- function(i, col){
  drop1 <- cv.glmnet(X[-i, ], y[-i], family = family, 
                     foldid = foldid[-i], lambda = lamint)
  
  x = orig$lambda
  y1 = orig$cvm
  y2 = drop1$cvm
  
  plot.new()
  plot.window(xlim = c(min(x), max(x)), ylim = c(0, 1))
  axis(1)
  axis(2)
  box()
  
  mtext("CVM", side = 2, line = 2.2)
  mtext("lambda", side = 1, line = 2.2)
  mtext(paste0("Drop", i), side = 3, line = 0.2, adj = 0)
  
  polygon(x = c(x, rev(x)),
          y = c(y1, rev(y2)),
          col = "grey", border = NA)
  
  lines(x, y1, col = "black")
  lines(x, y2, col =  col)
  legend("topleft", legend = c("orig", paste0("drop", i)),
         col = c("black", col), lty = "solid", inset = 0.04)
}

jpeg("./output/las.jpeg", width = 10, height = 6, units = "in", res = 300)
par(mfrow = c(2, 3), mar = c(4, 4, 2, 2))

a = 17
drop1.a <- cv.glmnet(X[-a, ], y[-a], family = family, 
                     foldid = foldid[-a], lambda = lamint)

b = 69
drop1.b <- cv.glmnet(X[-b, ], y[-b], family = family, 
                     foldid = foldid[-b], lambda = lamint)
 
c = 55
drop1.c <- cv.glmnet(X[-c, ], y[-c], family = family, 
                     foldid = foldid[-c], lambda = lamint)
#- get usual case
load("./output/lasso.cv.path.RData")
obj <- re.lasso.cv.path[[i]]
case <- order(obj$score, decreasing = F)[72 / 2 + 1]
d = case
drop1.d <- cv.glmnet(X[-d, ], y[-d], family = family, 
                     foldid = foldid[-d], lambda = lamint)

x = orig$lambda
y1 = abs(orig$cvm - drop1.a$cvm)

y2 = abs(orig$cvm - drop1.b$cvm)

y3 = abs(orig$cvm - drop1.c$cvm)

y4 = abs(orig$cvm - drop1.d$cvm)

col = c("#E495A5", "#ABB065", "#39BEB1", "#ACA4E2")

plot.i(a, col[1])
plot.i(b, col[2])
plot.i(c, col[3])
plot.i(d, col[4])

plot(x, y1, ylim = c(0, 0.2), type = "l", xlab = "", ylab = "", col = col[1])
lines(x, y2, col = col[2])
lines(x, y3, col = col[3])
lines(x, y4, col = col[4])

mtext("lambda", side = 1, line = 2.2)
mtext("|CVM(Orig) - CVM(drop1)|", side = 2, line = 2.2)

legend("topright", inset = 0.03, col = col,
       lty = "solid", legend = paste0("obs", c(a, b, c, d)))
dev.off()

#-------------------------------------------------------------------------------

y1 <- as.numeric(y) - 1
drop1 <- out.p$drop1
drop0 <- out.p$drop0

jpeg("./output/tmp.jpeg", width = 8, height = 6, units= "in", res = 300)
par(mfrow = c(3, 4), mar = c(3, 3, 2, 2))
for(i in c(1:10, 17, 69)){
  
  dif0 <- abs(y1 - drop0)
  dif1 <- abs(y1 - drop1[,i])
  pmin <- pmin(dif0, dif1)
  trans <- (dif0 - dif1) / pmin
  names(trans) = paste0("obs", 1:72)
  
  dat <- cbind(y1, drop0, obs17  = drop1[,i], 
               dif0 = dif0, dif1 = dif1, pmin, trans)
  
  dat[,-1] <- round(dat[,-1], 2)
  
  dat <- cbind(paste0("obs", 1:72), dat)
  write.csv(dat, "./output/tmp.csv")
  
  x <- split(trans, y1)
  boxplot(x, name = "-obs17")
  mtext(paste0("-obs", i), side = 3, line = 0.2)
  
  mapply(function(i){
    x_i   <- x[[i]]
    
    iqr <- quantile(x_i, 0.75) - quantile(x_i, 0.25)
    out <- which(x_i < (quantile(x_i, 0.25) - 1.5 * iqr) | 
                   x_i > (quantile(x_i, 0.75) + 1.5 * iqr))
    if (length(out)  > 0)
      text(rep(i, length(out)), x_i[out], labels = names(x_i)[out], 
           pos = 4, cex = 0.7)
  }, i = 1:2)
}
dev.off()

#-------------------------------------------------------------------------------
#-tmp, explore loglik

load("./output/cv.mstops.RData")
i = 5
family = type[[i]]
load(paste0(dir, d[i]))
X <- t(X)

y <- as.factor(y)

cv = cv.mstops[[d[i]]][["mstops"]]
mstop = max(cv)
obs.mstop = cv[-length(cv)]

re <- boost.pred.loglik(X, y, 
                        mstop = mstop, 
                        obs.mstop = obs.mstop,
                        family = family,
                        method = "cv",
                        ncores = NULL,
                        path = T,
                        measure = "p")

#---

drop1 <- re$drop1
drop0 <- re$drop0[cv[length(cv)]]
ld <- (drop1 - drop0) * 2

plot.new()
n = nrow(ld)
plot.window(xlim = c(0, n), ylim = c(0, max(ld)))
axis(1)
axis(2)
box()

for(i in 1:ncol(ld)){
  lines(1:n, ld[,i], col = "grey")
}

for(i in 1:72){
  
  points(cv[i], ld[cv[i],i], pch = 19)

}












