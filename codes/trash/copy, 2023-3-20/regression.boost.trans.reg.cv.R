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
#- regression models 
#-------------------------------------------------------------------------------

transregress <- function(X, y, mstop = NULL, family, cv, adaptive = T, ... ){
  
  #- get cv folds
  require(mboost)
  X <- scale(X)
  Z <- transdata(X, byrow = F, eps = 0.01, scale = F)
  Z <- scale(Z)
  
  modX <- glmboost(X, y, 
                   family = family,
                   control = boost_control(mstop = mstop, 
                                           nu = 0.1, 
                                           risk = "inbag"), 
                   center = T)
  
  cv.x <- cvrisk(modX, folds = cv)
  cv.x <- mstop(cv.x)
  coefmodX <- unique(modX$xselect()[1:cv.x])
  
  modZ <- glmboost(Z, y, 
                   family = family,
                   control = boost_control(mstop = mstop, 
                                           nu = 0.1, 
                                           risk = "inbag"), 
                   center = T) 
  
  if (adaptive){
    cv.z <- cvrisk(modZ, folds = cv)
    cv.z <- mstop(cv.z)
  }else{
    cv.z = cv.x 
  }
  coefmodZ <- unique(modZ$xselect()[1:cv.z])
  return(list(X = coefmodX, Z = coefmodZ, cv.x = cv.x, cv.z = cv.z))
}


mstops <- c(50, 818, 7904, 3660, 136, 262, 50, 218)
trans.mod.cv <- list()
load("./output/cvs.RData")
for(i in 1:8){
  
  print(i)
  load(paste0(dir, d[i]))
  X <- t(X)
  family = type[[i]]
  
  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv 
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  trans.mod.cv[[d[i]]] <- transregress(X, y, mstop = mstops[i], 
                                       family = family, cv = cvs[[i]])
}

#- save(list = "trans.mod.cv", file = "./output/trans.mod.cv.RData")

#- fixed tuning parameters

mstops <- c(50, 818, 7904, 3660, 136, 262, 50, 218)
trans.mod.fixed <- list()
load("./output/cvs.RData")
for(i in 1:8){
  
  print(i)
  load(paste0(dir, d[i]))
  X <- t(X)
  family = type[[i]]
  
  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv 
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  trans.mod.fixed[[d[i]]] <- transregress(X, y, mstop = mstops[i], 
                                       family = family, cv = cvs[[i]],
                                       adaptive = F)
}

#- save(list = "trans.mod.fixed", file = "./output/trans.mod.fixed.RData")

#- table -----------------------------------------------------------------------

getset <- function(obj){
  
  if("trans" %in% class(obj)) 
    stop("need the results from transregress")
  both  <- intersect(obj$X, obj$Z)
  onlyX <- setdiff(obj$X, obj$Z)
  onlyZ <- setdiff(obj$Z, obj$X)
  both  <- length(both)
  onlyX <- length(onlyX)
  onlyZ <- length(onlyZ)
  
  all <- length(unique(c(obj$X, obj$Z)))
  olp <-  round(both / all, 2)
  
  out   <-  c(onlyX, onlyZ, both, all, olp)
  names(out) <- c("onlyX", "onlyZ", "both", "all", "overlap")
  return(out)
}

#- adaptive
load("./output/trans.mod.cv.RData")
tab <- lapply(trans.mod.cv, getset)
tab <- rbind.data.frame(tab)
write.csv(t(tab), file = "./output/tab.csv")

#- fixed
load("./output/trans.mod.fixed.RData")
tab <- lapply(trans.mod.fixed, getset)
tab <- rbind.data.frame(tab)
write.csv(t(tab), file = "./output/tab.csv")

#-------------------------------------------------------------------------------
#- boosting iteration path
#-------------------------------------------------------------------------------

trans.path.plot <- function(X, y, name, mstop.x, mstop.z, 
                            family,
                            abs = T, 
                            boxplot = F){
  
  colnames(X) <- paste0("V", 1:ncol(X))
  Z <- transdata(X, byrow = F, eps = 0.01, scale = F)
  X <- scale(X)
  Z <- scale(Z)
  
  modX <- glmboost(X, y, 
                   family = family,
                   control = boost_control(mstop = mstop.x, 
                                           nu = 0.1, 
                                           risk = "inbag"), 
                   center = T)
  
  modZ <- glmboost(Z, y, 
                   family = family,
                   control = boost_control(mstop = mstop.z, 
                                           nu = 0.1, 
                                           risk = "inbag"), 
                   center = T)
  
  coefmodX <- names(coef(modX))
  coefmodZ <- names(coef(modZ))
  
  int  <- intersect(coefmodX, coefmodZ)
  difX <- setdiff(coefmodX, coefmodZ)
  difZ <- setdiff(coefmodZ, coefmodX)
  
  if (length(difX) == 0 & length(difZ) == 0)
    txt = names(which.max(abs(c(coef(modX), coef(modZ)))))
  
  if(abs)
    ylim = c(0, max(abs(c(coef(modX), coef(modZ))))) else
      ylim = c(min(c(coef(modX), coef(modZ))), max(c(coef(modX), coef(modZ))))
  
  plot.boostpath(modX, group = list(int, difX), col = c("grey", "red"),
                 ylim = ylim, abs = abs)
  sname <- gsub(".RData", "", name)
  mtext(paste0(sname, " (original)"), side = 3, line = 0.2, adj = 0)
  #- add the maximum change
  if (max(abs(coef(modX)[difX])) > max(abs(coef(modZ)[difZ]))) {
    ind <- which.max(abs(coef(modX))[difX])
    txt <- difX[ind]
    xind <- mstop.x
    yind <- max(abs(coef(modX))[difX])
    text(xind, yind, labels = txt, pos = 2)
  }
  
  plot.boostpath(modZ, group = list(int, difZ), col = c("grey", "red"),
                 ylim = ylim, abs = abs)
  mtext(paste0(sname, " (transformed)"), side = 3, 
        line = 0.2, adj = 0)
  #- add the maximum change
  if (max(abs(coef(modX)[difX])) < max(abs(coef(modZ)[difZ]))) {
    ind <- which.max(abs(coef(modZ))[difZ])
    txt <- difZ[ind]
    xind <- mstop.z
    yind <- max(abs(coef(modZ))[difZ])
    text(xind, yind, labels = txt, pos = 2)
  }
  # 
  #- whether to add boxplots for the largest changes
  if (boxplot){
    rownames(X) <- rownames(Z) <- 1:nrow(X)
    #- txt in above steps is the most changed var
    if (isTRUE(all.equal(family, Binomial()))) {
      
      x <- append(split(X[,txt], y), split(Z[,txt], y))
      ylimax <- max(unlist(x)) + 0.2 * (max(unlist(x)) - min(unlist(x)))
      ylim = c(min(unlist(x)), ylimax)
      at = c(1, 2, 4, 5)
      boxplot(x, at = at, names = c("g0", "g1", "g0", "g1"),
              ylim = ylim)
      mtext(txt, side = 3, line = 0.2, adj = 0)
      
      out.txt <- list()
      for(i in seq_along(x)){
        
        x_i <- x[[i]]
        if (outliers::grubbs.test(x_i)$p.value < 0.001) {
          outs = 100000 #- the initial drops
          p.value = 0
          p.val = c()
          while(p.value < 0.001){
            p.value <- outliers::grubbs.test(x_i[!names(x_i) %in% outs])$p.value
            if (p.value < 0.001){
              outcase <- as.numeric(names(p.value))
              outs <- c(outs, outcase)
              p.val <- c(p.val, p.value)
            }
          }
          outs <- outs[-1] #- remove 100000
          at = c(1, 2, 4, 5)
          
          outs.001 <- outs[p.val < 0.001]
          points(at[rep(i, length(outs.001))], x_i[match(outs.001, names(x_i))], pch = 19, 
                 col = "red", cex = 0.9)
          
          out.txt[[i]] <- outs[1:2]
        } else {
          out.txt[[i]] <- NA
        }
        
        txt <- out.txt[[2 - (i %% 2)]]
        text(rep(at[i], length(txt)), 
             x_i[match(txt, names(x_i))], 
             labels = txt, pos = 4, cex = 0.9)
      }
      text(1.5, ylimax, labels = "Original", pos = 1, cex = 1.2)
      text(4.5, ylimax, labels = "Transformed", pos = 1, cex = 1.2)
    }
    
    if (isTRUE(all.equal(family, CoxPH()))) {
      
      x <- list(X[,txt], Z[,txt])
      ylimax <- max(unlist(x)) + 0.2 * (max(unlist(x)) - min(unlist(x)))
      ylim = c(min(unlist(x)), ylimax)
      boxplot(x, at = c(1, 2), names = c("orig.", "trans"),
              ylim = ylim)
      mtext(txt, side = 3, line = 0.2, adj = 0)
      
      out.txt <- list()
      for(i in seq_along(x)){
        x_i <- x[[i]]
        if (outliers::grubbs.test(x_i)$p.value < 0.001) {
          outs = 100000 #- the initial drops
          p.value = 0
          p.val = c()
          while(p.value < 0.001){
            p.value <- outliers::grubbs.test(x_i[!names(x_i) %in% outs])$p.value
            if (p.value < 0.001){
              outcase <- as.numeric(names(p.value))
              outs <- c(outs, outcase)
              p.val <- c(p.val, p.value)
            }
          }
          
          at = c(1, 2, 4, 5)
          outs <- outs[-1] #- remove 100000
          
          outs.001 <- outs[p.val < 0.001]
          points(at[rep(i, length(outs.001))], x_i[match(outs.001, names(x_i))], 
                 pch = 19, col = "red", cex = 0.9)
          out.txt[[i]] <- outs[1:2]
        } else {
          out.txt[[i]] <- NA
          at = c(1, 2, 4, 5)
        }
        
        txt <- out.txt[[1]]
        text(rep(at[i], length(txt)), x_i[match(txt, names(x_i))], 
             labels = txt, pos = 4, cex = 0.9)
      }
      text(1, ylimax, labels = "Original", pos = 1, cex = 1.2)
      text(2, ylimax, labels = "Transformed", pos = 1, cex = 1.2)
    }
  }
}

#-------------------------------------------------------------------------------

load("./output/trans.mod.cv.R")
jpeg(paste0("./output/trans.path.main.jpeg"), width = 9, height = 8, 
     units = "in", res = 300)
par(mfrow = c(3, 3), mar = c(4, 4, 2, 2))

for(i in c(2, 3, 6)){
  
  print(i)
  load(paste0(dir, d[i]))
  X <- t(X)
  family = type[[i]]
  
  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv 
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  nam <- gsub(d[i], ".RData", "")
  mods <- trans.mod.cv[[d[i]]]
  
  trans.path.plot(X, y, name = d[i],
                  mstop.x = mods$cv.x, 
                  mstop.z = mods$cv.x,
                  family = family, boxplot = T)
}

dev.off()

jpeg(paste0("./output/trans.path.append.jpeg"), width = 9, height = 8, 
     units = "in", res = 300)
par(mfrow = c(3, 3), mar = c(4, 4, 2, 2))

for(i in c(2, 3, 6)){
  
  print(i)
  load(paste0(dir, d[i]))
  X <- t(X)
  family = type[[i]]
  
  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv 
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  nam <- gsub(d[i], ".RData", "")
  mods <- trans.mod.cv[[d[i]]]
  
  trans.path.plot(X, y, name = d[i],
                  mstop.x = mods$cv.x, 
                  mstop.z = mods$cv.x,
                  family = family, boxplot = T)
}
dev.off()

#-------------------------------------------------------------------------------
#- Leverage
#-------------------------------------------------------------------------------

hatplot <- function(X, set, ylim = NULL, name, type = c("cor", "BA")){
  
  X <- X[, set, drop = F]
  X <- scale(X)
  n <- nrow(X)
  p <- ncol(X)
  Z <- transdata(X, byrow = F, eps = 0.01, scale = F)
  Z <- scale(Z)
  
  levx <- hat(X, intercept = F)
  levz <- hat(Z, intercept = F)
  n = length(levx)
  if(is.null(ylim)) 
    ylim = c(0.8 * min(c(levx, levz)), pmax(5 * p / n, 1.2 * max(c(levx, levz)))) 
  
  x = (levx + levz) / 2
  y = levx - levz
  intx <- (max(x) - min(x)) * 0.05
  
  s <- which(levx >= (3 * p / n))

  if (type == "cor"){
    corrplot(levx, levz, ci.line = F, xlab = "lev(OD)", ylab = "lev(TD)",
             xlim = ylim,
             ylim = ylim, line = F)
    abline(v = 3 * p / n, h = 5 * p / n, lty = "dashed")
    # abline(v = 3 * p / n, h = 3 * p / n, lty = "dashed")
    lines(c(-100, 100), c(-100, 100))
    if (length(s) > 0)
      text(levx[s], levz[s], labels = s, pos = 1, cex = 0.8)
  }
  
  if (type == "BA"){
    blandartman(levx, levz, ci = F, ci.limit = F, spline = F,
                xlab = 'lev(OD)/2 + lev(TD)/2',
                ylab = 'lev(OD) - lev(TD)', 
                xlim = c(min(x) - intx, max(x) + intx), 
                ylim = c(-max(y) * 1.05, max(y) * 1.05),
                withinlimits = F, LOA = F)
    
    sdy <- sd(y)
    # proportion beyond the limits
    z <- qnorm(0.99 + (1 - 0.99) / 2)
    abline(h = mean(y) - (z * sdy), lty = "longdash", lwd = 1)
    abline(h = mean(y) + (z * sdy), lty = "longdash", lwd = 1)
    
    if (length(s) > 0)
      text(x[s], y[s], labels = s, pos = 1, cex = 0.8)
  }
  mtext(name, side = 3, line = 0.2, adj = 0)
}

jpeg(paste0("./output/lev.bland.jpeg"), width = 9, height = 10, units = "in", 
     res = 300)
mat1 <- matrix(1:6, nrow = 2, ncol = 3, byrow = F)
mat2 <- matrix(7:12, nrow = 2, ncol = 3, byrow = F)

layout(rbind(mat1, mat2))
par(mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 1))

#- 1, 6
#- 2, 3
#- 4, 5
#- 7, 8
adaptive = F

for(i in c(8)){
  
  load(paste0(dir, d[i]))
  X <- t(X)
  
  if (adaptive){
    load("./output/trans.mod.cv.RData")
    sets  <- trans.mod.cv[[d[i]]]
  } else{
    load("./output/trans.mod.fixed.RData")
    sets  <- trans.mod.fixed[[d[i]]]
  }
  
  onlyX <- setdiff(sets$X, sets$Z)
  onlyZ <- setdiff(sets$Z, sets$X)
  both  <- intersect(sets$X, sets$Z)
  
  name = gsub(".RData", "", d[i])
  par(mar = c(3.5, 3.5, 1.5, 1.5))
  hatplot(X, onlyX, name = paste0("Only in OD(", name, ")"), 
          type = "cor")
  hatplot(X, onlyX, name = "Only in OD", type = "BA")
  hatplot(X, onlyZ, name = paste0("Only in TD(", name, ")"), 
          type = "cor")
  hatplot(X, onlyZ, name = "Only in TD", type = "BA")
  hatplot(X, both, name = paste0("In Both(", name, ")"),
          type = "cor")
  hatplot(X, both, name = "In Both", type = "BA")
}

dev.off()

#-------------------------------------------------------------------------------
#- Transformation on Prediction changes
#-------------------------------------------------------------------------------

trans.pred <- function(X, y, cv.x, cv.z, family, ... ){
  
  require(mboost)
  X <- scale(X)
  Z <- transdata(X, byrow = F, eps = 0.01, scale = F)
  Z <- scale(Z)
  
  modX <- glmboost(X, y, 
                   family = family,
                   control = boost_control(mstop = cv.x, 
                                           nu = 0.1, 
                                           risk = "inbag"), 
                   center = T) 
  
  modZ <- glmboost(Z, y, 
                   family = family,
                   control = boost_control(mstop = cv.z, 
                                           nu = 0.1, 
                                           risk = "inbag"), 
                   center = T) 
  
  pred.x <- predict(modX, X, type = "link")
  pred.z <- predict(modZ, Z, type = "link")

  return(list(X = pred.x, Z = pred.z))
}

labeltoblandartman <- function(x1, x2, alpha = 0.95){
  
  x <- (x1 + x2) / 2
  y <- x1 - x2
  sdy <- sd(y)
  stdy <- sdy / sqrt(length(y))
  
  z.alpha <- qnorm(alpha + (1 - alpha) / 2)
  out <- which(abs(y - mean(y)) > z.alpha * sdy)
  if (length(out) > 0 )
    text(x[out], y[out], labels = out, cex = 0.6, pos = 4)
}

# load("./output/trans.mod.cv.RData")
jpeg("./output/pred.trans.jpeg", width = 9, height = 5, units = "in", 
     res = 300)
layout(matrix(1:6, nrow = 2, byrow = F))
par(mar = c(3.5, 3.5, 1.5, 1.5))
adaptive = F

for(i in c(2, 3, 6)){
  
  load(paste0(dir, d[i]))
  family = type[[i]]
  X <- t(X)
  
  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv 
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  if (!adaptive){
    load("./output/trans.mod.fixed.RData")
    obj <- trans.mod.fixed[[d[i]]]
  }
  
  if (adaptive){
    load("./output/trans.mod.cv.RData")
    obj <- trans.mod.cv[[d[i]]]
  }
  
  #- fixed tuning   
  pred <- trans.pred(X, y, cv.x = obj$cv.x, cv.z = obj$cv.z, family)
  
  #- corrplot
  xlim <- ylim <- range(min(c(pred$X, pred$Z)), max(c(pred$X, pred$Z))) * 1.1
  corrplot(pred$X, pred$Z, xlab = "Pred(orig)", ylab = "Pred(trans)", 
           line = F, ci.line = F, xlim = xlim,
           ylim = ylim)
  lines(c(-100, 100), c(-100, 100), lwd = 1.5)
  mtext(gsub(".RData", "", d[i]), side = 3, line = 0.2, adj = 0)
  
  dif <- pred$X - pred$Z
  alpha = 0.99
  ylim <- c(pmin(min(dif) * 1.2, - sd(dif) * qnorm(alpha + (1 - alpha) / 2)),
            pmax(max(dif) * 1.2, sd(dif) * qnorm(alpha + (1 - alpha) / 2)))
  
  #- blandartman
  blandartman(as.numeric(pred$X), as.numeric(pred$Z), 
              xlab = "1/2 Pred(orig) + 1/2 Pred(trans)",
              ylab = "Pred(orig) - Pred(trans)",  
              ci.limit = FALSE, alpha = alpha, withinlimits = T,
              xlim = range(pred$X/2 + pred$Z/2) * 1.1, 
              ylim = ylim, 
              spline = F)
  
  labeltoblandartman(pred$X, pred$Z, alpha = alpha)
}

dev.off()

#-------------------------------------------------------------------------------

predtrans <- function(X, y, cv.x, cv.z, family, ... ){
  
  require(mboost)
  X <- scale(X)
  Z <- transdata(X, byrow = F, eps = 0.01, scale = F)
  Z <- scale(Z)
  
  #- 

  modX <- glmboost(X, y, 
                   family = family,
                   control = boost_control(mstop = cv.x, 
                                           nu = 0.1, 
                                           risk = "inbag"), 
                   center = T)

  modZ <- glmboost(Z, y, 
                   family = family,
                   control = boost_control(mstop = cv.z, 
                                           nu = 0.1, 
                                           risk = "inbag"), 
                   center = T)
  
  predx <- predict(modX, newdata = X, type = "response")
  predz <- predict(modZ, newdata = X, type = "response")
  
  out <- list(x = predx, z = predz)
  class(out) <- append(class(out), "predtrans")
  return(out)
}


jpeg("./output/transpred.jpeg", width = 9, height = 6, units = "in",
     res = 300)
par(mar = c(3.5, 3.5, 1.5, 1.5), mfrow = c(2, 3))
for (i in 1:6) {
  load(paste0(dir, d[i]))
  X <- t(X)
  y <- as.factor(y)
  re <- predtrans(X, y, cv.x, cv.z, family = Binomial())
  corrplot(re$x, re$z, p.col = "white")
  
  px <- split(re$x, y)
  py <- split(re$z, y)
  
  points(px[[1]], py[[1]], pch = 19, col = "steelblue4")
  points(px[[2]], py[[2]], pch = 19,  col = "tomato2")
  
  ind <- order(abs(re$x - re$z), decreasing = T)[1:5]
  text(re$x[ind], re$z[ind], labels = ind, cex = 0.7, pos = 4)
  mtext("Original model", side = 1, line =2.2)
  mtext("Transformed model", side = 2, line = 2.2)
  mtext(gsub(".RData", "", d[i]), side = 3, line = 0.2, adj = 0)
}
dev.off()

#-------------------------------------------------------------------------------

boxp <- function(vars, name){
  
  colnames(X) <- `if`(is.null(colnames(X)), paste0("V", 1:ncol(X)), colnames(X))
  x <- X[,vars]
  names(x) <- 1:length(x)
  x <- split(x, y)
  
  colnames(Z) <- `if`(is.null(colnames(Z)), paste0("V", 1:ncol(Z)), colnames(Z))
  z <- Z[,vars]
  names(z) <- 1:length(z)
  z <- split(z, y)
  ylim =  c(min(c(unlist(x), unlist(z))), 
            max(c(unlist(x), unlist(z))))
  
  boxplot(x, ylim = ylim)
  mapply(function(i){
    x_i   <- x[[i]]
    
    iqr <- quantile(x_i, 0.75) - quantile(x_i, 0.25)
    out <- which(x_i < (quantile(x_i, 0.25) - 1.5 * iqr) | 
                   x_i > (quantile(x_i, 0.75) + 1.5 * iqr))
    if (length(out)  > 0)
      text(rep(i, length(out)), x_i[out], labels = names(x_i)[out], 
           pos = 4, cex = 0.7)
  }, i = 1:2)
  mtext(paste0(name, "(original)"), side = 3, line = 0.2, adj = 0, cex = 0.8)
  
  boxplot(z, ylim = ylim)
  mapply(function(i){
    x_i <- z[[i]]
    iqr <- quantile(x_i, 0.75) - quantile(x_i, 0.25)
    out <- which(x_i < (quantile(x_i, 0.25) - 1.5 * iqr) | 
                   x_i > (quantile(x_i, 0.75) + 1.5 * iqr))
    if (length(out)  > 0)
      text(rep(i, length(out)), x_i[out], labels = names(x_i)[out], 
           pos = 4, cex = 0.7)
  }, i = 1:2)
  mtext(paste0(name, "(transformed)"), side = 3, line = 0.2, adj = 0, cex = 0.8)
}


#- explore scherzer07 -#
mstop = 100
s = "scherzer07.RData"
load(paste0(dir, s))
X <- t(X)
y <- as.factor(y)
Z <- transdata(X, byrow = F, eps = 0.01, scale = F)
X <- scale(X)
Z <- scale(Z)
modX <- glmboost(X, y, 
                 family = Binomial(),
                 control = boost_control(mstop = mstop, 
                                         nu = 0.1, 
                                         risk = "inbag"), 
                 center = T) 

modZ <- glmboost(Z, y, 
                 family = Binomial(),
                 control = boost_control(mstop = mstop, 
                                         nu = 0.1, 
                                         risk = "inbag"), 
                 center = T) 

betax <- coef(modX)
betaz <- coef(modZ)
onlyx <- setdiff(names(betax), names(betaz))
onlyz <- setdiff(names(betaz), names(betax))
beta <- betaz[onlyz][which.max(abs(betaz[onlyz]))]

jpeg("./output/box-trans.jpeg", width = 8, height = 6, units = "in",
     res = 300)
colname <- `if`(is.null(colnames(X)), paste0("V", 1:ncol(X)), colnames(X))
name = paste0("V", match(names(beta), colname))
par(mar = c(4, 4, 2, 1), mfrow = c(2, 2))
trans.path.plot(X, y, name = s, abs = T, family = Binomial())
boxp(names(beta), name)
dev.off()

#- the rest variables - 

png("./output/box1.png", width = 8, height = 10, units = "in", res = 300)
par(mfrow = c(5, 4), mar = c(3, 3, 1.5, 1))
mapply(boxp, vars = onlyx, name = onlyx)
dev.off()

png("./output/box2.png", width = 8, height = 6, units = "in", res = 300)
par(mfrow = c(3, 4), mar = c(3, 3, 1.5, 1))
mapply(boxp, vars = onlyz, name = onlyz)
dev.off()

#-------------------------------------------------------------------------------
#- boxplot for the variables of only in X, only in Z, and both
#-------------------------------------------------------------------------------

pbox <- function(vars, X, y, name, ylim1 = NULL, ylim2 = NULL){
  
  colnames(X) <- `if`(is.null(colnames(X)), paste0("V", 1:ncol(X)), colnames(X))
  rownames(X) <- 1:nrow(X)
  np <- length(vars)
  
  x1 <- X[y == levels(y)[1], vars]
  x2 <- X[y == levels(y)[2], vars]
  colnames(x1) <- colnames(x2) <- vars

  pbox <- function(x, ylim1, x.labels = F){
    x.labels <- ifelse(x.labels, "s", "n")
    boxplot(x, xaxt = x.labels, las = 2, col = "white", ylim = ylim1)

    mapply(function(i){
      x_i <- x[,i]
      if (outliers::grubbs.test(x_i)$p.value < 0.001) {
        outs = 100000 #- the initial drops
        p.value = 0
        p.val = c()
        while(p.value < 0.001){
          p.value <- outliers::grubbs.test(x_i[!names(x_i) %in% outs])$p.value
          if (p.value < 0.001){
            outcase <- as.numeric(names(p.value))
            outs <- c(outs, outcase)
            p.val <- c(p.val, p.value)
          }
        }
        outs <- outs[-1] #- remove 100000
        #- outs.01 <- outs[p.val >= 0.001]
        #- points(rep(i, length(outs.01)), x_i[match(outs.01, names(x_i))], 
        #         pch = 19, col = "black", cex = 0.9)
        
        outs.001 <- outs[p.val < 0.001]
        points(rep(i, length(outs.001)), x_i[match(outs.001, names(x_i))], 
               pch = 19, col = "red", cex = 0.9)
        
        out.txt <- outs[1:2]
        text(rep(i, length(out.txt)), x_i[match(out.txt, names(x_i))], 
             labels = out.txt, pos = 4, cex = 0.7)
      }
    }, i = 1: np)
  }
  
  par(mar = c(0, 4, 3.5, 1))
  if(is.null(ylim1)) ylim1 = c(min(x1, x2), max(x1, x2))
  pbox(x1, ylim1, x.labels = F)
  mtext("Group 1", side = 2, line = 2.3)
  mtext(name, side = 3, line = 0.2, adj = 0)
  
  par(mar = c(4, 4, 0, 1))
  if(is.null(ylim2)) ylim2 = c(min(x1, x2), max(x1, x2))
  pbox(x2, ylim1, x.labels= T)
  mtext("Group 2", side = 2, line = 2.3)
}

get.layoutmatrix <- function(nrow, ncol){
  m <- m0 <- matrix(1: (2 * ncol), ncol = ncol, byrow = F)
  if (nrow > 1)
    for(i in 2:nrow){
      add <- (i - 1) * ncol * 2 
      m <- rbind(m, m0 + add)
    }
  return(m)
}

jpeg("./output/box-trans-reg.jpeg", width = 7, height = 9, units = "in",
     res = 300)
layout(get.layoutmatrix(3, 2))
par(mar = c(4, 4, 2, 2))
for(s in d[4:6]){
  
  load(paste0(dir, s))
  rownames(X) <- paste0("V", 1:nrow(X))
  X <- t(X)
  y <- as.factor(y)
  Z <- transdata(X, byrow = F, eps = 0.01, scale = F)
  X <- scale(X)
  Z <- scale(Z)
  modX <- glmboost(X, y, 
                   family = Binomial(),
                   control = boost_control(mstop = 100, 
                                           nu = 0.1, 
                                           risk = "inbag"), 
                   center = T) 
  
  modZ <- glmboost(Z, y, 
                   family = Binomial(),
                   control = boost_control(mstop = 100, 
                                           nu = 0.1, 
                                           risk = "inbag"), 
                   center = T) 
  
  betax <- coef(modX)
  betaz <- coef(modZ)
  onlyx <- setdiff(names(betax), names(betaz))
  onlyz <- setdiff(names(betaz), names(betax))
  both  <- intersect(names(betax), names(betaz))
  pbox(onlyx, X, y, name = paste0(gsub(".RData", "", s), "(only in OD)"))
  pbox(onlyz, X, y, name = paste0(gsub(".RData", "", s), "(only in TD)"))
}
dev.off()

#- box plot for survival in regression

pboxsurv <- function(vars, X, name){
  
  colnames(X) <- `if`(is.null(colnames(X)), paste0("V", 1:ncol(X)), colnames(X))
  rownames(X) <- 1:nrow(X)
  np <- length(vars)
  
  x <- X[, vars, drop = F]
  colnames(x) <- vars
  
  pbox <- function(x, x.labels = F){
    boxplot(x,  las = 2, col = "white", show.names = T)
    mapply(function(i){
      x_i <- x[,i]
      if (outliers::grubbs.test(x_i)$p.value < 0.001) {
        outs = 100000 #- the initial drops
        p.value = 0
        p.val = c()
        
        while(p.value < 0.001){
          p.value <- outliers::grubbs.test(x_i[!names(x_i) %in% outs])$p.value
          if (p.value < 0.001){
            outcase <- as.numeric(names(p.value))
            outs <- c(outs, outcase)
            p.val <- c(p.val, p.value)
          }
        }
        outs <- outs[-1] #- remove 100000
#        outs.01 <- outs[p.val >= 0.001]
#        points(rep(i, length(outs.01)), x_i[match(outs.01, names(x_i))], 
#               pch = 19, col = "black", cex = 0.9)
        
        outs.001 <- outs[p.val < 0.001]
        points(rep(i, length(outs.001)), x_i[match(outs.001, names(x_i))], 
               pch = 19, col = "red", cex = 0.9)
        
        out.txt <- outs[1:2]
        text(rep(i, length(out.txt)), x_i[match(out.txt, names(x_i))], 
             labels = out.txt, pos = 4, cex = 0.7)
      }
    }, i = 1: np)
  }
  pbox(x, x.labels = F)
  mtext(name, side = 3, line = 0.2, adj = 0)
}

jpeg("./output/box-trans-reg-surv.jpeg", width = 7, height = 6, units = "in",
     res = 300)
par(mar = c(4, 4, 2, 2), mfrow = c(2, 2))
for(i in 7:8){
  
  load(paste0(dir, d[i]))
  rownames(X) <- paste0("V", 1:nrow(X))
  X <- t(X)
  family = type[[i]]
  
  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv 
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  Z <- transdata(X, byrow = F, eps = 0.01, scale = F)
  X <- scale(X)
  Z <- scale(Z)
  modX <- glmboost(X, y, 
                   family = family,
                   control = boost_control(mstop = 100, 
                                           nu = 0.1, 
                                           risk = "inbag"), 
                   center = T) 
  
  modZ <- glmboost(Z, y, 
                   family = family,
                   control = boost_control(mstop = 100, 
                                           nu = 0.1, 
                                           risk = "inbag"), 
                   center = T) 
  
  betax <- coef(modX)
  betaz <- coef(modZ)
  onlyx <- setdiff(names(betax), names(betaz))
  onlyz <- setdiff(names(betaz), names(betax))
  both  <- intersect(names(betax), names(betaz))
  
  pboxsurv(onlyx, X, name = paste0(gsub(".RData", "", d[i]), "(only in OD)"))
  pboxsurv(onlyz, X, name = paste0(gsub(".RData", "", d[i]), "(only in TD)"))
}
dev.off()

#-------------------------------------------------------------------------------
#- effects of transformation based on iteration path changes
#-------------------------------------------------------------------------------

trans.path.changes <- function(X, y, mstop = 100,
                               family){
  
  require(mboost)
  colnames(X) <- paste0("V", 1:ncol(X))
  X <- scale(X)
  Z <- transdata(X, byrow = F, eps = 0.01, scale = F)
  Z <- scale(Z)

  modX <- glmboost(X, y, 
                   family = family,
                   control = boost_control(mstop = mstop, 
                                           nu = 0.1, 
                                           risk = "inbag"), 
                   center = T) 
  tmp <- coef(modX, aggregate = "cumsum", off2int = F)
  pathX <- Reduce(rbind, tmp)
  rownames(pathX) <- names(tmp)
  
  modZ <- glmboost(Z, y, 
                   family = family,
                   control = boost_control(mstop = mstop, 
                                           nu = 0.1, 
                                           risk = "inbag"), 
                   center = T) 
  tmp <- coef(modZ, aggregate = "cumsum", off2int = F)
  pathZ <- Reduce(rbind, tmp)
  rownames(pathZ) <- names(tmp)
  
  allvar <- unique(c(rownames(pathX), rownames(pathZ)))
  matZ <- matX <- matrix(0, nrow = length(allvar), ncol = mstop, 
                         dimnames = list(allvar, 1:mstop))
  
  matX[rownames(pathX),] <- pathX
  matZ[rownames(pathZ),] <- pathZ
  
  dismat <- abs(matX - matZ)
  disout <- apply(dismat, 1, sum)
  names(disout) <- rownames(dismat)

  varX <- setdiff(rownames(pathX), rownames(pathZ))
  varZ <- setdiff(rownames(pathZ), rownames(pathX))
  varboth <- intersect(rownames(pathX), rownames(pathZ))
  out <- list(mat = dismat, score = disout, 
              varX = varX, varZ = varZ, varboth = varboth)
  class(out) <- append(class(out), "trans.path.changes")
  return(out)
}

plot.trans.path.changes <- function(obj, path = T,
                                    cex.axis = 0.5,
                                    cex.text = 0.5
                                    ){
  
  if(path == T) {
    score <- obj$score
  } else{
    score <- obj$mat[,ncol(obj$mat)]
    names(score) <- rownames(obj$mat)
  }
  
  score <- sort(score, decreasing = T)
  varX  <- obj$varX
  varboth <- obj$varboth
  varZ  <- obj$varZ
  
  #- plot 
  score <- c(score[names(score) %in% varX],
             score[names(score) %in% varboth],
             score[names(score) %in% varZ])
  axisx <- 1:length(varX)
  axisboth <- (length(varX) + 1.25) : (length(c(varX, varboth)) + 0.25)
  axisz <- (length(c(varX, varboth)) + 1.5) : (length(c(varX, varboth, varZ)) + 0.5)
  
  varname <- names(score)
  y = score
  x = c(axisx, axisboth, axisz)
  plot.new()
  plot.window(ylim = c(0, max(y)* 1.2), xlim = c(1, max(x)),
              yaxs="i")
  axis(1, las = 2, labels = varname, at = x, cex.axis = cex.axis, 
       tick = F, line = -0.5)
  axis(2)
  box()
  
  points(x, y, cex = 0.6, pch = 19)
  sapply(x, function(i){
    lines(c(i, i), c(0, y[i]))
  })
  
  abline(v = c((length(varX) + 0.5), (length(c(varX, varboth)) + 0.75)), 
         lty = "dashed")
  axis(3, at = c(mean(axisx), mean(axisboth), mean(axisz)),
       labels = c("in O", "in Both", "in T"), tick = F, line = -1)
  
  ord <- order(y, decreasing = T)[1:5]
  text(x[ord], y[ord], cex = cex.text, pos = 3, labels = varname[ord])
}  

#- path = TRUE
pdf("./output/trans.path.pdf", width = 6.5, height = 5)
par(mar = c(4, 3.5, 3, 2))
for(i in 1:6){
  load(paste0(dir, d[[i]]))
  X <- t(X)
  
  if (isTRUE(all.equal(type[[i]], CoxPH())))
    y <- surv 
  if (isTRUE(all.equal(type[[i]], Binomial())))
    y <- as.factor(y)
  
  obj <- trans.path.changes(X, y, family = type[[i]])
  plot.trans.path.changes(obj, path = T)
  sname <- paste0("bold(", gsub(".RData", "", i), ")")
  mtext(parse(text = sname), line = 1.2, side = 3, lwd = 2)
  mtext(bquote(Delta ~ beta), line = 2.2, side = 2, lwd = 2)
}
dev.off()

#- path = FALSE
pdf("./output/trans.beta.pdf", width = 6.5, height = 5)
par(mar = c(4, 3.5, 3, 2))
for(i in 1:6){
  load(paste0(dir, d[i]))
  X <- t(X)
  if (isTRUE(all.equal(type[[i]], CoxPH())))
    y <- surv 
  if (isTRUE(all.equal(type[[i]], Binomial())))
    y <- as.factor(y)
  
  obj <- trans.path.changes(X, y, family = type[[i]])
  plot.trans.path.changes(obj, path = F)
  sname <- paste0("bold(", gsub(".RData", "", i), ")")
  mtext(parse(text = sname), line = 1.2, side = 3, lwd = 2)
  mtext(bquote("|"~Delta~beta~"|"), line = 2.2, side = 2, lwd = 2)
}
dev.off()

#- box-trans2 for scherzer07 --
jpeg("./output/box-trans-2.jpeg", width = 6, height = 5, units = "in",
     res = 300)
layout(matrix(c(1, 1, 1, 1, 2, 3, 4, 5), nrow = 2, ncol = 4, byrow = T))
par(mar = c(4, 3.5, 2, 2))
i = "scherzer07.RData"
load(paste0(dir, i))
X <- t(X)
colnames(X) <- NULL
obj <- trans.path.changes(X, y)
plot.trans.path.changes(obj, path = F)
mtext(bquote("|"~Delta~beta~"|"), line = 2.2, side = 2, lwd = 2)

#- boxplot --
X <- scale(X) 
Z <- transdata(X, byrow = F, eps = 0.01, scale = F)
Z <- scale(Z)
y <- as.numeric(as.factor(y)) - 1
vars <- names(sort(obj$score, decreasing = T))[1:2]

par(mar = c(3, 2, 1, 2))
boxp(vars[1], vars[1])
boxp(vars[2], vars[2])

dev.off()

#-------------------------------------------------------------------------------
#- add plot
#-------------------------------------------------------------------------------

hatplot <- function(X, set, ylim = NULL, name){
  
  X <- X[, set, drop = F]
  X <- scale(X)
  n <- nrow(X)
  p <- ncol(X)
  Z <- transdata(X, byrow = F, eps = 0.01, scale = F)
  Z <- scale(Z)
  
  levx <- hat(X, intercept = F)
  levz <- hat(Z, intercept = F)
  n = length(levx)
  if(is.null(ylim)) ylim = c(0, max(levx) * 1.2)
  
  plot(1:n, levx, xlab = "", ylab = "", pch = 1, cex = 0.7, 
       ylim = ylim)
  s <- which(levx >= (5 * p / n))
  if(length(s) > 0)
    text(s, levx[s], labels = s, pos = 3, cex = 0.9)
  mtext("Leverage", side = 2, line = 2.2)
  mtext("Observations", side = 1, line = 2.2)
  mtext(name, side = 3, line = 0.2, adj = 0)
  
  plot(1:n, levz, xlab = "", ylab = "", pch = 1, cex = 0.7,
       ylim = ylim)
  s <- which(levz >= (5 * p / n))
  if(length(s) > 0)
    text(s, levz[s], labels = s, pos = 3, cex = 0.9)
  mtext("Leverage", side = 2, line = 2.2)
  mtext("Observations", side = 1, line = 2.2)
  mtext(name, side = 3, line = 0.2, adj = 0)
}

jpeg(paste0("./output/lev.jpeg"), width = 8, height = 10, 
     units = "in", res = 300)

mat <- rbind(1:4, matrix(5:12, ncol = 4, byrow = F),
            matrix(13:20, ncol = 4, byrow = F))
layout(mat, widths = c(0.5, 4, 4, 4), heights = c(0.5, 4, 4, 4, 4))
par(mar = c(0, 0, 0, 0))
plot.new()
plot.new()
text(0.5, 0.5, label = "Vars only in OD", cex = 1.5)
plot.new()
text(0.5, 0.5, label = "Vars only in TD", cex = 1.5)
plot.new()
text(0.5, 0.5, label = "Vars in Both", cex = 1.5)

#- 1, 4
#- 2, 3
#- 5, 6
#- 7, 8

for(i in c(7, 8)){
  
  load(paste0(dir, d[i]))
  X <- t(X)
  
  if (isTRUE(all.equal(type[[i]], CoxPH())))
    y <- surv 
  if (isTRUE(all.equal(type[[i]], Binomial())))
    y <- as.factor(y)
  
  sets  <- transregress(X, y, family = type[[i]])
  onlyX <- setdiff(sets$X, sets$Z)
  onlyZ <- setdiff(sets$Z, sets$X)
  both  <- intersect(sets$X, sets$Z)
  
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(0.5, 0.5, label = "Original Data", cex = 1.5, srt = 90)
  plot.new()
  text(0.5, 0.5, label = "Transformed Data", cex = 1.5, srt = 90)
  
  name = gsub(".RData", "", d[i])
  hatplot(X, onlyX, name = name)
  hatplot(X, onlyZ, name = name)
  hatplot(X, both, name = name)
}
dev.off()

#-------------------------------------------------------------------------------
#- Tietjen-Moore Test
#-------------------------------------------------------------------------------

tm <- function(x, k){
  
  n = length(x)
  
  ## Compute the absolute residuals.
  r = abs(x - mean(x))
  
  ## Sort data according to size of residual.
  df = data.frame(x, r)
  dfs = df[order(df$r),]
  
  ## Create a subset of the data without the largest k values.
  klarge = c((n-k+1):n)
  subx = dfs$x[-klarge]
  
  ## Compute the sums of squares.
  ksub = (subx - mean(subx))**2
  all = (df$x - mean(df$x))**2
  
  ## Compute the test statistic.
  ek = sum(ksub)/sum(all)
  return(ek)
}


#----

plotfun.addtop10 <- function(obj, range = 100, topn = 10, 
                             labels = F, labels.topn = 1){
  rank1 <- obj$rank1
  rank2 <- obj$rank2
  
  ord <- order(abs(rank1 - rank2), decreasing = T)[1:topn]
  ex <- which(rank1[ord] > range | rank2[ord] > range)
  exrank1 <- pmin(rank1[ord], range)[ex]
  exrank2 <- pmin(rank2[ord], range)[ex]
  points(exrank1, exrank2, pch = 19, col = "red")
  points(rank1[ord], rank2[ord], col = "red", pch = 19)
  if(labels){
    ord <- order(abs(rank1 - rank2), decreasing = T)[1:labels.topn]
    exrank1 <- pmin(rank1[ord], 100)
    exrank2 <- pmin(rank2[ord], 100)
    text(exrank1, exrank2, labels = obj$genes[ord], pos = 2, cex = 0.8)
  }
}

#- Boxplot of two groups

pbox <- function(obj, X, name, topn = NULL, ylim1 = NULL, ylim2 = NULL){
  
  X <- apply(X, 2, scale)
  rownames(X) <- 1:nrow(X)
  
  rank1 <- obj$rank1
  rank2 <- obj$rank2
  
  d_r <- abs(rank1 - rank2)
  g   <- order(d_r, decreasing = T)[1:topn]
  g   <- obj$genes[g]
  
  x1 <- X[y == levels(y)[1], g, drop = F]
  x2 <- X[y == levels(y)[2], g, drop = F]
  colnames(x1) <- colnames(x2) <- g
  
  par(mar = c(0, 4, 3.5, 1))
  if(is.null(ylim1)) ylim1 = c(min(x1, x2), max(x1, x2))
  boxplot(x1, xaxt = "n", las = 2, col = "white", ylim = ylim1)
  
  #- get critical
  fun.cut <- function(x, k, p){
    test = c(1:10000)
    for (i in 1:10000){
      xx = rnorm(length(x))
      test[i] = tm(xx, k)
      }
    quantile(test, p)
  }
  
  c1 = c2 = c() #- cutoffs
  for (i in 1:ncol(x1)) {
    
    x_i <- x1[,i]
    j = 1
    if(!is.numeric(c1[j])) 
      c1[j] <- fun.cut(x_i, j, p = 0.001)
    
    while (tm(x_i, j) < c1[j]) {
      j = j + 1
      if(is.na(c1[j]))
        c1[j] = fun.cut(x_i, j, p = 0.001)
    }
    
    j = j - 1
    r = abs(x_i - mean(x_i))
    if (j >= 1)
      ind = order(r, decreasing = T)[1:j] else
        ind = NULL
    
    out <- x_i[ind]
    points(rep(i, length(out)), out, pch = 19, col = "red", cex = 0.7)
    
    out.txt <- ind[1:2]
    text(rep(i, length(out.txt)), out[1:2], 
         labels = out.txt, pos = 4, cex = 0.7)
  }
  
  mtext("Group 1", side = 2, line = 2.3)
  mtext(name, line = 0.2, side = 3, adj = 0)
  
  par(mar = c(4, 4, 0, 1))
  if(is.null(ylim2)) ylim2 = c(min(x1), max(x2))
  boxplot(x2, las = 2, col = "white", ylim = ylim2)
  mtext("Group 2", side = 2, line = 2.3)
  
  for(i in 1:ncol(x2)){
    x_i <- x2[, i]
    j = 1
    if(!is.numeric(c2[j])) 
      c2[j] <- fun.cut(x_i, j, p = 0.001)
    
    while (tm(x_i, j) < c2[j]) {
      j = j + 1
      if(is.na(c2[j]))
        c2[j] = fun.cut(x_i, j, p = 0.001)
    }
    
    j = j - 1
    r = abs(x_i - mean(x_i))
    if (j >= 1)
      ind = order(r, decreasing = T)[1:j] else
        ind = NULL
    
    out <- x_i[ind]
    points(rep(i, length(out)), out, pch = 19, col = "red", cex = 0.7)
    
    out.txt <- ind[1:2]
    text(rep(i, length(out.txt)), out[1:2], 
         labels = out.txt, pos = 4, cex = 0.7)
  }
}

#- get the matrix for layout
get.layoutmatrix <- function(nrow, ncol){
  m <- m0 <- matrix(1: (3 * ncol), ncol = ncol, byrow = F)
  if (nrow > 1)
    for(i in 2:nrow){
      add <- (i - 1) * ncol * 3 
      m <- rbind(m, m0 + add)
    }
  return(m)
}

get.height <- function(nrow) rep(c(2.3, 1.5, 1.5), nrow)

#- plot all data of binary outcomes

load("./output/re.transrank.RData")

#- main text figure
s = c("ressom06.RData", "scherzer07.RData", "petricoin02.RData")
#- appendix figure
s = c("singh02.RData", "wang05.RData", "golub99.RData")

jpeg("./output/rank-box-all.jpeg", width = 10, height = 6, 
     units = "in", res = 300)
layout(get.layoutmatrix(1, 3), heights = get.height(1))

for(i in 1:3){
  par(mar = c(3, 4, 2.5, 1))
  rk <- re.transrank[[s[i]]]
  plotfun(rk$rank1, rk$rank2, range = 100, labels = rk$genes, n.labels = 1)
  plotfun.addtop10(rk, range = 100, topn = 10, 
                   labels = F, labels.topn = 1)
  mtext(paste0(letters[i], "), ", gsub(".RData", "", s[i])), 
        side = 3, line = 0.2, adj = 0)
  load(paste0(dir, s[i]))
  X <- t(X)
  pbox(rk, X = X, name = gsub(".RData", "", s[i]), topn = 10)
}

dev.off()


#- explan the problem with this issue
s = c("ressom06.RData", "scherzer07.RData", "petricoin02.RData")

rk <- re.transrank[[s[2]]]
load(paste0(dir, s[2]))
X <- t(X)

X <- apply(X, 2, scale)
rownames(X) <- 1:nrow(X)

rank1 <- rk$rank1
rank2 <- rk$rank2

d_r <- abs(rank1 - rank2)
g   <- order(d_r, decreasing = T)[1:10]
g   <- rk$genes[g]

x1 <- X[y == levels(y)[1], g, drop = F]
x2 <- X[y == levels(y)[2], g, drop = F]
colnames(x1) <- colnames(x2) <- g

jpeg("./output/tmp.jpeg", width = 12, height = 9, units = "in", res = 300)
par(mfrow = c(3, 4), mar = c(2, 2, 2, 2))
for(i in 1:ncol(x1)){
  
  x = x1[,i]
  hist(x, main = colnames(x1)[i])
  #plot(density(x), xlim = c(-4, 4), main = colnames(x1)[i])
  
}

dev.off()

#-------------------------------------------------------------------------------

