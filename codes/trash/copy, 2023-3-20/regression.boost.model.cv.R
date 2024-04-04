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
#- IPs based on models
#-------------------------------------------------------------------------------

load("./output/cv.mstops.RData")

drop1boost <- function(s, family, method = c("fixed", "cv")){
  
  load(paste0(dir, s))
  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv 
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  cv = cv.mstops[[s]][["mstops"]]
  mstop = cv[length(cv)]
  obs.mstop = cv[-length(cv)]
  
  n <- ncol(X)
  p <- nrow(X)
  X <- t(X)
  
  if (method == "fixed") mstops = rep(mstop, (n + 1))
  if (method == "cv") mstops = c(obs.mstop, mstop)  
  
  cores <- detectCores() - 4
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  #- (n + 1) is the boosting without removing obs
  reboost <- foreach(i = 1:(n + 1), j = mstops, 
                     .packages = c("mboost", "survival")) %dopar%{
      
                       boost <- glmboost(X[-i,], y[-i], 
                                         family = family,
                                         control = boost_control(mstop = j, 
                                                                 nu = 0.1, 
                                                                 risk = "inbag"), 
                                         center = T) 
                       return(names(coef(boost)))
                     }
  registerDoSEQ()
  stopCluster(cl)
  return(reboost)
}

#- transformed data -
boost.mod.fix <- mapply(drop1boost, s = d, family = type, 
                        MoreArgs = list(method = "fixed"), SIMPLIFY = F)
#- save(list = "boost.mod.fix", file = "./output/boost.mod.fix.RData")

boost.mod.cv <- mapply(drop1boost, s = d, family = type, 
                       MoreArgs = list(method = "cv"), SIMPLIFY = F)
#- save(list = "boost.mod.cv", file = "./output/boost.mod.cv.RData")

#- make barplot, leave-one-out--------------------------------------------------
#- ("M1", "M2") --

method = "M2"
load("./output/boost.mod.fix.RData")
load("./output/boost.mod.cv.RData")

#- bar, barplot
#- lol, lollipop plot

plotm <- function(obj, method, name, xlim, plot = c("bar", "lol"),
                  center = F, refline = NULL){
    #- pog, compare to the original model
    n = length(obj) - 1
    boost_all   <- obj[[n + 1]]   
    boost_drop1 <- obj[-(n + 1)]
    
    if (method == "M1") {
      re.dif <- mapply(ham, sj = boost_drop1, MoreArgs = list(si = boost_all))
      names(re.dif) <- paste0("obs", 1:n)
    }
    
    if (method == "M2") {
      re.dif <- moddist(mod = boost_drop1, std = F, L = 1)
      names(re.dif) <- paste0("obs", 1:n)
    }
    
    if (center)
      re.dif.sd <- re.dif / sd(re.dif)
    else
      re.dif.sd <- re.dif / sd(re.dif)
#   
    if (plot == "bar")
      hor.barplot(re.dif.sd, digit = 2, decreasing = T, topn = 20, xlim = xlim,
                  score = T)
    
    if (plot == "lol")
      lollipop(re.dif.sd, decreasing = T, topn = 5, ylim = xlim, 
               refline = refline, ylab = method)
    
    mtext(name, side = 3, line = 0.2, adj = 0)
}

jpeg(paste0("./output/boost.mod.m1.main.jpeg"), width = 9, height = 6, 
     units = "in", res = 300)

layout(matrix(1:6, nrow = 2, byrow = F))
par(mar = c(3.5, 3.5, 2, 1))

for(s in d[c(4, 5, 8)]){
  
  plotm(boost.mod.fix[[s]], method = "M1", 
        name = paste0(gsub(".RData", "", s), "(fixed)"),
        xlim = c(-3, 7), plot = "lol", center = T, 
        refline = c(-2, 2))
  
  plotm(boost.mod.cv[[s]],  method = "M1", 
        name = paste0(gsub(".RData", "", s), "(adaptive)"),
        xlim = c(-3, 7), plot = "lol", center = T,
        refline = c(-2, 2))
}

dev.off()

jpeg(paste0("./output/boost.mod.m1.append.jpeg"), width = 9, height = 6, 
     units = "in", res = 300)

layout(matrix(1:6, nrow = 2, byrow = F))
par(mar = c(3.5, 3.5, 2, 1))
for(s in d[c(2, 3, 6)]){
  
  plotm(boost.mod.fix[[s]], method = "M1", 
        name = paste0(gsub(".RData", "", s), "(fixed)"),
        xlim = c(-3, 7), plot = "lol", center = T,
        refline = c(-2, 2))
  
  plotm(boost.mod.cv[[s]],  method = "M1", 
        name = paste0(gsub(".RData", "", s), "(adaptive)"),
        xlim = c(-3, 7), plot = "lol", center = T,
        refline = c(-2, 2))
}

dev.off()

#- m2
jpeg(paste0("./output/boost.mod.m2.main.jpeg"), width = 9, height = 6, 
     units = "in", res = 300)
layout(matrix(1:6, nrow = 2, byrow = F))
par(mar = c(3.5, 3.5, 2, 1))
for(s in d[c(4, 5, 8)]){
  plotm(boost.mod.fix[[s]], method = "M2", 
        name = paste0(gsub(".RData", "", s), "(fixed)"),
        xlim = c(0, 12), plot = "lol", center = F,
        refline = NULL)
  plotm(boost.mod.cv[[s]],  method = "M2", 
        name = paste0(gsub(".RData", "", s), "(adaptive)"),
        xlim = c(0, 12), plot = "lol", center = F,
        refline = NULL)
}
dev.off()

jpeg(paste0("./output/boost.mod.m2.append.jpeg"), width = 9, height = 6, 
     units = "in", res = 300)

layout(matrix(1:6, nrow = 2, byrow = F))
par(mar = c(3.5, 3.5, 2, 1))
for(s in d[c(2, 3, 6)]){
  
  plotm(boost.mod.fix[[s]], method = "M2", 
        name = paste0(gsub(".RData", "", s), "(fixed)"),
        xlim = c(0, 12), plot = "lol", center = F,
        refline = NULL)
  
  plotm(boost.mod.cv[[s]],  method = "M2", 
        name = paste0(gsub(".RData", "", s), "(adaptive)"),
        xlim = c(0, 12), plot = "lol", center = F,
        refline = NULL)
}
dev.off()

#- compare M1 and M2 -----------------------------------------------------------

plot.m1vsm2 <- function(obj, plot.type = c("blandartman", "corrplot"), name){
  
  plot.type <- match.arg(plot.type, choices = c("blandartman", "corrplot"))
  n = length(obj) - 1
  boost_all   <- obj[[n + 1]]   
  boost_drop1 <- obj[-(n + 1)]
  
  re.dif1 <- mapply(ham, sj = boost_drop1, MoreArgs = list(si = boost_all))
  re.dif1 <- re.dif1 / sd(re.dif1)
  
  names(re.dif1) <- paste0("obs", 1:n)
  
  re.dif2 <- moddist(mod = boost_drop1, std = F, L = 1)
  re.dif2 <- re.dif2 / sd(re.dif2)
  names(re.dif2) <- paste0("obs", 1:n)
  
  if (plot.type == "blandartman"){
    blandartman(re.dif1, re.dif2, spline = F, ci = T, ci.limit = F)
    mtext("(M1 + M2) / 2", side = 1, line = 2.2)
    mtext("M1 - M2", side = 2, line = 2.2)
  }
  
  if (plot.type == "corrplot"){
    corrplot(re.dif1, re.dif2, ylim = getlimits(re.dif2),
             xlim = getlimits(re.dif1),
             ci.line = F)
    mtext("M1", side = 1, line = 2.2)
    mtext("M2", side = 2, line = 2.2)
    
    ord1 <- order(re.dif1, decreasing = T)[1:3]
    ord2 <- order(re.dif2, decreasing = T)[1:3]
    ord  <- unique(ord1, ord2)
    
    text(re.dif1[ord], re.dif2[ord], labels = ord, pos = 2, cex = 0.9)
  }
  mtext(name, side =3, line = 0.2, adj = 0)
}

load("./output/boost.mod.fix.RData")
load("./output/boost.mod.cv.RData")

jpeg("./output/m1vsm2.jpeg", width = 9, height = 3, units = "in", res = 300)
par (mfrow = c(1, 3), mar = c(4, 3.5, 1.5, 1.5))
for(s in d[c(4, 5, 8)]){
  obj <- boost.mod.cv[[s]]
  plot.m1vsm2(obj, plot.type = "c", name = gsub(".RData", "", s))
}
dev.off()

#-------------------------------------------------------------------------------
#- M1 and M2, using adaptive tuning
#-------------------------------------------------------------------------------

# load("./output/boost.mod.fix.RData")
load("./output/boost.mod.cv.RData")
jpeg("./output/m1m2.jpeg", width = 8, height = 5, units = "in", res = 300)

mat <- matrix(1:6, nrow = 2, byrow = F)
layout(mat)
par(mar = c(3.5, 3.5, 1.5, 1.5), oma = c(1, 1, 1, 1))
for(s in d[c(4, 5, 8)]){
  
  obj <- boost.mod.cv[[s]]
  
   # plotm(obj, method = "M1", 
   #       name = paste0(gsub(".RData", "", s), "(adaptive)"),
   #       xlim = c(0, 12), plot = "lol", center = T)
   plotm(obj,  method = "M2",
         name = paste0(gsub(".RData", "", s), "(adaptive)"),
         xlim = c(0, 12), plot = "lol", center = F)
   
   plot.m1vsm2(obj, plot.type = "c", name = "")
}

dev.off()


load("./output/boost.mod.cv.RData")
load("./output/boost.mod.fix.RData")

jpeg("./output/m1m2.append.jpeg", width = 9, height = 5, 
     units = "in", res = 300)

mat <- matrix(1:6, nrow = 2, byrow = F)

layout(mat)
par(mar = c(3.5, 3.5, 1.5, 1.5), oma = c(1, 1, 1, 1))
for(s in d[c(2, 3, 6)]){
  
  obj <- boost.mod.cv[[s]]
  plotm(obj,  method = "M2", 
        name = paste0(gsub(".RData", "", s), "(adaptive)"),
        xlim = c(0, 12), plot = "lol", center = F)
  
  plot.m1vsm2(obj, plot.type = "c", name = "")
}

dev.off()

# - compare fixed and adaptive tuning in M2 ------------------------------------

plot.fix.vs.ada <- function(obj.fix, obj.ada,
                            method = c("M1", "M2"),
                            plot.type = c("blandartman", "corrplot"), name){
  
  plot.type <- match.arg(plot.type, choices = c("blandartman", "corrplot"))

  n = length(obj.fix) - 1
  mods.fix <- obj.fix[-(n + 1)]
  mods.ada <- obj.ada[-(n + 1)]
  
  if (method == "M2"){
    re.fix <- moddist(mod = mods.fix, std = F, L = 1)
    re.fix <- re.fix / sd(re.fix)
    re.ada <- moddist(mod = mods.ada, std = F, L = 1)
    re.ada <- re.ada / sd(re.ada)
    names(re.fix) <- names(re.ada) <- paste0("obs", 1:n)
  }
  
  if (method == "M1"){
    
    boost_all   <- obj.fix[[n + 1]]   
    boost_drop1 <- obj.fix[-(n + 1)]
    
    re.fix <- mapply(ham, sj = boost_drop1, MoreArgs = list(si = boost_all))
    re.fix <- scale(re.fix)
    
    boost_all   <- obj.ada[[n + 1]]   
    boost_drop1 <- obj.ada[-(n + 1)]
    
    re.ada <- mapply(ham, sj = boost_drop1, MoreArgs = list(si = boost_all))
    re.ada <- scale(re.ada)
    names(re.fix) <- names(re.ada) <- paste0("obs", 1:n)
  }
  
  # 
  # if (plot.type == "blandartman"){
  #   blandartman(re.fix, re.ada, spline = F, ci = T, ci.limit = F)
  #   mtext(xlab, side = 1, line = 2.2)
  #   mtext(ylab, side = 2, line = 2.2)
  # }

  if (plot.type == "corrplot"){
    corrplot(re.fix, re.ada, ylim = getlimits(re.ada),
             xlim = getlimits(re.fix),
             ci.line = F)
    mtext("adaptive", side = 2, line = 2.2)
    mtext("fixed", side = 1, line = 2.2)

    ord1 <- order(re.fix, decreasing = T)[1:3]
    ord2 <- order(re.ada, decreasing = T)[1:3]
    ord  <- unique(c(ord1, ord2))

    text(re.fix[ord], re.ada[ord], labels = ord, pos = 2, cex = 1)
  }
  mtext(name, side =3, line = 0.2, adj = 0)
}

#-------------------------------------------------------------------------------

load("./output/boost.mod.cv.RData")
load("./output/boost.mod.fix.RData")

jpeg("./output/fixvsadaM1.jpeg", width = 9, height = 6, 
     units = "in", res = 300)

par(mfrow = c(2, 3), oma = c(1, 1, 1, 1), mar = c(3.5, 3.5, 1.5, 1.5))

for(s in d[c(2:6, 8)]){
  
  obj.fix <- boost.mod.fix[[s]]
  obj.ada <- boost.mod.cv[[s]]
  
  plot.fix.vs.ada(obj.fix, obj.ada, plot.type = "c", method = "M1",
                  name = gsub(".RData", "", s))
}

dev.off()

#-------------------------------------------------------------------------------
#- outlierness of obs100 in petricoin02
#-------------------------------------------------------------------------------

getpred <- function(X, y, drop, mstop0, mstop1){
  boost <- glmboost(X, y, 
                    family = Binomial(),
                    control = boost_control(mstop = mstop0, 
                                            nu = 0.1, 
                                            risk = "inbag"), 
                    center = T)
  p <- predict(boost, newdata = X)
  
  #- remove m
  boost.m <- glmboost(X[-drop,], y[-drop], 
                      family = Binomial(),
                      control = boost_control(mstop = mstop1, 
                                              nu = 0.1, 
                                              risk = "inbag"), 
                      center = T) 
  p.m <- predict(boost.m, newdata = X)
  return(list(p = p, p.m = p.m))
}

#- predict p values in box plot
plot.p.box <- function(X, y, cv.mstops){
  boost <- glmboost(X, y, 
                    family = Binomial(),
                    control = boost_control(mstop = cv.mstops, 
                                            nu = 0.1, 
                                            risk = "inbag"), 
                    center = T)
  p <- predict(boost, newdata = X)
  
  pp <- exp(p) / (1 + exp(p))
  yy = as.numeric(y) - 1
  
  pgroup <- split(pp, yy)
  boxplot(pgroup, ylim = c(0, 1))
  res <- order(pgroup$`0`, decreasing = T)[1]
  text(rep(1, 1), pgroup$`0`[res], labels = res, cex = 0.9, 
       pos = c(4, 4, 4, 2))
  mtext("probability", side = 2, line = 2.2)
  mtext("Groups", side = 1, line = 2.2)
}

load("./output/cv.mstops.RData")
s = d[5]
load(paste0(dir, s))
X <- t(X)
y <- as.factor(y)
cv.mstops <- cv.mstops[[s]][["mstops"]]
cv.mstops <- cv.mstops[length(cv.mstops)]
plot.p.box(X, y, cv.mstops = cv.mstops)


#- aggrement between original model and L-1-O model 
plot.p.cor <- function(X, y, drop, cv.mstop){

  re.p <- getpred(X, y, drop = drop, mstop0, mstop1)
  top <- function(p) exp(p) / (1 + exp(p))
  px <- top(re.p$p.m)
  py <- top(re.p$p)
  corrplot(px, py, p.col = "white")
  pxx <- split(px, y)
  pyy <- split(py, y)
  
  points(pxx$normal, pyy$normal, pch = 19, col = "steelblue4") # 
  points(pxx$cancer, pyy$cancer, pch = 19,  col = "tomato2") #
  text(px[m], py[m], pos = 4, cex = 0.8, labels = m)
  mtext(paste0("Model(-obs", m, ")"), side = 1, line =2.2)
  mtext("Original model", side = 2, line = 2.2)
  mtext("B", side = 3, line = 0.2, adj = 0)
  
  #- plot the median influential obseravtion
  load("./output/reboost.RData")
  reboost <- reboost[["petricoin02.RData"]]
  infscore <- moddist(reboost)[-length(reboost)]
  m <- order(infscore, decreasing = T)[floor(length(infscore)/2)] 
  re.p <- getpred(X, y, drop = m)
  top <- function(p) exp(p) / (1 + exp(p))
  px <- top(re.p$p.m)
  py <- top(re.p$p)
  corrplot(px, py, p.col = "white")
  
  pxx <- split(px, y)
  pyy <- split(py, y)
  points(pxx$normal, pyy$normal, pch = 19, col = "steelblue4")  
  points(pxx$cancer, pyy$cancer, pch = 19,  col = "tomato2")  
  text(px[m], py[m], pos = 4, cex = 0.8, labels = m)
  mtext(paste0("Model(-obs", m, ")"), side = 1, line =2.2)
  mtext("Original model", side = 2, line = 2.2)
  mtext("C", side = 3, line = 0.2, adj = 0)
}

#-------------------------------------------------------------------------------
#- hierarchical clustering
#-------------------------------------------------------------------------------

#- results in main text
load("./output/boost.mod.cv.RData")
s <- c("petricoin02.RData", "golub99.RData", "veer02.RData")
jpeg("./output/cluster.jpeg", width = 12, height = 4, units = "in",
     res = 300)
par(mfrow = c(1, length(s)), mar = c(0, 4, 2, 1))
for (i in s){
  boost  <- boost.mod.cv[[i]]
  boost  <- boost[-length(boost)]
  plotcluster(boost, dist = "binary")
  mtext(gsub(".RData", "", i), side = 3, line = 0)
}
dev.off()

#- results in appendix
jpeg("./output/cluster-append.jpeg", width = 12, height = 4, units = "in",
     res = 300)
par(mar = c(1, 4, 2, 3), mfrow = c(1, 3))
for(i in d[c(2, 3, 6)]){
  boost  <- boost.mod.cv[[i]]
  boost  <- boost[-length(boost)]
  plotcluster(boost)
  mtext(gsub(".RData", "", i), side = 3, line = 0)
}
dev.off()

#- hierachical clustering based on Kulczynski---

dist.kul <- function(X){
  
  kul <- function(v1, v2){
    
    s1 <- sum(v1 %in% v2)
    s2 <- length(v1)
    s3 <- length(v2)
    1 - (s1 / (2 * s2) + s1 / (2 * s3))
  }
  
  n <- ncol(X)
  mat <- matrix(NA, n, n)
  for(i in 1:n){
    for(j in i:n){
      v1 <- which(X[,i] != 0)
      v2 <- which(X[,j] != 0)
      mat[i, j] <- mat[j, i] <- kul(v1, v2)
    }
  }
  
  mat <- mat[lower.tri(mat)]
  
  attr(mat, "Size") <- n
  attr(mat, "Diag") <- F
  attr(mat, "Upper") <- F
  attr(mat, "method") <- "kul"
  return(mat)
}

plotcluster.kul <- function(boost, nlabs = 5){
  
  allvars <- unique(unlist(boost))
  boost <- mapply(function(x){
    as.numeric(allvars %in% x)
  }, x =  boost)
  
  d <- dist.kul(boost)
  d[is.na(d)] <- 0
  
  fit <- hclust(d, method = "average")
  n = length(fit$height) + 1
  points <- t(fit$merge)
  points <- rev(-points[points < 0]) 
  lab <- rep("", n)
  lab[points[1:nlabs]] <- points[1:nlabs]
  plot(fit, cex = 0.9,  main = "", labels = lab, frame.plot = F)
  
}

load("./output/boost.mod.cv.RData")
s <- c("petricoin02.RData", "golub99.RData", "veer02.RData")
jpeg("./output/cluster.jpeg", width = 12, height = 4, units = "in",
     res = 300)
par(mfrow = c(1, length(s)), mar = c(0, 4, 2, 1))
for (i in s){
  boost  <- boost.mod.cv[[i]]
  boost  <- boost[-length(boost)]
  plotcluster.kul(boost)
  mtext(gsub(".RData", "", i), side = 3, line = 0)
}
dev.off()


jpeg("./output/cluster-append.jpeg", width = 8, height = 8, units = "in",
     res = 300)
par(mar = c(1, 4, 2, 3), mfrow = c(3, 2))
for(i in d[-c(4, 5, 8)]){
  print(i)
  boost  <- boost.mod.cv[[i]]
  boost  <- boost[-length(boost)]
  plotcluster.kul(boost)
  mtext(gsub(".RData", "", i), side = 3, line = 0)
}
dev.off()

#-------------------------------------------------------------------------------
#-hierarchical clustering in appendix, example
#-------------------------------------------------------------------------------

load("./output/boost.mod.cv.RData")
s <- "golub99.RData"

#- using the 1-10 models as an example
mod <- boost.mod.cv[[s]][1:10]

allvars <- unique(unlist(mod))
boost <- mapply(function(x){
  as.numeric(allvars %in% x)
}, x =  mod)

d <- dist.kul(boost)
d[is.na(d)] <- 0


jpeg("./output/hier.jpeg", width = 5, height = 4, units = "in", res = 300)
par(mar = c(1, 4, 1, 1))

nlabs = 10
fit <- hclust(d, method = "average")
n = length(fit$height) + 1
points <- t(fit$merge)
points <- rev(-points[points < 0]) 
lab <- rep("", n)
lab[points[1:nlabs]] <- points[1:nlabs]
plot(fit, cex = 0.9,  main = "",  frame.plot = F, xlab = "")

dev.off()

dat <- matrix(NA, 10, 10)
dat[lower.tri(dat)] <- d
write.csv(dat, "./output/d.csv")

#-------------------------------------------------------------------------------
#- tuning parameters on models 
#-------------------------------------------------------------------------------

load("./output/boost.mod.cv.RData")
load("./output/cv.mstops.RData")
load("./output/boost.mod.fix.RData")

#- get the top3 observations
plot.tuning <- function(i){
  obj <- boost.mod.cv[[d[i]]]
  boost_drop1 <- obj[-length(obj)]
  boost_all   <- obj[length(obj)]
  re.dif <- moddist(mod = boost_drop1, std = F, L = 1)
  m2.cv <- order(re.dif, decreasing = T)[1:3]
  
  obj <- boost.mod.fix[[d[i]]]
  boost_drop1 <- obj[-length(obj)]
  boost_all   <- obj[length(obj)]
  re.dif <- moddist(mod = boost_drop1, std = F, L = 1)
  m2.fixed <- order(re.dif, decreasing = T)[1:3]
  obs <- unique(c(m2.cv, m2.fixed))
  
  #- plot the mboost steps
  load(paste0(dir, d[i]))
  family = type[[i]]
  n <- ncol(X)
  p <- nrow(X)
  X <- t(X)
  X <- scale(X)
  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv 
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  cv = cv.mstops[[d[i]]][["mstops"]]
  mstop = cv[length(cv)]
  obs.mstop = cv[-length(cv)]
  max.steps = max(obs.mstop)
  
  boost <- glmboost(X, y, 
                    family = family,
                    control = boost_control(mstop = max.steps, 
                                            nu = 0.1, 
                                            risk = "inbag"), 
                    center = T)
  plot.boostpath(boost)
  abline(v = obs.mstop[obs], lty = "dashed")
  mtext(text = obs, at = obs.mstop[obs], line = 0.2, side = 3, cex = 0.8)
  abline(v = mstop, lty = "dashed")
  mtext(text = "orig", at = mstop, line = 0.2, side = 3, cex = 0.8)
  mtext(gsub(".RData", "", d[i]), adj = 0, line = 0.2, side = 3)
}

jpeg("./output/mod.tuning.main.jpeg", units = "in", 
     width = 13, height = 4, res = 300)
par(mfrow = c(1, 3), mar = c(3.5, 3.5, 2, 2))
plot.tuning(4)
plot.tuning(5)
plot.tuning(8)
dev.off()

#-------------------------------------------------------------------------------
#- tuning parameters bar plot, for slides in DAGSTAT
#-------------------------------------------------------------------------------

load("./output/boost.mod.cv.RData")
load("./output/cv.mstops.RData")
load("./output/boost.mod.fix.RData")

bar.tuning <- function(obj){
  
  boost_drop1 <- obj[-length(obj)]
  boost_all   <- obj[length(obj)]
  m2 <- moddist(mod = boost_drop1, std = F, L = 1)
  # m2 <- m2 / sd(m2)
  names(m2) <- paste0("obs", seq_along(m2))
  
  hor.barplot(m2, decreasing = T, topn = 10, score = F)
}

par(mfcol = c(2, 2), mar = c(1, 4, 1, 1))
for (i in c(4, 5)){
  
  obj <- boost.mod.fix[[i]]
  bar.tuning(obj)
  
  obj <- boost.mod.cv[[i]]
  bar.tuning(obj)
  
}

#-------------------------------------------------------------------------------
#- The effects of fixed and adaptive lambda on M2
#-------------------------------------------------------------------------------

boostmod.path <- function(X, y, family, obs.mstop = NULL, 
                          ncores = NULL){
  
  if (class(family) == "boost_family_glm")
    y = as.factor(y)
  if (class(family) == "boost_family")
    y = surv
  
  require(mboost)
  require(doParallel)
  
  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X)

  if (is.null(ncores)) 
    ncores <- detectCores() - 1
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  mstops = max(obs.mstop)
  reboost <- foreach(i = 1:n, .packages = c("mboost")) %dopar% {
    
                       boost <- glmboost(X[-i,], y[-i], 
                                         family = family,
                                         control = boost_control(mstop = mstops, 
                                                                 nu = 0.1, 
                                                                 risk = "inbag"), 
                                         center = F)
                       
                       mod <- boost$xselect()
                       return(mod)
                     }
  registerDoSEQ()
  stopCluster(cl)
  
  scores <- c()
  for(i in 1:mstops){
    mod <- lapply(reboost, function(x){
     unique(as.character(x[1:i])) 
    })
    scores <- rbind(scores, moddist(mod))
  }

  return(out = list(scores = scores, mstops =  obs.mstop))
}

plot.boostmod.path <- function(obj){
  
  mat <- obj$scores
  maxstop <- nrow(mat)
  plot.new()
  plot.window(xlim = c(1, maxstop), ylim = c(min(mat), max(mat)))
  axis(1)
  axis(2)
  box()
  
  x <- 1:maxstop
  n <- ncol(mat)
  apply(mat, 2, function(y) lines(x, y, col = "grey"))
  
  mstops <- obj$mstops
  score.cv <- mapply(function(case, m){
    return(mat[m, case])
  }, case = 1:ncol(mat), m = mstops[-length(mstops)])
  
  org <- mstops[length(mstops)]
  score.fix <- mat[org,]
  
  # add top2 lines
  # top2 <- unique(order(score.cv, decreasing = T)[1:2], 
  #                order(score.fix, decreasing = T)[1:2])

  top2 <- order(score.fix, decreasing = T)[1:2]
  
  for(i in top2){
    lines(x, mat[,i])
  }
  
  # add top2 points
  points(mstops[top2], score.cv[top2], col = "red", pch = 19)
  #points(rep(org, 2), score.fix[top2], pch = 19)

  #- add orig line
  abline(v = mstops[length(mstops)], lty = "dashed")
  mtext("orig.", line = 0.2, side = 3, at = mstops[length(mstops)])

  #- add the other points
  y = setdiff(1:n, top2)
  points(mstops[y], score.cv[y], cex = 0.7, pch = 19)
  #- add line nunbers
  at <- mat[nrow(mat), top2]
  mtext(top2, side = 4, at = at, line = 0, las = 2)
  axis(side = 3, at = mstops[top2], tick = T, cex = 0.6,  tcl = 0.3, labels = F)
  mtext(top2, side = 3, at = mstops[top2], line = 0)
  
  #- add labels
  mtext("Boosting iterations", line = 2.2, side = 1)
  mtext("M2", line = 2.2, side = 2)
}

load("./output/cv.mstops.RData")

re.boostmod.path <- list()
for(i in c(4, 5, 8)){
  
  print(d[i])
  load(paste0("../data/", d[i]))
  X <- t(X)
  
  family = type[[i]]
  cv = cv.mstops[[d[i]]][["mstops"]]
  
  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  re.boostmod.path[[d[i]]] <- boostmod.path(X, y, family, cv)
}

#- save(list = "re.boostmod.path", file = "./output/re.boostmod.path.RData")
load("./output/re.boostmod.path.RData")
jpeg("./output/tuning.mod.jpeg", width = 12, height = 3.5, units = "in",
     res = 300)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 2))
for(i in names(re.boostmod.path)){
  
  plot.boostmod.path(re.boostmod.path[[i]])
  mtext(gsub(".RData", "", i), side = 3, line = 0.2, adj = 0)
}

dev.off()

#-------------------------------------------------------------------------------
#- table, IF observations
#-------------------------------------------------------------------------------

load("./output/boost.mod.cv.RData")
load("./output/cv.mstops.RData")

out <- c()
for(i in c(4, 5, 8)){
  
  obj <- boost.mod.cv[[i]]
  #- get IF
  n = length(obj)
  
  m1 <- mapply(ham, sj = obj[-n], MoreArgs = list(si = obj[[n]]))
  m2 <- moddist(mod = obj, std = F, L = 1)

  obs.if <- order(m2, decreasing = T)[1:3]
  obs.m  <- which.min(abs(m2 - median(m2)))
  mod.orig <- obj[[n]]
  nam <- gsub(".RData", "", names(boost.mod.cv)[i])
  
  re <- c()
  for(j in c(obs.if, obs.m)){
    mod <- obj[[j]]
    tmp <- c(nam, j, 
            sum(!mod %in% mod.orig),
            sum(!mod.orig %in% mod),
            sum(mod %in% mod.orig),
            length(union(mod.orig, mod)),
            m1[j], round(m2[j], 1))
    re <- rbind(re, tmp)
  }

  colnames(re) <- c("name", "Obs", "inL", "inO", "inboth", "union", "m1", "m2")
  out <- rbind(out, re)
}      

write.csv(out, "./output/mod.ip.out.csv")


#-------------------------------------------------------------------------------
#- BA plot, between M1 and M2
#- level, tmp
#-------------------------------------------------------------------------------

method = "M2"
load("./output/boost.mod.fix.RData")
load("./output/boost.mod.cv.RData")

plotm <- function(obj, plot = c("cor", "BA")){
  
  n = length(obj) - 1
  boost_all   <- obj[[n + 1]]   
  boost_drop1 <- obj[-(n + 1)]
  
  m1 <- mapply(ham, sj = boost_drop1, MoreArgs = list(si = boost_all))
  m1 <- scale(m1)
  
  m2 <- moddist(mod = boost_drop1, std = F, L = 1)
  m2 <- m2 / sd(m2)
  names(m2) <- names(m1) <- paste0("obs", 1:n)
  
  if (plot == "cor")
    corrplot(abs(m1), m2, ci.line = F, xlab = "M1", ylab = "M2",
             xlim = c(min(c(abs(m1), m2)), max(c(abs(m1), m2)) * 1.1))
  
  if (plot == "BA")
    blandartman(abs(m1), m2, ci = F, spline = F, ci.limit = F, 
                xlab = "(M1 + M2) / 2", ylab = "M1 - M2")
  
}

jpeg("./output/BA.in.mod.jpeg", width = 9, height = 6, res = 300, units = "in")
par(mfcol = c(2, 3), mar = c(3.5, 3.5, 1, 1), oma = c(1, 1, 1, 1))

for(i in c(4, 5, 8)){
  
  nam <- d[i]
  plotm(boost.mod.cv[[i]], plot = "cor")
  mtext(gsub(".RData", "", nam), side = 3, line = 0.2)
  plotm(boost.mod.cv[[i]], plot = "BA")
}

dev.off()

#-------------------------------------------------------------------------------
#- BA plot. Ada and fix tuning
#-------------------------------------------------------------------------------

load("./output/boost.mod.fix.RData")
load("./output/boost.mod.cv.RData")

plotm <- function(nam, plot = c("cor", "BA"), method = c("M1", "M2")){
  
  fix <- boost.mod.fix[[nam]]
  ada <- boost.mod.cv[[nam]]
  
  #- fix
  n <- length(fix) - 1
  boost_all   <- fix[[n + 1]]   
  boost_drop1 <- fix[-(n + 1)]
  
  if (method == "M1"){
    fix <- mapply(ham, sj = boost_drop1, MoreArgs = list(si = boost_all))
    fix <- scale(fix)
  }
  
  if (method == "M2"){
    fix <- moddist(mod = boost_drop1, std = F, L = 1)
    fix <- fix / sd(fix)
  }
  
  #- ada
  boost_all   <- ada[[n + 1]]   
  boost_drop1 <- ada[-(n + 1)]
  
  if (method == "M1"){
    ada <- mapply(ham, sj = boost_drop1, MoreArgs = list(si = boost_all))
    ada <- scale(ada)
  }
  
  if (method == "M2"){
    ada <- moddist(mod = boost_drop1, std = F, L = 1)
    ada <- ada / sd(ada)
  }
  
  ind <- c(order(fix, decreasing = T)[1:3], order(ada, decreasing = T)[1:3])
  
  if (plot == "cor"){
    corrplot(fix, ada, ci.line = F, xlab = "Fix", ylab = "Ada",
             xlim = c(min(c(fix, ada)), max(c(fix, ada)) * 1.1),
             ylim = c(min(c(fix, ada)), max(c(fix, ada)) * 1.1))
    text(fix[ind], ada[ind], labels = ind, pos = 3)
  }
  
  
  if (plot == "BA"){
    blandartman(fix, ada, ci = F, spline = F, ci.limit = F, 
                xlab = "(Fix + Ada) / 2", ylab = "Fix - Ada",
                ylim = c(min(fix - ada) * 1.1, max(fix - ada) * 1.1),
                xlim = c(min(fix/2 + ada/2), max(fix/2 + ada/2) * 1.1))
    text((fix/2 + ada/2)[ind], (fix - ada)[ind], labels = ind, pos = 3)
  }
}

jpeg("./output/BA.in.tuning.jpeg", width = 9, height = 6, res = 300, units = "in")
par(mfcol = c(2, 3), mar = c(3.5, 3.5, 1, 1), oma = c(1, 1, 1, 1))

for(i in c(4, 5, 8)){
  nam <- d[i]
  plotm(nam, plot = c("cor"), method = "M2")
  mtext(gsub(".RData", "", nam), side = 3, line = 0.2)
  plotm(nam, plot = c("BA"),  method = "M2")
}

dev.off()

#-------------------------------------------------------------------------------











































