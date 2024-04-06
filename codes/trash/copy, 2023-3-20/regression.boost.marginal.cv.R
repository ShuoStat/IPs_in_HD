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
#- marginal correlations
#-------------------------------------------------------------------------------

drop1.margin <- function(s, family){
  
  load(paste0(dir, s))
  X <- t(X)
  X <- scale(X)
  
  n = nrow(X)
  p = ncol(X)
  
  if (family == "cox")
    y <- surv 
  if (family == "binomial")
    y <- as.factor(y)
  
  library(doParallel)
  ncors <- detectCores() - 2
  cl <- makeCluster(ncors)
  registerDoParallel(cl)
  
  if (family == "binomial"){
    beta <- foreach(i = 1:(n + 1), 
                    .combine = "cbind") %dopar% {
                      beta <- apply(X[-i,], 2, function(x){
                        glm.fit(x, y[-i], family = binomial())$coefficients
                      })
                      return(beta)
                    }
  }
  
  if (family == "cox"){
    beta <- foreach(i = 1:(n + 1), .export = "survival",
                    .combine = "cbind") %dopar% {
                      beta <- apply(X[-i,], 2, function(x){
                        coxph(y[-i] ~ x)$coefficients
                      })
                      return(beta)
                    }
  }
  
  registerDoSEQ()
  stopCluster(cl)
  return(beta)
}

f <- c("binomial", "binomial", "binomial", "binomial", 
       "binomial", "binomial", "cox", "cox")
re.drop1.margin <- list()
for(i in 1:8){
  s = d[i]
  family = f[i]
  re.drop1.margin[[s]] <- drop1.margin(s = s, family = family)
}
#- save(list = "re.drop1.margin", file = "./output/re.drop1.margin.RData")
#-------------------------------------------------------------------------------

load("./output/re.drop1.margin.RData")

jpeg("./output/drop1margin.main.jpeg", width = 9, height = 2.8, units = "in", 
     res = 300)
par(mar = c(3.5, 3.5, 2, 1.5), mfrow = c(1, 3))
for(i in d[c(4, 5, 8)]){
  
  beta = re.drop1.margin[[i]]
  n <- ncol(beta) - 1
  p <- nrow(beta)
  score <- apply(beta[,-(n + 1)], 2, function(x) 
    sum((x - beta[,(n + 1)])^2) / p)
  score <- score / sd(score)
  names(score) <- paste0("obs", 1:n)
  
  lollipop(score, decreasing = T, topn = 5, ylim = c(0, 12), refline = 0, 
           ylab = "HIM")
  mtext(gsub(".RData", "", i), line = 0.2, side = 3, adj = 0)
}
dev.off()


jpeg("./output/drop1margin.append.jpeg", width = 9, height = 3, units = "in", 
     res = 300)
par(mar = c(3.5, 3.5, 2, 1.5), mfrow = c(1, 3))
for(i in d[c(2, 3, 6)]){
  beta = re.drop1.margin[[i]]
  n <- ncol(beta) - 1
  p <- nrow(beta)
  score <- apply(beta[,-(n + 1)], 2, function(x) 
    sum((x - beta[,(n + 1)])^2) / p)
  score <- score / sd(score)
  names(score) <- paste0("obs", 1:n)
  lollipop(score, decreasing = T, topn = 5, ylim = c(0, 15), refline = 0, 
           ylab = "HIM")
  mtext(gsub(".RData", "", i), line = 0.2, side = 3, adj = 0)
}
dev.off()

#-------------------------------------------------------------------------------


