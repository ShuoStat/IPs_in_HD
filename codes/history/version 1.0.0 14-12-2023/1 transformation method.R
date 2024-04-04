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

#- 50 if optimal mstops is too small
mstopSet <- c(818, 3230, 1272, 136, 262, 218)
nu       <- c(0.1, 0.3, 0.3, 0.1, 0.1, 0.1)

#-------------------------------------------------------------------------------

transMod.ada <- list()
transMod.fix <- list()

for (i in 1:6) {
  
  print(i)
  load(paste0("../data/", d[i]))
  n = nrow(X)
  X <- t(X)
  family = type[[i]]
  
  if (link[i] == "binomial") 
    y = as.factor(y)
  if (link[i] == "cox")
    y = surv
  
  # get cv
  set.seed(1)
  nfolds = 10
  foldid <- sample(rep(1:nfolds, length = nrow(X)), nrow(X)) 
  cv <- sapply(1:max(foldid), function(x) as.numeric(foldid != x))
  
  transMod.ada[[d[i]]] <- transregress(X, y, family = family, 
                                       nu = nu[i],
                                       mstop = mstopSet[i], 
                                       cv = cv,
                                       adaptive = T)
  
  transMod.fix[[d[i]]] <- transregress(X, y, family = family, 
                                       nu = nu[i],
                                       mstop = mstopSet[i], 
                                       cv = cv,
                                       adaptive = F)
}

save(list = c("transMod.ada", "transMod.fix"), file = "../output/transMod.RData")

#- table -----------------------------------------------------------------------

load("../output/transMod.RData")

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
tab <- lapply(transMod.ada, getset)
tab <- cbind.data.frame(tab)
write.csv(t(tab), file = "../results/transMod.ada.csv")

#- fixed
tab <- lapply(transMod.fix, getset)
tab <- crbind.data.frame(tab)
write.csv(t(tab), file = "../results/transMod.fix.csv")

#-------------------------------------------------------------------------------
#- Leverage
#-------------------------------------------------------------------------------

load("../output/transMod.RData")

fix = T

if (fix == T){
  obj = transMod.fix
} else {
  obj = transMod.ada
}

jpeg(paste0("../results/lev.", ifelse(fix, "fix", "ada"), ".jpeg"),
     width = 9, height = 6, res = 300, units = "in")

par(mfrow = c(3, 3), mar = c(3.5, 3.5, 2, 1.5), oma = c(1, 3, 1, 1))

for(i in c(3, 4, 6)){
  
  hat.ada.x <- diag(obj[[i]]$hatX)
  hat.ada.z <- diag(obj[[i]]$hatZ)
  
  lollipop(hat.ada.x, ylim = c(0, max(hat.ada.x)* 1.2), ylab = "Leverage")
  mtext(text = paste0(gsub(".RData", "", d[i])), side = 2, line = 4, cex = 1.2)
  mtext(text = "Before transformation", side = 3, line = 0.2,  adj = 0)
  
  lollipop(hat.ada.z, ylim = c(0, max(hat.ada.x)* 1.2), ylab = "Leverage")
  mtext(text = "After transformation", side = 3, line = 0.2,  adj = 0)
  
  xlim = c(0, max(hat.ada.x + hat.ada.z) / 2 * 1.2)
  ylim = max(abs(hat.ada.x - hat.ada.z)) * c(-1.2, 1.2)
  blandartman(hat.ada.x, hat.ada.z, spline = F, ci = T, ci.limit = F,
              xlab = "lev(OD)/2 + lev(TD)/2",
              ylab = "lev(OD) - lev(TD)",
              xlim = xlim,
              ylim = ylim)
  
  blandartmanAddmarks(hat.ada.x, hat.ada.z, pos = 4, cex = 0.8)
}
  
dev.off()

#-------------------------------------------------------------------------------












