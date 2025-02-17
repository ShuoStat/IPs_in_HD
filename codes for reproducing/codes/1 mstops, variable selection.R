
# load help function
source("./helper.R")
# loading packages 
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

# prior mstop of each data
# double of the optimized mstops in original model
mstopSet <- c(818, 3230, 1272, 136, 262, 218)
# a large nu for the data with large mstops
nu      <- c(0.1, 0.3, 0.3, 0.1, 0.1, 0.1)
# Note, large RAM is needed, therefore, consider balance the CPU cores and 
# RAM 
# The user need to define the Ncores based on the computer.
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
                         time2 <- Sys.time()
                         
                         cvboost <- cvrisk(boost, folds = cv[-k,])
                         out <- mstop(cvboost)
                         time3 <- Sys.time()
                         
                         out <- list(cvboost = cvboost,
                                     "id" = k,
                                     "mstop" = out, 
                                     "timeFit" = as.numeric(
                                       difftime(time2, 
                                                time1, 
                                                unit = "secs")),
                                     "timeCV" = as.numeric(
                                       difftime(time3, 
                                                time2, 
                                                unit = "secs")))
                         
                         return(out)
                       }
  
  mstops <- list(mstops = re.mstops, foldid = foldid)
  save(list = "mstops", file = paste0("../output/mstop_", d[i]))
  
  registerDoSEQ()
  stopCluster(cl)
}

#-------------------------------------------------------------------------------
#- Computation Cost
#-------------------------------------------------------------------------------

# timeMat, the total time used for model fit process(columns Fit) and tuning 
# parameter selection process(CV columns)

timeMat <- matrix(NA, length(d), 2,
                  dimnames = list(gsub(".RData", "", d),
                                  c("Fit", "CVs")))
nSample <- c()

for(i in d) {
  
  load(paste0("../output/mstop_", i))
  obj = mstops$mstops
  timeFit <- unlist(lapply(obj, `[[`, "timeFit"))
  timeCV <- unlist(lapply(obj, `[[`, "timeCV"))
  timeMat[gsub(".RData", "", i), ] <- c(sum(timeFit), sum(timeCV))
  nSample <- c(nSample, length(timeFit) - 1)
}

# units, from seconds to minutes
timeMat <- timeMat / 60 

# visualize time cost on a barplot
jpeg("../results/time.jpeg", width = 6, height = 3, units = "in", res = 300)

par(mar = c(4, 6, 1, 9), cex = 0.5)
# barplot, time cost 
barplot(t(timeMat), 
        horiz = TRUE, 
        xlim = c(0, max(timeMat) * 1.2),
        las = 2, 
        legend.text = c("Model Fiting Procedure", "CV Procedure"),
        col = c("white", "gray"),   
        cex.axis = 0.7, 
        xaxt = "n")

# add sample size information
mtext(nSample,
      at = seq(0.75, length(nSample) + 0.75, length = length(nSample)),
      side = 4,
      las = 2,
      line = 0.5,
      cex = 0.5)

# add number of iterations
mtext(mstopSet,
      at = seq(0.75, length(mstopSet) + 0.75, length = length(mstopSet)),
      side = 4,
      las = 2,
      line = 4.5,
      cex = 0.5)

# add columns name for sample size and iterations
mtext(c("N(sample)", "N(iterations)"),
      at = length(mstopSet) + 1.5,
      side = 4,
      las = 2,
      line = c(0.5, 4.5),
      cex = 0.5)


axis(1)
mtext("Time, min", side = 1, line = 2, cex = 0.7)
box()

dev.off()

#-------------------------------------------------------------------------------
#- IPs on variable selection
#-------------------------------------------------------------------------------

drop1vs <- function(datName, 
                    family, 
                    nu, 
                    method = c("fixed", "adaptive"), 
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
                                         control = boost_control(
                                           mstop = mstops[i], 
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

load("../output/vs.RData")

# Figure 1, ada_m1vsm2 
jpeg("../results/ada_m1vsm2.jpeg", width = 8, height = 7, units = "in", res = 300)

mat <- matrix(1:9, nrow = 3, byrow = F)
layout(mat)
par(mar = c(3.5, 3.5, 1.5, 1.5), oma = c(1, 1, 1, 1))
for(s in d[c(3, 4, 6)]){
  
  obj <- vsAda[[s]]
  
  plot_vs(obj,  method = "M1",
          name = paste0(gsub(".RData", "", s), "(adaptive)"),
          xlim = c(0, 12), plot = "lol", center = F)
  
  plot_vs(obj,  method = "M2",
          name = "",
          xlim = c(0, 12), plot = "lol", center = F)
  
  plot_m1vsm2(obj, plot.type = "c", name = "")
}

dev.off()

# ada_m1vsm2 in Appendix

jpeg("../results/ada_m1vsm2_Appendix.jpeg", width = 8, height = 7, units = "in", res = 300)

mat <- matrix(1:9, nrow = 3, byrow = F)
layout(mat)
par(mar = c(3.5, 3.5, 1.5, 1.5), oma = c(1, 1, 1, 1))

for(s in d[-c(3, 4, 6)]){
  obj <- vsAda[[s]]
  plot_vs(obj,  method = "M1",
          name = paste0(gsub(".RData", "", s), "(adaptive)"),
          xlim = c(0, 12), plot = "lol", center = F)
  
  plot_vs(obj,  method = "M2",
          name = paste0(gsub(".RData", "", s), "(adaptive)"),
          xlim = c(0, 12), plot = "lol", center = F)
  
  plot_m1vsm2(obj, plot.type = "c", name = "")
}

dev.off()

#-------------------------------------------------------------------------------
#- Plot: golub99 Data, show the difference between M1 and M2
#-------------------------------------------------------------------------------

load("../output/vs.RData")

m1_m2_Plot <- function(data, dropID, type = c("M1", "M2"), title = "") {
  
  obj <- vsAda[[paste0(data, ".RData")]]
  n = length(obj) - 1
  boost_all   <- obj[[n + 1]]   
  boost_drop1 <- obj[-(n + 1)]
  
  # calculate the selection frequency
  n <- length(boost_drop1)
  allvar <- unique(unlist(boost_drop1))
  f <- table(unlist(boost_drop1))[allvar]
  p <-  f / n  # inclusion frequency
  p <-  sort(p, decreasing = TRUE)
  
  x <- seq_along(p)
  y <- as.numeric(p)
  labels <- names(p)
  
  p0 <- boost_all
  p1 <- boost_drop1[[dropID]]
  
  if (type == "M1") {
    # M1 
    plot.new()
    plot.window(xlim = c(0, length(y)), ylim = c(0, 1))
    axis(1, xaxt = "n"); axis(2, at = c(0, 1), tick = c(0, 1), labels = c(0, 1)); box()
    # points(x, y, pch = 20)
    abline(h = c(0, 1), lty = "dashed")
    
    y0 <- as.numeric(labels %in% p0)
    y1 <- as.numeric(labels %in% p1)
    x <- seq_along(labels)
    # original model
    points(x, y0, col = "black", pch = 1, cex = 0.9)
    # p32 model
    points(x, y1, col = "red", pch = 4, cex = 0.9)
    
    for(i in seq_along(y1)){
      lines(c(i, i), c(y1[i], y0[i]), col = "black", pch = 20)
    }
    
    legend("topright", pch = c(1, 4), col = c("black", "red"), legend = c("orig", paste0("-obs", dropID)), inset = 0.03, bty = "n")
    mtext(title, side = 3, line = 0.2, adj = 0)
    mtext("Variable Selection Status", side = 2, line = 2.2)
    
  }
  # M2
  
  if (type == "M2") {
    plot.new()
    plot.window(xlim = c(0, length(y)), ylim = c(0, 1))
    axis(1, xaxt = "n"); axis(2); box()
    points(x, y, pch = 1, cex = 0.9)
    abline(h = c(0, 1), lty = "dashed")
    
    x1 <- match(p1, labels)
    x0  <- match(p0, labels)
    
    points(x, y1, col = "red", pch = 4, cex = 0.9)
    
    y1 <- as.numeric(labels %in% p1)
    for(i in seq_along(y1)){
      lines(c(i, i), c(y1[i], y[i]), col = "black")
    }
    
    legend("topright", pch = c(1, 4), col = c("black", "red"), legend = c("IF", paste0("-obs", dropID)), inset = 0.03, bty = "n")
    mtext(title, side = 3, line = 0.2, adj = 0)
    mtext("Inclusion Frequency", side = 2, line = 2.2)
  }
}


jpeg("../results/m1_m2_compare.jpeg", width = 10, height = 6, units = "in", res = 300)
par(mfrow = c(2, 2), mar = c(2, 4, 1, 1), oma = c(1, 1, 1, 1))
# layout(mat = matrix(c(1, 2, 3, 4), byrow = F, nrow = 2))

m1_m2_Plot(data = "golub99", dropID = 32, type = c("M1"), title = "A, M1 (-obs32)")
m1_m2_Plot(data = "golub99", dropID = 32, type = c("M2"), title = "B, M2 (-obs32)")

m1_m2_Plot(data = "golub99", dropID = 69, type = c("M1"), title = "C, M1 (-obs69)")
m1_m2_Plot(data = "golub99", dropID = 69, type = c("M2"), title = "D, M2 (-obs69)")

dev.off()

# Compare Fixed and adaptive models based on M2 ---------------------------------

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
  
  vsPathlist[[d[i]]] <- vs_path(X, y, family, obs.mstop = mstops, nu[i], nCores = 4)
}

# save(list = "vsPathlist", file = "../output/vsPathlist.RData")
load("../output/vsPathlist.RData")
load("../output/vs.RData")

# based on M1
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

# based M2

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

#-------------------------------------------------------------------------------
# 
# 
# # Compare Fixed and adaptive models based on M1 ---------------------------------
# 
# vsPathlist <- list()
# 
# for(i in c(3, 4, 6)){
#   
#   print(d[i])
#   load(paste0("../output/mstop_", d[i]))
#   load(paste0("../data/", d[i]))
#   X <- t(X)
#   
#   family = type[[i]]
#   mstops = unlist(lapply(mstops$mstops, `[[`, "mstop"))
#   
#   if (class(family) == "boost_family_glm")
#     y = as.factor(y)
#   if (class(family) == "boost_family")
#     y = surv
#   
#   vsPathlist[[d[i]]] <- vs_path(X, y, family, obs.mstop = mstops, nu[i], nCores = 20, method = "M1")
# }
# 
# # save(list = "vsPathlist", file = "../output/vsPathlist_M1.RData")
# load("../output/vsPathlist_M1.RData")
# load("../output/vs.RData")
# 
# # based on M1
# jpeg("../results/vsPath_M1.jpeg", width = 8, height = 5, units = "in", res = 300)
# par(mfcol = c(2, 3), mar = c(3.5, 3.5, 2, 2))
# 
# for(i in names(vsPathlist)){
#   
#   # pathlist
#   plot.vs_path(vsPathlist[[i]], ylab = "M1")
#   mtext(gsub(".RData", "", i), side = 3, line = 0.2, adj = 0)
#   
#   # corplot to compare fixed and adaptive
#   fixObj <- vsFixed[[i]]
#   adaObj <- vsAda[[i]]
#   
#   plot_fixvsada(fixObj, adaObj, method = "M1", name = gsub(".RData", "", i))
# }
# 
# dev.off()
# 


#-------------------------------------------------------------------------------

