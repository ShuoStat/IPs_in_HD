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

#- Univariate analysis ---------------------------------------------------------
#- output all results

fun <- c("t.test", "t.test", "t.test", "t.test", "t.test", "t.test", 
         "cox", "cox")

#- summary all results
re.transrank <- list()
for (i in 1:length(d)){
  load(paste0(dir, d[i]))
  if (fun[i] == "cox")
    y = surv
  re.transrank[[d[i]]] <- transrank(X, y, range = 100, eps = 0.01, 
                                    scale = TRUE, FUN = fun[i],
                                    plot = F, cut = c(-5, -2, -1, 1, 2, 5))
}
#- save(list = "re.transrank", file = "./output/re.transrank.RData")


#- how many 100+ ---------------------------------------------------------------

load("./output/re.transrank.RData")
lapply(re.transrank, nrow)

#--plot all --------------------------------------------------------------------

plotfun.addtop10 <- function(obj, range = 100, topn = 10, 
                             labels = c("weighted", "unweighted"),
                             method = c("weighted", "unweighted"),
                             labels.topn = 1){
  
  rank1 <- obj$rank1
  rank2 <- obj$rank2
  
  ord <- order(abs(rank1 - rank2), decreasing = T)[1:topn]
  rank1[rank1 > range] <- range
  rank2[rank2 > range] <- range
  
  if ("unweighted" %in% method){
    points(rank1[ord], rank2[ord], col = "red", pch = 19)
  }

  if ("weighted" %in% method){
    ord1 <- order(abs(obj$delta), decreasing = T)[1:topn]
    points(rank1[ord1], rank2[ord1], col = "blue", pch = 19)
  }
  
  if (setequal(method, c("weighted", "unweighted"))){
    s <- intersect(ord1, ord)
    points(rank1[s], rank2[s], col = "green", pch = 19)
  }
  
  if ("unweighted" %in% labels) {
    ordt <- ord[1:labels.topn]
    text(rank1[ordt], rank2[ordt], labels = obj$genes[ordt], 
         pos = 2, cex = 0.9)
  }
  
  if ("weighted" %in% labels) {
    ord1t <- ord1t[1:labels.topn]
    text(rank1[ord1t], rank2[ord1t], labels = obj$genes[ord1t], 
         pos = 2, cex = 0.9)
  }
}

#- Boxplot of two groups

pbox <- function (obj, X, name, topn = NULL, ylim1 = NULL, ylim2 = NULL,
                  weighted = F, col = "red"){
  
  X <- apply(X, 2, scale)
  rownames(X) <- 1:nrow(X)
  
  rank1 <- obj$rank1
  rank2 <- obj$rank2
  
  d_r <- abs(rank1 - rank2)
  
  if (!weighted)
    g <- order(d_r, decreasing = T)[1:topn]
  else
    g <- order(abs(obj$delta), decreasing = T)[1:topn]
  g <- obj$genes[g]
  
  x1 <- X[y == levels(y)[1], g, drop = F]
  x2 <- X[y == levels(y)[2], g, drop = F]
  colnames(x1) <- colnames(x2) <- g
  
  par(mar = c(0, 4, 3, 1))
  if(is.null(ylim1)) ylim1 = c(min(x1, x2), max(x1, x2))
  boxplot(x1, xaxt = "n", las = 2, col = "white", ylim = ylim1)
  
  mapply(function(i){
    
    x_i <- x1[,i]
    if (outliers::grubbs.test(x_i)$p.value < 0.001) {
      outs = 100000 #- the initial drops
      p.value = 0
      p.val = c() # collect the p values
      while(p.value < 0.001){
        p.value <- outliers::grubbs.test(x_i[!names(x_i) %in% outs])$p.value
        if (p.value < 0.001){
          outcase <- as.numeric(names(p.value))
          outs <- c(outs, outcase)
          p.val <- c(p.val, p.value)
        }
      }
      
      outs <- outs[-1] #- remove 100000
#     outs.01 <- outs[p.val >= 0.001]
#     points(rep(i, length(outs.01)), x_i[match(outs.01, names(x_i))], pch = 19, 
#            col = "black", cex = 0.7)
      
      outs.001 <- outs[p.val < 0.001]
      points(rep(i, length(outs.001)), x_i[match(outs.001, names(x_i))], pch = 19, 
             col = col, cex = 0.7)
      
      out.txt <- outs[1:2]
      text(rep(i, length(out.txt)), x_i[match(out.txt, names(x_i))], 
           labels = out.txt, pos = 4, cex = 0.9)
    }
  }, i = 1: ncol(x1))
  
  mtext("Group 1", side = 2, line = 2.3)
  mtext(name, line = 0.2, side = 3, adj = 0)
  
  par(mar = c(3, 4, 0, 1))
  if(is.null(ylim2)) ylim2 = c(min(x1), max(x2))
  boxplot(x2, las = 2, col = "white", ylim = ylim2)
  mtext("Group 2", side = 2, line = 2.3)
  
  mapply(function(i){
    
    x_i <- x2[,i]
    if (outliers::grubbs.test(x_i)$p.value < 0.001) {
      
      outs = 100000 #- the initial drops
      p.value = 0
      p.val = c()   #- collect the p values
      while(p.value < 0.001){
        p.value <- outliers::grubbs.test(x_i[!names(x_i) %in% outs])$p.value
        if (p.value < 0.001){
          outcase <- as.numeric(names(p.value))
          outs <- c(outs, outcase)
          p.val <- c(p.val, p.value)
        }
      }
      
      outs <- outs[-1] 
      #- remove 100000
      #- outs.01 <- outs[p.val >= 0.001]
      #- points(rep(i, length(outs.01)), x_i[match(outs.01, names(x_i))], pch = 19, 
      #-        col = "black", cex = 0.7)
      
      outs.001 <- outs[p.val < 0.001]
      points(rep(i, length(outs.001)), x_i[match(outs.001, names(x_i))], pch = 19, 
             col = col, cex = 0.7)
      
      out.txt <- outs[1:2]
      text(rep(i, length(out.txt)), x_i[match(out.txt, names(x_i))], 
           labels = out.txt, pos = 4, cex = 0.9)
    }
  }, i = 1: ncol(x2))
}

obj <- re.transrank[[s]]

pbox.surv <- function(obj, X, name, topn = NULL, ylim, weighted = F, col = "red") {
  
  rk <- obj

  #- boxplot
  X <- apply(X, 2, scale)
  rownames(X) <- 1:nrow(X)
  rank1 <- rk$rank1
  rank2 <- rk$rank2
  
  d_r <- abs(rank1 - rank2)

  if (!weighted)
    g <- order(d_r, decreasing = T)[1:topn]
  else
    g <- order(abs(obj$delta), decreasing = T)[1:topn]
  g <- obj$genes[g]
  
  x <- X[, g]
  colnames(x) <- g
  boxplot(x, las = 2, ylim = ylim)
  
  pfun <- function(i, col){
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
      outs <- outs[-1] 
      
      outs.001 <- outs[p.val < 0.001]
      points(rep(i, length(outs.001)), x_i[match(outs.001, names(x_i))], pch = 19, 
             col = col, cex = 0.9)
      
      out.txt <- outs[1:2]
      text(rep(i, length(out.txt)), x_i[match(out.txt, names(x_i))], 
           labels = out.txt, pos = 4, cex = 0.7)
    }
  }
  
  mapply(pfun, i = 1: ncol(x), col = col)
  mtext(name, line = 0.2, side = 3, adj = 0)
}

#- get the matrix for layout
get.layoutmatrix <- function(nrow, ncol){
  
  m <- m0 <- matrix(1: (5 * ncol), ncol = ncol, byrow = F)
  if (nrow > 1)
    for(i in 2:nrow){
      add <- (i - 1) * ncol * 3
      m <- rbind(m, m0 + add)
    }
  return(m)
}

get.height <- function(nrow) rep(c(2.5, 1.5, 1.5, 1.5, 1.5), nrow)

#- plot all data of binary outcomes

load("./output/re.transrank.RData")

#- main text figure
s = c( "petricoin02.RData", "golub99.RData", "veer02.RData")
ylim <- list(c(-3, 6), c(-1, 7), c(-6, 8))

s = c("singh02.RData", "ressom06.RData","wang05.RData")
#- appendix figure
ylim <- list(c(-2, 8), c(-2, 8), c(-3, 8))

jpeg("./output/rank-box-all.jpeg", width = 10, height = 8, 
     units = "in", res = 300)

mat <- get.layoutmatrix(1, 3)
if ("veer02.RData" %in% s)
  mat[2:5, 3] <- c(12, 12, 13, 13)

layout(mat, heights = get.height(1))
par(oma = c(1, 1, 1, 1))

for(i in 1:3){
  
  par(mar = c(3, 4, 2.5, 1))
  rk <- re.transrank[[s[i]]]
  plotfun(rk$rank1, rk$rank2, range = 100, labels = rk$genes, n.labels = 1)
  plotfun.addtop10(rk, range = 100, topn = 10, 
                   labels = T, labels.topn = 1)
  mtext(paste0(letters[i], "), ", gsub(".RData", "", s[i])), 
        side = 3, line = 0.2, adj = 0)
  load(paste0(dir, s[i]))
  X <- t(X)
  
  if (s[i] != "veer02.RData" ) {
    pbox(rk, X = X, name = paste0(gsub(".RData", "", s[i]), "(RD)"), 
         topn = 10, weighted = F, col = "red", ylim1 = ylim[[i]],
         ylim2 = ylim[[i]])
    
    pbox(rk, X = X, name = paste0(gsub(".RData", "", s[i]), "(RD-W)"), 
         topn = 10, weighted = T, col = "blue", ylim1 = ylim[[i]],
         ylim2 = ylim[[i]])
  } else{
    pbox.surv(rk, X = X, name = paste0(gsub(".RData", "", s[i]), "(RD)"), 
         topn = 10, weighted = F, col = "red", ylim = ylim[[i]])
    
    pbox.surv(rk, X = X, name = paste0(gsub(".RData", "", s[i]), "(RD-W)"), 
         topn = 10, weighted = T, col = "blue", ylim = ylim[[i]])
  }
}

dev.off()

#- correlation matrix in ressom06 ----------------------------------------------

i = "ressom06.RData"
obj <- re.transrank[[i]]
rank1 <- obj$rank1
rank2 <- obj$rank2
#- absolute changes
# g <- abs(rk$rank1 - rk$rank2)
# g <- rk$genes[order(g, decreasing = T)[1:topn]]

#- absolute changes
topn <- 10
g <- rk$genes[order(abs(obj$delta), decreasing = T)[1:topn]]

load(paste0(dir, i))
X <- apply(X, 1, scale)
x <- X[, g]

cormat <- round(cor(x), 3)
cormat[upper.tri(cormat)] <- ""
colnames(cormat) <- rownames(cormat) <- g

#-- scate plot for ressom06 ----------------------------------------------------

i = "ressom06.RData"
load(paste0(dir, i))
#- 2679 and 2677； 9572 and 9571
ressom <- function(x, y){
  
  xx <- scale(X[x, ])
  yy <- scale(X[y, ])
  
  plot.new()
  plot.window(xlim = range(c(xx, yy)), ylim = range(c(xx, yy)))
  axis(2)
  axis(1)
  box()
  lines(c(-10^10, 10^10), c(-10^10, 10^10))
  
  points(xx, yy)
  legend("topleft", legend = paste0("r = ", round(cor(xx, yy), 3)), border = F)
  z <- lm(yy ~ xx)
  lines(z$model$xx, z$fitted.values, lty = "dashed")
  mtext(x, side = 1, line = 2.2)
  mtext(y, side = 2, line = 2.2)
}

# 2679 and 2677； 9572 and 9571
jpeg("./output/ressom.cor.jpeg", width = 8, height = 4, res = 300, units = "in")
par(mfrow = c(1, 2), mar = c(3.5, 3.5, 1, 1))
ressom(2679, 2677)
ressom(9572, 9571)
dev.off()



write.csv(cormat, "./output/tmp.csv")

#-------------------------------------------------------------------------------



#- grubbs for two points test -#
#- wang05, 11110; golub99 2505; ressom06 9574; scherzer07 22256
grubbs.2p <- function(s, feature){
  load(paste0(dir, s))
  X <- apply(X, 1, scale)
  rownames(X) <- 1:nrow(X)
  x = X[y == levels(y)[1],feature]
  outliers::grubbs.test(x, type = 20)
}

grubbs.2p("wang05.RData", 11110)
grubbs.2p("golub99.RData", 2505)
grubbs.2p("ressom06.RData", 9574)
grubbs.2p("scherzer07.RData", 22256)




















