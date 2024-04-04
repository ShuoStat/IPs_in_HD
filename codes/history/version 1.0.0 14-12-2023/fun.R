
#' @title  Data transformation
#' @description
#' transdata is a data transformation method aiming at coping with extreme values.
#' This programs can automatically transform vectors, matrix, data.frame, and list
#' and return the transformed data in the same type of data. Note that, only continuous
#' data can be transformed. Missing data and categorical variables will be kept the same
#' without any tranformation.
#'
#' @inheritParams replace.value
#' @return This function return the transformed data with the same format.
#' @param x, a vectors, matrix, data.frame or list.
#' @param eps, parameter adjusting the curve shape between original data and transformed data.
#' With larger epi, extreme values are more concentrated or shrunk after transformation.
#' The default value for epi is 0.01.
#' @param scale, a logical value. Setting scale  = TRUE will scale transformed values to [0, 1]. Default as TRUE
#' @param byrow, a logical value. For data.frame and matrix, Setting byrow = TRUE will tranform data by rows. Default as FALSE.
#' @param formula, a formula. For data.frame and matrix, response and categorical variables (e.g. factor(x)) indicated by a formula will not be transformed.
#' @param categorical, a vector of variables. For data.frame and matrixs, categorical indicates the categorical variables which do not need transformation.
#' @...
#' @seealso rankfun, plotrank
#' @author Shuo Wang
#' @details The data transformation strategy was proposed by Royston and Sauerbrier(2007)
#' to reduce the leverage of large values. ...
#'
#'
#' @export
#' @examples
#'
#'cars[,"group"] <- sample(c(1, 2), nrow(cars), replace = TRUE)
#' #1. effects of eps
#'original <- x$dist
#'tran1 <- transdata(original, eps = 0.01)
#'tran2 <- transdata(original, eps = 0.10)
#'tran3 <- transdata(original, eps = 0.20)
#'tran4 <- transdata(original, eps = 0.50)
#'
#'par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
#'plot(original, tran1, main = "eps = 0.01",
#'     xlab = "Original", ylab = "Transformed")
#'plot(original, tran2, main = "eps = 0.10",
#'     xlab = "Original", ylab = "Transformed")
#'plot(original, tran3, main = "eps = 0.20",
#'     xlab = "Original", ylab = "Transformed")
#'plot(original, tran4, main = "eps = 0.50",
#'     xlab = "Original", ylab = "Transformed")
#'
#' #2. transformation using categorical
#' trandat <- transdata(x = cars, eps = 0.01, scale = F,
#'                     categorical = "group")
#'head(trandat)
#' #3. transformation using formula
#'formula1 <- as.formula(dist~speed + factor(group))
#'trandat <- transdata(x = cars, eps = 0.01, scale = F,
#'                    formula = formula1)
#'
#'

#-------------------------------------------------------------------------------

transdata <- function(x, eps = 0.01, scale = TRUE, byrow = FALSE,
                      formula = NULL, categorical = NULL,
                      ...){
  
  arg <- mget(names(formals()))
  if (!is.logical(scale)) stop("Scale should be logical")
  if (!is.logical(byrow)) stop("Byrow should be logical")
  
  type <- class(x)[1] #use the first class type
  if (!type %in% c('matrix', 'data.frame', 'list', 'other')) type = 'other'
  
  if ((!type %in% c("matrix", "data.frame")) & byrow == TRUE)
    warning("byrow = TURE only work for matrix and data.frame")
  
  # transformation for matrix
  mfun <- function(x, eps, scale, byrow, ...){
    margin    <- 2 - as.numeric(byrow)
    newmatrix <- apply(x, margin, function(x) transfo(x, eps = eps,
                                                      scale = scale))
    if (byrow) return(t(newmatrix)) else
      return(newmatrix)
  }
  
  # transformation for data.frame
  dfun <- function(x, eps, scale, byrow, ...){
    
    if(byrow == FALSE){
      newframe <- lapply(x, function(x) transfo(x, eps = eps,
                                                scale = scale))
    }
    
    if(byrow == TRUE){
      x <-  as.data.frame(t(x))
      newframe <- lapply(x, function(x) transfo(x, eps = eps,
                                                scale = scale))
      
      #- if characters rows existed, all rows would be character vectors
      #- In such situation, no transformation is returned
      isnumeric <- lapply(x, is.numeric)
      if(all(unlist(isnumeric)) == FALSE)
        warning("No transformation, because of character rows. \n\r byrow = T is not recommended in data.frame.")
    }
    
    newframe <- as.data.frame(newframe)
    if (byrow) return(t(newframe)) else
      return(newframe)
  }
  
  # transformation for list
  lfun <- function(x, eps, scale, ...){
    newlist <- lapply(x, function(x) transfo(x, eps = eps,
                                             scale = scale))
    return(newlist)
  }
  
  transd <- switch(type,
                   'matrix'     = do.call(mfun, arg),
                   'data.frame' = do.call(dfun, arg),
                   'list'       = do.call(lfun, arg),
                   'other'      = do.call(transfo, arg)
  )
  
  # using formula to guide data transformation
  if (!is.null(formula)) {
    if (!byrow & is.data.frame(x) & inherits(formula, "formula")) {
      modf <- model.frame(formula, x)
      attrmod <- attr(modf, "terms")
      factx   <- match("factor", attr(attrmod, "dataClasses"))
      allvar  <- all.vars(formula)
      untranvar <- allvar[c(1, factx)]
      untranvar <- untranvar[!is.na(untranvar)]
      transd[, untranvar] <- x[, untranvar]
      transd <- transd[, allvar]
    } else{
      stop("Checking the following condition: \n
                 1, byrow = TURE is not allowed\n
                 2, x should be data.frame\n
                 3, whether formula is correctly defined")
    }
  }
  
  if (!is.null(categorical)) {
    if (byrow)
    { if (!categorical %in% rownames(x))
      stop("Variable is not found in data") else
        transd[categorical, ] <- x[categorical, ]}
    if (!byrow)
    {if (!categorical %in% colnames(x))
      stop("Variable is not found in data") else
        transd[, categorical] <- x[, categorical]}
  }
  return(transd)
}

#- transformation method

transfo <- function(x, eps, scale, ...) {
  
  if (is.character(x) | is.factor(x)) stop(return(x))
  
  sdx   <- sd(x, na.rm = T)
  meanx <- mean(x, na.rm = T)
  z     <- (x - meanx) / sdx
  eps2  <- log((1 + eps) / eps)
  
  phiz <- pnorm(z)
  g    <- (log((phiz + eps) / (1 - phiz + eps)) + eps2) / (2*eps2)
  if (scale) {
    return(g)
  } else {
    return(scale(g)*sdx + meanx)
  }
}

#- loading packages

load_packages <- function(packages){
  
  # Check if packages are installed
  not_installed <- packages[!(packages %in% installed.packages()[,"Package"])]
  
  # Install packages if they are not already installed
  if (length(not_installed) > 0) {
    install.packages(not_installed)
  }
  
  # Load packages
  lapply(packages, require, character.only = TRUE)
}

# function to obtain the selected variables and leverages before and after transformation

transregress <- function(X, y, family, mstop = NULL, nu = 0.1, cv, 
                         adaptive = T, ... ){
  
  # cv, matrix for cv folds
  require(mboost)
  X <- scale(X)
  Z <- transdata(X, byrow = F, eps = 0.01, scale = F)
  Z <- scale(Z)
  
  modX <- glmboost(X, y, 
                   family = family,
                   control = boost_control(mstop = mstop, 
                                           nu = nu, 
                                           risk = "inbag"), 
                   center = T)
  
  cv.x <- cvrisk(modX, folds = cv)
  cv.x <- mstop(cv.x)
  coefmodX <- unique(modX$xselect()[1:cv.x])
  
  hatX <- glmboost(X, y, 
                   family = family,
                   control = boost_control(mstop = cv.x, 
                                           nu = nu, 
                                           risk = "inbag"), 
                   center = T)
  
  hats <- hatX$hatvalues()[!sapply(hatX$hatvalues(), is.null)]
  hatX  <- Reduce(`+`, hats)
  
  modZ <- glmboost(Z, y, 
                   family = family,
                   control = boost_control(mstop = mstop, 
                                           nu = nu, 
                                           risk = "inbag"), 
                   center = T) 
  
  if (adaptive){
    cv.z <- cvrisk(modZ, folds = cv)
    cv.z <- mstop(cv.z)
  }else{
    cv.z = cv.x 
  }
  coefmodZ <- unique(modZ$xselect()[1:cv.z])
  
  hatZ <- glmboost(Z, y, 
                   family = family,
                   control = boost_control(mstop = cv.z, 
                                           nu = nu, 
                                           risk = "inbag"), 
                   center = T)
  
  hats <- hatZ$hatvalues()[!sapply(hatZ$hatvalues(), is.null)]
  hatZ <- Reduce(`+`, hats)
  
  return(list(X = coefmodX, 
              Z = coefmodZ, 
              cv.x = cv.x, 
              cv.z = cv.z,
              hatX = hatX,
              hatZ = hatZ))
}

#- distance method
#-------------------------------------------------------------------------------

moddist <- function(mod, std = F, L = 1) {
  
  #- std, whether p1 need to be standardized
  #- L, norm, 1, manhattan distance; 2, Euclidean distance
  n = length(mod)
  allvar <- unique(unlist(mod))
  f <- table(unlist(mod))[allvar]
  p1 <-  f / n
  p2 <- 1 - p1
  
  score <- mapply(function(x) {
    s <-  as.numeric(allvar %in% x) 
    sel <- p1 < 1
    s  <- s[sel]
    q1 <- p1[sel]
    q2 <- p2[sel]
    ifelse(std, sum(abs((s - q1)^L)) / (n * q1 * q2), 
           sum(abs((s - q1)^L)))
  }, mod)
  return(score)
  
}

#- hamming distance, not correct

ham <- function(si, sj){
  
  if(!is.list(sj)) sj = list(sj)
  ham <- mapply(function(sj){
    s1 <- length(setdiff(si, sj))
    s2 <- length(setdiff(sj, si))
    s1 + s2
  }, sj = sj)
  return(mean(ham))
}

#-------------------------------------------------------------------------------

hor.barplot <- function(x, digit = 2, decreasing = T, topn = 20, xlim = NULL, 
                        score = T, cex.names = par("cex.axis")){
  
  ord <- order(x, decreasing = decreasing)
  h = rev(x[ord[1:topn]])
  if (is.null(xlim)) xlim = c(0, 1.2 *max(h)) 
  p <- barplot(h, horiz = T, las = 2, xpd = F, xaxt = "n", 
               xlim = xlim, cex.names = cex.names)
  if (score == T) 
    text(x = h, y = p, labels = format(round(h, digit), nsamll = digit), 
         pos = 4, cex = 0.7)
}

#-------------------------------------------------------------------------------

lollipop <- function(x, decreasing = T, topn = 5, ylim = NULL, 
                     refline = 0, ylab = "", xlab = "Observations"){

  y = x
  x = seq_along(x)
  
  if (is.null(ylim))
    ylim = c(min(y), max(y))
  
  if (is.null(refline))
    refline = 0
  
  plot(x, y, xlab = "", ylab = "", cex = 0.8, ylim = ylim, pch = 19,
       yaxs = "i")
  abline(h = refline, lty = "dashed")
  
  for(i in x){
    lines(c(x[i], x[i]), c(0, y[i]), col = "grey")
  }
  
  ord <- order(abs(y), decreasing = T)[1:topn]
  text(x[ord], y[ord], pos = 3, cex = 0.9, labels = ord)
  mtext(xlab, side = 1, line = 2.2)

  mtext(ylab, side = 2, line = 2.2)
  abline(h = 0)
  
}

plotcluster <- function(boost, nlabs = 5, dist = "manhattan") {
  
  allvars <- unique(unlist(boost))
  boost <- mapply(function(x){
    as.numeric(allvars %in% x)
  }, x =  boost)
  d <- dist(t(boost), method = dist)
  
  fit    <- hclust(d, method = "average")
  n = length(fit$height) + 1
  points <- t(fit$merge)
  points <- rev(-points[points < 0]) 
  lab <- rep("", n)
  lab[points[1:nlabs]] <- points[1:nlabs]
  plot(fit, cex = 0.9,  main = "", labels = lab) 
}


#- variable selection ----------------------------------------------------------

plot_vs <- function(obj, method, name, xlim, plot = c("bar", "lol"),
                    center = F, refline = NULL){
  
  #- plot, compare to the original model
  
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

plot_m1vsm2 <- function(obj, plot.type = c("blandartman", "corrplot"), name){
  
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

# corplot to compare fixed and adaptive mstops
plot_fixvsada <- function(obj.fix, obj.ada, method = c("M1", "M2"), name){
  
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
  
  # plot
  corrplot(re.fix, re.ada, ylim = getlimits(re.ada),
           xlim = getlimits(re.fix),
           ci.line = F)
  mtext("Adaptive", side = 2, line = 2.2)
  mtext("Fixed", side = 1, line = 2.2)
  
  ord1 <- order(re.fix, decreasing = T)[1:3]
  ord2 <- order(re.ada, decreasing = T)[1:3]
  ord  <- unique(c(ord1, ord2))
  
  text(re.fix[ord], re.ada[ord], labels = ord, pos = 2, cex = 1)
  mtext(name, side =3, line = 0.2, adj = 0)
}

# boost path

vs_path <- function(X, y, family, obs.mstop = NULL, nu, 
                    nCores = NULL){
  
  X <- scale(X)
  require(mboost)
  require(doParallel)
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (is.null(nCores)) nCores <- detectCores() - 1
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  
  mstops = max(obs.mstop)
  reboost <- foreach(i = 1:n, 
                     .packages = c("mboost")) %dopar% {
                       
                       boost <- glmboost(X[-i,], y[-i], 
                                         family = family,
                                         control = boost_control(mstop = mstops, 
                                                                 nu = nu, 
                                                                 risk = "inbag"), 
                                         center = T)
                       
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
  
  out = list(scores = scores, mstops = obs.mstop)
  class(out) <- "vs_path"
  
  return(out)
}

plot.vs_path <- function(obj){
  
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


#- prediction ------------------------------------------------------------------

plot_predPath <- function(obj, top = 3, name = ""){
  
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
  
  # add title
  mtext(name, adj = 0, line = 0.2, side = 3)
}

plot_predScore <- function(obj, top = 3, ylim = NULL, name = ""){
  
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
  
  # add name
  mtext(name, adj = 0, line = 0.2, side = 3)
}

#-------------------------------------------------------------------------------

boostbeta <- function(X, y, family, mstop = NULL, obs.mstop = NULL, 
                      method = c("path", "cv"), 
                      measure = c("averaged", "original"),
                      ncores = NULL){
  
  if (isTRUE(all.equal(family, CoxPH())))
    y <- surv 
  if (isTRUE(all.equal(family, Binomial())))
    y <- as.factor(y)
  
  X <- scale(X)
  require(mboost)
  require(doParallel)
  
  n <- nrow(X)
  p <- ncol(X)

  if (is.null(ncores)) ncores <- detectCores() - 1
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  if (method == "path") mstops = rep(mstop, (n + 1))
  if (method == "cv")   mstops = c(obs.mstop, mstop)
  reboost <- foreach(i = 1:(n + 1), j = mstops, 
                     .packages = c("mboost")) %dopar% {
                       boost <- glmboost(X[-i,], y[-i], 
                                         family = family,
                                         control = boost_control(mstop = j, 
                                                                 nu = 0.1, 
                                                                 risk = "inbag"), 
                                         center = T)
                       
                       if (method == "path"){
                         tmp <- coef(boost, aggregate = "cumsum", off2int = F)
                         out <- Reduce(rbind, tmp)
                         rownames(out) <- names(tmp)
                       }
                       
                       if (method == "cv"){
                         out <- coef(boost)
                       }
                       
                       return(out)
                     }
  registerDoSEQ()
  stopCluster(cl)
  
  boost_drop1 <- reboost[-(n + 1)]
  boost_original <- reboost[[(n + 1)]]
  
  if (measure == "averaged") {
    
    difbeta <- function(beta){
      m <- unique(names(unlist(beta)))
      b <- rep(0, length(m))
      names(b) <- m
      matb <- mapply(function(x){
        b[names(x)] <- x
        return(b)
      }, x = beta)
      if (is.vector(matb)) matb <- matrix(matb, nrow = 1)
      meanb <- rowMeans(matb)
      apply(matb, 2, function(x) sum(abs(x - meanb)))
    }  
    
    if (method == "path"){ 
      out <- mapply(function(i){
        stopbeta <- lapply(boost_drop1, function(x) x[,i])
        difbeta(stopbeta)
      }, i = 1:mstop)
      class(stepbeta) <- append(class(out), "regpath")
    }
    
    if (method == "cv"){
      cvbeta <- difbeta(boost_drop1)
      out <- list(cvbeta = cvbeta, cv = cv)
      class(out) <- append(class(out), "cv")
    }
  }
  
  if (measure == "original") {
    
    difbeta <- function(beta, beta0){ 
      # beta0, beta in original model
      difbeta <- mapply(function(beta){
        m0 <- unique(c(names(beta), names(beta0)))
        b  <- rep(0, length(m0))
        names(b) <- m0
        b1 <- b2 <- b
        b1[names(beta)] <- beta
        b2[names(beta0)] <- beta0
        sum(abs(b1 - b2))
      }, beta = beta)
    }
    
    if (method == "path"){
      out <- mapply(function(i){
        stopbeta <- lapply(boost_drop1, function(x) x[,i])
        beta0    <- boost_original[,i]
        difbeta(stopbeta, beta0)
      }, i = 1:mstop)
      class(out) <- append(class(out), "regpath")
    }
    
    if (method == "cv"){
      cvbeta <- difbeta(boost_drop1, boost_original)
      out <- list(cvbeta = cvbeta, mstops = mstops)
      class(out) <- append(class(out), "cv")
    }
  }
  return(out)
}

plot.regpath <- function(obj, topn = 1, cases = NULL){
  
  n = ncol(obj)
  obj <- cbind(0, obj)
  plot.new()
  plot.window(ylim = c(0, max(obj)), xlim = c(0, n))
  if(n >= 10) 
    axis(1)
  if(n < 10) {
    axis(1, at = 0:n)
  }
  
  axis(2)
  box()
  apply(obj, 1, function(x){
    lines(0:n, x, col = "grey")
  })
  
  score <- apply(obj, 1, sum)
  if (is.null(cases))
    ind <- order(score, decreasing = T)[1:topn]
  if (!is.null(cases))
    ind <- cases
  
  apply(obj[ind,,drop = F], 1, function(x){
    lines(0:n, x, col = "black")
  })
  mtext("Boosting iterations", side = 1, line = 2.2)
  # mtext(bquote(Sigma ~ "|" ~ Delta ~ beta ~ "|"), side = 2, line = 2.2)
  mtext("DFBETA", side = 2, line = 2.2)
}

#-------------------------------------------------------------------------------

plot.boostpath <- function(mod, group = NULL, col = NULL, ylim = NULL,
                           abs = F){
  
  stepmod <- function(mod){
    steps <- mod$mstop()
    tmp  <- coef(mod, aggregate = "cumsum", off2int = F)
    beta <- Reduce(rbind, tmp)
    rownames(beta) <- names(tmp)
    return(beta)
  }
  
  beta <- stepmod(mod)
  beta <- cbind(0, beta)
  if (abs) beta <- abs(beta)
  
  plot.new()
  if (is.null(ylim)) ylim = c(min(beta), max(beta))
  plot.window(xlim = c(0, ncol(beta)), ylim = ylim)
  axis(1)
  axis(2)
  box()
  
  if (is.null(col)) col = "black"
  if (is.null(group)) group = list(rownames(beta))
  mapply(function(x, col){
    if (is.numeric(x)) x = paste0("V", x)
    D <- beta[x, ,drop = F]
    apply(D, 1, function(x){
      lines(0:(ncol(D) - 1), x, col = col)})
  }, x = group, col = col)
  abline(h = 0)
  mtext("Boosting iterations", side =  1, line = 2.2)
  if (abs) 
    mtext(bquote("|"~beta~"|"), side =  2, line = 2.2) else
      mtext(bquote(beta), side =  2, line = 2.2) 
}

#-------------------------------------------------------------------------------

loglik.coxph <- function(pred, surv.time, surv.event, verbose = FALSE) {	
  
  n <- length(pred)
  r <- rank(surv.time)
  ita <- pred
  epita <- exp(ita)
  d <- rep(0, n)
  dono <- rep(0, n)
  for(i in 1:n) {
    d[i] <- sum(surv.event[r == r[i]])
    dono[i] <- sum(epita[r >= r[i]])
  }

  lik <- (ita - log(dono)) * surv.event
  return(lik)
}


boost.pred.loglik <- function(X, y, family, 
                              mstop = NULL, 
                              obs.mstop = NULL, 
                              method = c("fixed", "cv"),
                              ncores = NULL,
                              path = F,
                              measure = c("loglik", "cook")){
  
  n <- nrow(X)
  p <- ncol(X)
  
  #- (n + 1) is the boosting without removing obs
  if (is.null(ncores)) ncores <- detectCores() - 1
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  if (method == "fixed") mstops = rep(mstop, (n + 1))
  if (method == "cv")    mstops = c(obs.mstop, mstop)
  if (path == T & measure == "cook") 
    stop("cook'D can only be calculated when path = F")
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
                       
                       m <- ifelse(path, 1, j)
                       out <- c()
                       
                       for (k in m:j){
                         if (isTRUE(all.equal(family, Binomial()))){
                           y0   <- as.numeric(y) - 1
                           pred <- predict(boost[k], newdata = X, type = "response")
                           loglik <- y0 * log(pred) + (1 - y0)*log(1 - pred)
                         }
                         
                         if (isTRUE(all.equal(family, CoxPH()))){
                           pred   <- predict(boost[k], newdata = X)
                           loglik <- loglik.coxph(pred, 
                                                  surv.time  = y[,1], 
                                                  surv.event = y[,2])	
                         }
                         if (measure == "cook") {
                           out <- loglik}else{
                             out <- c(out, sum(loglik))
                           }
                       }
                       return(out)
                     }
  registerDoSEQ()
  stopCluster(cl)
  
  if (!path){
    if(measure == "cook"){
      drop1 = reboost[,-(n + 1)]
      drop0 = reboost[,(n + 1)]
    }else{
      drop1 = reboost[-(n + 1)]
      drop0 = reboost[(n + 1)]
    }
  }
  
  if (path){
    drop1 = reboost[,-(n + 1)]
    drop0 = reboost[,(n + 1)]
  }
  
  out <- list(drop1 = drop1, drop0 = drop0, y = y, family = family, 
              mstop = mstops, method = method, measure = measure)
  return(out)
}

#- not necessary -- 

drop1.pred.loglik.plot <- function(obj, ylim = NULL, name = ""){
  
  drop1 <- obj$drop1
  drop0 <- obj$drop0
  y <- 2 * (drop0 - drop1)
  
  if (is.null(ylim)) 
    ylim <- c(min(y), max(y))
  
  n <- ncol(drop1)
  x <- 1:nrow(drop1)
  
  tmp.par <- par()
  plot.new()
  plot.window(xlim = c(1, nrow(drop1)), ylim = c(min(y), max(y)))
  axis(1)
  axis(2)
  box()
  
  for(i in 1:n){
    lines(x, y[,i], lwd = 0.8, col = "grey")
  }
  mtext(name, side = 3, line = 0.2, adj = 0)
  
  score <- colSums(y) 
  score <- (score - min(score)) / (max(score) - min(score))
  lines(x, y[,which.max(score)], lwd = 0.8, col = "black")
  
  mtext("number of iterations", side = 1, line = 2.2)
  mtext("Deviance changes", side = 2, line = 2.2)
  
  mar0 <- tmp.par$mar
  mar <- mar0 * c(0.75, 1, 1, 1)
  par(mar = mar)
  names(score) <- paste0("obs ", 1:n)
  hor.barplot(score, xlim = c(0, max(score) * 1.2))
  par(mar = mar0)
}

#-------------------------------------------------------------------------------

library(pheatmap)
library(grid)
library(gtable)

heatmap_motor <- function (matrix, border_color, cellwidth, cellheight, tree_col, 
    tree_row, treeheight_col, treeheight_row, filename, width, 
    height, breaks, color, legend, annotation_row, annotation_col, 
    annotation_colors, annotation_legend, annotation_names_row, 
    annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
    hjust_col, vjust_col, angle_col, fmat, fontsize_number, number_color, 
    gaps_col, gaps_row, labels_row, labels_col, ...) 
{
    lo = pheatmap:::lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), 
        ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, 
        treeheight_col = treeheight_col, treeheight_row = treeheight_row, 
        legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, 
        annotation_colors = annotation_colors, annotation_legend = annotation_legend, 
        annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, 
        main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
        fontsize_col = fontsize_col, angle_col = angle_col, gaps_row = gaps_row, 
        gaps_col = gaps_col, ...)
    res = lo$gt
    mindim = lo$mindim
    if (!is.na(filename)) {
        if (is.na(height)) {
            height = convertHeight(gtable_height(res), "inches", valueOnly = T)
        }
        if (is.na(width)) {
            width = convertWidth(gtable_width(res), "inches", valueOnly = T)
        }
        r = regexpr("\\.[a-zA-Z]*$", filename)
        if (r == -1) 
            stop("Improper filename")
        ending = substr(filename, r + 1, r + attr(r, "match.length"))
        f = switch(ending, pdf = function(x, ...) pdf(x, ...), 
            png = function(x, ...) png(x, units = "in", res = 300, 
                ...), jpeg = function(x, ...) jpeg(x, units = "in", 
                res = 300, ...), jpg = function(x, ...) jpeg(x, 
                units = "in", res = 300, ...), tiff = function(x, 
                ...) tiff(x, units = "in", res = 300, compression = "lzw", 
                ...), bmp = function(x, ...) bmp(x, units = "in", 
                res = 300, ...), stop("File type should be: pdf, png, bmp, jpg, tiff"))
        f(filename, height = height, width = width)
        gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, 
            border_color = border_color, tree_col = tree_col, 
            tree_row = tree_row, treeheight_col = treeheight_col, 
            treeheight_row = treeheight_row, breaks = breaks, 
            color = color, legend = legend, annotation_col = annotation_col, 
            annotation_row = annotation_row, annotation_colors = annotation_colors, 
            annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, 
            annotation_names_col = annotation_names_col, filename = NA, 
            main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
            fontsize_col = fontsize_col, hjust_col = hjust_col, 
            vjust_col = vjust_col, angle_col = angle_col, fmat = fmat, 
            fontsize_number = fontsize_number, number_color = number_color, 
            labels_row = labels_row, labels_col = labels_col, 
            gaps_col = gaps_col, gaps_row = gaps_row, ...)
        grid.draw(gt)
        dev.off()
        return(gt)
    }
    if (mindim < 3) 
        border_color = NA
    if (!is.na(main)) {
        elem = pheatmap:::draw_main(main, fontsize = 1.3 * fontsize, ...)
        res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main", 
            clip = "off")
    }
    if (!pheatmap:::is.na2(tree_col) & treeheight_col != 0) {
        elem = pheatmap:::draw_dendrogram(tree_col, gaps_col, horizontal = T)
        res = gtable_add_grob(res, elem, t = 2, l = 3, name = "col_tree")
    }
    if (!pheatmap:::is.na2(tree_row) & treeheight_row != 0) {
        elem = pheatmap:::draw_dendrogram(tree_row, gaps_row, horizontal = F)
        res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
    }
    elem = pheatmap:::draw_matrix(matrix, border_color, gaps_row, gaps_col, 
        fmat, fontsize_number, number_color)
    res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
        name = "matrix")
    if (length(labels_col) != 0) {
        pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, 
            hjust_col = hjust_col, vjust_col = vjust_col, angle_col = angle_col, 
            ...)
        elem = do.call(pheatmap:::draw_colnames, pars)
        res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", 
            name = "col_names")
    }
    if (length(labels_row) != 0) {
        pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, 
            ...)
        elem = do.call(pheatmap:::draw_rownames, pars)
        res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
            name = "row_names")
    }
    if (!pheatmap:::is.na2(annotation_col)) {
        converted_annotation = convert_annotations(annotation_col, 
            annotation_colors)
        elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
            gaps_col, fontsize, horizontal = T)
        res = gtable_add_grob(res, elem, t = 3, l = 3, clip = "off", 
            name = "col_annotation")
        if (annotation_names_col) {
            elem = pheatmap:::draw_annotation_names(annotation_col, fontsize, 
                horizontal = T)
            res = gtable_add_grob(res, elem, t = 3, l = 4, clip = "off", 
                name = "col_annotation_names")
        }
    }
    if (!pheatmap:::is.na2(annotation_row)) {
        converted_annotation = convert_annotations(annotation_row, 
            annotation_colors)
        elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
            gaps_row, fontsize, horizontal = F)
        res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", 
            name = "row_annotation")
        if (annotation_names_row) {
            elem = pheatmap:::draw_annotation_names(annotation_row, fontsize, 
                horizontal = F, hjust_col = hjust_col, vjust_col = vjust_col, 
                angle_col = angle_col)
            res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", 
                name = "row_annotation_names")
        }
    }
    annotation = c(annotation_col[length(annotation_col):1], 
        annotation_row[length(annotation_row):1])
    annotation = annotation[unlist(lapply(annotation, function(x) !pheatmap:::is.na2(x)))]
    if (length(annotation) > 0 & annotation_legend) {
        elem = pheatmap:::draw_annotation_legend(annotation, annotation_colors, 
            border_color, fontsize = fontsize, ...)
        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, 
            clip = "off", name = "annotation_legend")
    }
    if (!pheatmap:::is.na2(legend)) {
        elem = pheatmap:::draw_legend(color, breaks, legend, fontsize = fontsize, 
            ...)
        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 5, b = 5, 
            clip = "off", name = "legend")
    }
    return(res)
}

# Modified pheatmap:::lo    
lo <- function (rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, 
    treeheight_col, treeheight_row, legend, annotation_row, annotation_col, 
    annotation_colors, annotation_legend, annotation_names_row, 
    annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
    angle_col, gaps_row, gaps_col, ...) 
{
    if (!is.null(coln[1]) | (!pheatmap:::is.na2(annotation_row) & annotation_names_row)) {
        if (!is.null(coln[1])) {
            t = coln
        }
        else {
            t = ""
        }
        tw = strwidth(t, units = "in", cex = fontsize_col/fontsize)
        if (annotation_names_row) {
            t = c(t, colnames(annotation_row))
            tw = c(tw, strwidth(colnames(annotation_row), units = "in"))
        }
        longest_coln = which.max(tw)
        gp = list(fontsize = ifelse(longest_coln <= length(coln), 
            fontsize_col, fontsize), ...)
        coln_height = unit(1, "grobheight", textGrob(t[longest_coln], 
            rot = angle_col, gp = do.call(gpar, gp))) + unit(10, 
            "bigpts")
    }
    else {
        coln_height = unit(5, "bigpts")
    }
    if (!is.null(rown[1])) {
        t = rown
        tw = strwidth(t, units = "in", cex = fontsize_row/fontsize)
        if (annotation_names_col) {
            t = c(t, colnames(annotation_col))
            tw = c(tw, strwidth(colnames(annotation_col), units = "in"))
        }
        longest_rown = which.max(tw)
        gp = list(fontsize = ifelse(longest_rown <= length(rown), 
            fontsize_row, fontsize), ...)
        rown_width = unit(1, "grobwidth", textGrob(t[longest_rown], 
            rot = 0, gp = do.call(gpar, gp))) + unit(10, "bigpts")
    }
    else {
        rown_width = unit(5, "bigpts")
    }
    gp = list(fontsize = fontsize, ...)
    if (!pheatmap:::is.na2(legend)) {
        longest_break = which.max(nchar(names(legend)))
        longest_break = unit(1.1, "grobwidth", 
            textGrob(as.character(names(legend))[longest_break], 
            gp = do.call(gpar, gp)))
        title_length = unit(1.1, "grobwidth", textGrob("Scale", 
            gp = gpar(fontface = "bold", ...)))
        legend_width = unit(12, "bigpts") + longest_break * 1.2
        legend_width = max(title_length, legend_width)
    }
    else {
        legend_width = unit(0, "bigpts")
    }
    if (is.na(main)) {
        main_height = unit(0, "npc")
    }
    else {
        main_height = unit(1.5, "grobheight", textGrob(main, 
            gp = gpar(fontsize = 1.3 * fontsize, ...)))
    }
    textheight = unit(fontsize, "bigpts")
    if (!pheatmap:::is.na2(annotation_col)) {
        annot_col_height = ncol(annotation_col) * (textheight + 
            unit(2, "bigpts")) + unit(2, "bigpts")
        t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col))
        annot_col_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
            gp = gpar(...))) + unit(12, "bigpts")
        if (!annotation_legend) {
            annot_col_legend_width = unit(0, "npc")
        }
    }
    else {
        annot_col_height = unit(0, "bigpts")
        annot_col_legend_width = unit(0, "bigpts")
    }
    if (!pheatmap:::is.na2(annotation_row)) {
        annot_row_width = ncol(annotation_row) * (textheight + 
            unit(2, "bigpts")) + unit(2, "bigpts")
        t = c(as.vector(as.matrix(annotation_row)), colnames(annotation_row))
        annot_row_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
            gp = gpar(...))) + unit(12, "bigpts")
        if (!annotation_legend) {
            annot_row_legend_width = unit(0, "npc")
        }
    }
    else {
        annot_row_width = unit(0, "bigpts")
        annot_row_legend_width = unit(0, "bigpts")
    }
    annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)
    treeheight_col = unit(treeheight_col, "bigpts") + unit(5, 
        "bigpts")
    treeheight_row = unit(treeheight_row, "bigpts") + unit(5, 
        "bigpts")
    if (is.na(cellwidth)) {
        mat_width = unit(1, "npc") - rown_width - legend_width - 
            treeheight_row - annot_row_width - annot_legend_width
    }
    else {
        mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * 
            unit(4, "bigpts")
    }
    if (is.na(cellheight)) {
        mat_height = unit(1, "npc") - main_height - coln_height - 
            treeheight_col - annot_col_height
    }
    else {
        mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * 
            unit(4, "bigpts")
    }
    gt = gtable(widths = unit.c(treeheight_row, rown_width,  
        mat_width, treeheight_row, legend_width, annot_legend_width), 
        heights = unit.c(main_height, treeheight_col, annot_col_height, 
            mat_height, coln_height), vp = viewport(gp = do.call(gpar, 
            gp)))
    cw = convertWidth(mat_width - (length(gaps_col) * unit(4, 
        "bigpts")), "bigpts", valueOnly = T)/ncol
    ch = convertHeight(mat_height - (length(gaps_row) * unit(4, 
        "bigpts")), "bigpts", valueOnly = T)/nrow
    mindim = min(cw, ch)
    res = list(gt = gt, mindim = mindim)
    return(res)
}

# Modified pheatmap:::draw_rownames      
draw_rownames <- function (rown, gaps, ...) 
{
    coord = pheatmap:::find_coordinates(length(rown), gaps)
    y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
    res = textGrob(rown, x = unit(-3, "bigpts"), y = y, vjust = 0.5, 
        hjust = 1, gp = gpar(...))
    return(res)
}

assignInNamespace(x = "draw_rownames", value = draw_rownames, ns="pheatmap")
assignInNamespace(x = "lo", value = lo, ns = "pheatmap")
assignInNamespace(x = "heatmap_motor", value = heatmap_motor, ns = "pheatmap") 

#-------------------------------------------------------------------------------

#' BlandArtman plot
#'
#' This funcition return a BlandArtman Plot
#'
#' @inheritParams blandartman
#' @return A plot
#' @param x1,x2 two vectos for comparison
#' @param xlab, xlab
#' @param ylab, ylab
#' @param ci, ci for the mean
#' @param ci.limit, ci for the limits
#' @param spline, logical for spline
#' @param ci.spline, confidence interval for spline
#' @param xlim, xlim
#' @param ylim ylim
#' @seealso correlation plot
#'
#' @export
#' @examples

blandartman <- function(x1, x2, xlab = "", ylab = "", ci = TRUE, ci.limit = TRUE,
                        spline = TRUE, ci.spline = TRUE, xlim = c(), ylim = c(),
                        alpha = 0.95, beta = 0.95, withinlimits = T, LOA = T){
  require(dplyr)
  y <- x1 - x2
  sdy <- sd(y)
  stdy <- sdy / sqrt(length(y))
  x <- (x1 + x2)/2
  # proportion beyond the limits
  
  z.alpha <- qnorm(alpha + (1 - alpha) / 2)
  z.beta  <- qnorm(beta + (1 - beta) / 2)
  
  out <- length(which(abs(y - mean(y)) > z.alpha * sdy))
  out <- round(1 - out / length(y), 3) * 100
  
  if(is.null(xlim)){xlim = c(min(x), max(x))}
  if(is.null(ylim)){ylim = c(-max(abs(y)), max(abs(y)))}
  
  #trim observations in x axis
  d <- data.frame(x,y)
  din <- d %>% filter(x > xlim[1], x < xlim[2], y < ylim[2], y > ylim[1])
  dout <- d %>% filter(x < xlim[1] | x > xlim[2] | y < ylim[1] | y > ylim[2])
  ifextreme<-length(dout)>0
  if(ifextreme){
    dout[dout$x < xlim[1], "x"] <- xlim[1]
    dout[dout$x > xlim[2], "x"] <- xlim[2]
    dout[dout$y < ylim[1], "y"] <- ylim[1]
    dout[dout$y > ylim[2], "y"] <- ylim[2]
  }
  #
  plot.new()
  plot.window(xlim,ylim)
  axis(1)
  axis(2)
  box()
  points(din$x, din$y, pch = 19, col = "grey45", cex = 0.8)
  mtext(xlab, side = 1, line = 2.2)
  mtext(ylab, side = 2, line = 2.2)
  abline(h = mean(y), col = "black", lwd = 1) # bias
  
  if(ifextreme){
    points(dout$x,dout$y, pch = 1)
  }
  
  if(ci == TRUE){
    abline(h = mean(y) - (z.beta * stdy), lty = "dashed", lwd = 1.2)
    abline(h = mean(y) + (z.beta * stdy), lty = "dashed", lwd = 1.2)
  }
  
  if (LOA) {
    abline(h = mean(y) - (z.alpha * sdy), lty = "dashed", lwd = 1)
    abline(h = mean(y) + (z.alpha * sdy), lty = "dashed", lwd = 1)
  }
  # Limits of agreement with 95% CI
  
  n1 = length(y)
  t = qt(p = 0.975, df = (n1 - 1))
  sdy_limit <- sqrt((1 / n1 + (1.96^2 / (2 * (n1-1)))) * sdy^2)
  limit1 = mean(y) + (1.96 * sdy)
  limit2 = mean(y) - (1.96 * sdy)
  
  if(ci.limit == TRUE){
    rect(-abs(xlim[1]) * 1.1,  limit1 - (t * sdy_limit), abs(xlim[2]) * 1.1, 
         limit1 + (t * sdy_limit), col ="NA", border = TRUE, lty = "dashed")
    rect(-abs(xlim[1]) * 1.1,  limit2 - (t * sdy_limit), abs(xlim[2]) * 1.1, 
         limit2 + (t * sdy_limit), col ="NA", border = TRUE, lty = "dashed")
  }
  
  # abline(h = mean(y) + (1.96 * sdy), col = "black", lwd = 1, lty = "dashed")
  # abline(h = mean(y) - (1.96 * sdy), col = "black", lwd = 1, lty = "dashed")
  
  ord <- order(x, decreasing = FALSE)
  x <- x[ord]
  y <- y[ord]
  alpha = 0.05
  if(spline == TRUE){
    re <- predict(loess(y ~ x, span = 0.40), se = TRUE)
    lines(x,re$fit)
    if(ci.spline == TRUE){
      upper <- re$fit + qnorm(1 - alpha / 2) * re$se.fit
      lower <- re$fit - qnorm(1 - alpha / 2) * re$se.fit
      lines(x, upper, lty = "dashed")
      lines(x, lower, lty = "dashed")
    }
  }
  
  if (withinlimits){
    ltext<-paste0("Within limits: ", out, "%")
    legend("top", legend = ltext, bty = "n", cex = 1.1)
  }
}

# help function for blandartman
blandartmanAddmarks <- function(x1, x2, topn = 5, 
                                alpha = 0.95, pos= 4, cex = 1){
  
  y <- x1 - x2
  sdy <- sd(y)
  stdy <- sdy / sqrt(length(y))
  x    <- (x1 + x2) / 2
  
  # proportion beyond the limits
  z.alpha <- qnorm(alpha + (1 - alpha) / 2)
  ind.out <- which(abs(y - mean(y)) > z.alpha * sdy)
  
  tops <- order(abs(y), decreasing = T)[1:topn]
  out  <- intersect(ind.out, tops)
  
  if (length(out) == 0)
    print("No case outside the LOA")
  
  if (length(out) > 0)
    text(x[out], y[out], labels = out, cex = cex, pos = pos)
}


#' Correlation plot
#'
#' This funcition return a correlation plot
#'
#' @inheritParams corrplot
#' @return A plot
#' @param x1,x2 two vectos for comparison
#' @param xlab, xlab
#' @param ylab, ylab

#' @param xlim, xlim
#' @param ylim ylim
#' @seealso correlation plot
#'
#' @export
#' @examples
#' 
corrplot <- function(x, y, 
                     xlab = "", ylab = "", 
                     xlim = c(), ylim = c(),
                     p.col = "grey60",
                     line = T, 
                     ci.line = T){
  
  require(dplyr)
  ord <- order(x, decreasing = FALSE)
  x <- x[ord]
  y <- y[ord]
  
  if(is.null(xlim) & is.null(ylim)){
    max1 <- max(c(x, y))
    min1 <- min(c(x, y))
    ylim = c(min1, max1)
    xlim = c(min1, max1)
  }
  if(is.null(xlim)) {xlim = ylim}
  if(is.null(ylim)) {ylim = xlim}
  
  d <- data.frame(x, y)
  din <-d %>% filter(x < xlim[2], x > xlim[1], y < ylim[2], y > ylim[1])
  dout <- d %>% filter(x < xlim[1] | x > xlim[2] | y < ylim[1] | y > ylim[2])
  mod  <- lm(y ~ x)
  beta <- coef(mod)[2]
  addextreme = nrow(dout) > 0
  if(addextreme){
    dout <- dout %>% mutate(y1 = ifelse(x > x[1], (xlim[2] - x) * beta + y, 
                                        (xlim[1] - x) * beta + y),
                            x1 = ifelse(x > x[1], xlim[2], xlim[1]),
                            x2 = ifelse(y > y[1], (ylim[2]-y) * (1 / beta) + x, 
                                        (ylim[1] - y) * (1 / beta) + x),
                            y2 = ifelse(y > y[1], ylim[2], ylim[1]))
    # 4 possibilities
    dout <- dout %>% mutate(c1 = ifelse(abs(y1) + abs(x1) - (x[2] - x[1]) < 0, 
                                        sign(x1 + y1) * 1, sign(x2 + y2) * 2))
    xyfun<-function(x, dout){
      switch(x,"1"= {temp <- dout[dout$c1 == 1, "y1"]; 
      return(list(xy = c(xlim[2], mean(temp)), n = length(temp)))},
      "-1" = {temp <- dout[dout$c1 == -1, "y1"]; 
      return(list(xy = c(xlim[1], mean(temp)), n = length(temp)))},
      "2" = {temp <- dout[dout$c1 == 2, "x2"]; 
      return(list(xy = c(mean(temp), xlim[2]), n = length(temp)))},
      "-2" = {temp <- dout[dout$c1 == -2, "x2"]; 
      return(list(xy = c(mean(temp), xlim[1]), n = length(temp)))})
    }
    xy <- lapply(as.character(unique(dout$c1)), function(x) xyfun(x, dout))
  }
  #plot
  plot.new()
  plot.window(xlim, ylim)
  axis(1)
  axis(2)
  mtext(xlab, side = 1,line = 2.1)
  mtext(ylab, side = 2,line = 2.1)
  box()
  points(y ~ x, data = din, pch = 19, col = p.col)
  k <- lm(y~x)
  d <- vcov(k)
  ci = sqrt(d[1,1] + d[2,2] * x^2 + 2*d[1,2]*x)
  #select x, draw the points and lines in limits
  xin <- x %in% din$x
  if (line)
    lines(k$fitted.values~x, type="l", xlab="", ylab="")
  if (ci.line) {
    upper <- k$fitted.values + 1.96 * ci
    lower <- k$fitted.values - 1.96 * ci
    lines(x,upper, lty = "dashed")
    lines(x,lower, lty = "dashed")
  }
  if(addextreme){
    for(i in 1:length(xy)){
      temp <- xy[[i]]
      points(temp$xy[1], temp$xy[2], pch = as.character(temp$n), cex = 1.2)
    }
  }
}

getlimits <- function(x, p = 0.05){
  
  lower <- min(x) - p * (max(x) - min(x))
  upper <- max(x) + p * (max(x) - min(x))
  return(c(lower, upper))
}

#-------------------------------------------------------------------------------
#- CVrisk for boosting
#-------------------------------------------------------------------------------


#- a simple modification from cvrisk.boost
#- achieve the prediction of the deleted case

cvrisk1 <- function (object, folds = cv(model.weights(object)),
                     grid = 0:mstop(object), papply = mclapply,
                     fun = NULL, mc.preschedule = FALSE,
                     ...) {
  
  papply  <- match.fun(papply)
  weights <- model.weights(object)
  
  if (any(weights == 0))
    warning("zero weights")
  
  if (is.null(folds)) {
    
    folds <- rmultinom(25, length(weights), weights/sum(weights))
    attr(folds, "type") <- "25-fold bootstrap"
    
  } else {
    stopifnot(is.matrix(folds) && nrow(folds) == length(weights))
  }

  fitfct  <- object$update
  oobrisk <- matrix(0, nrow = ncol(folds), ncol = length(grid))
  
  if (!is.null(fun))
    stopifnot(is.function(fun))
  fam_name <- object$family@name
  call <- deparse(object$call)
  
  if (is.null(fun)) {
    
    dummyfct <- function(weights, oobweights) {
      
      mod <- fitfct(weights = weights, oobweights = oobweights)
      mstop(mod) <- max(grid)
      ## return all risk values in grid (+ 1 as 0 is included)
      risk(mod)[grid + 1]
    }
    if (fam_name == "Cox Partial Likelihood" && all(colSums(folds == 0) == 1))
      stop("Leave-one-out cross-validation cannot be used with ", 
           sQuote("family = CoxPH()"))
    
  } else { ## !is.null(fun)
    
    dummyfct <- function(weights, oobweights) {
      
      mod <- fitfct(weights = weights, oobweights = oobweights)
      mod[max(grid)]
      ## make sure dispatch works correctly
      class(mod) <- class(object)
      fun(mod)
    }
  }

  #- Modification begin
  ## use case weights as out-of-bag weights (but set inbag to 0)
  # OOBweights <- matrix(rep(rep(1, length(weights)), ncol(folds)), 
  #                      ncol = ncol(folds))
  
  OOBweights <- abs(folds - 1) 
  
  # #---------
  # 
  # k = 8
  # w <- rep(1, n)
  # if (k <= n) w[k] <- 0
  # 
  # object <- glmboost(X, y, 
  #                    family = family,
  #                    weights = w,
  #                    control = boost_control(mstop = no.mstop, 
  #                                            nu = 0.1, 
  #                                            risk = "inbag"), 
  #                    center = F)
  # 
  # fitfct  <- object$update
  # weights <- model.weights(object)
  # 
  # #-
  # i = 1
  # s <- fitfct(weights = folds[, i] * weights, 
  #             oobweights = OOBweights[, i])
  # 
  # out <- OOBweights[, i]
  # mstop(s)
  # pout <- s$predict(newdata = X[out == 1,])
  # yout <- as.numeric(y[out == 1]) - 1
  # 
  # Binomial()@risk(pout, yout)
  # Binomial()@risk(s$fitted(), as.numeric(y) - 1)
  # 
  # 
  # #-
  # 
  # re <- c()
  # for(i in 1:10){
  #   tmp <- dummyfct(weights = folds[, i] * weights, oobweights = OOBweights[, i])
  #   re <- cbind(re, tmp)
  # }
  # 
  # sum(re8[76 + 1,]) - sum(re0[68 + 1,])
  # sum(re1[66 + 1,]) - sum(re0[68 + 1,])
  
  #-------
  
  if (identical(papply, mclapply)) {
    
    oobrisk <- papply(1:ncol(folds),
                      function(i) try(dummyfct(weights = folds[, i] * weights,
                                               oobweights = OOBweights[, i]),
                                      silent = TRUE, ...),
                      mc.preschedule = mc.preschedule)
  } else {
    
    dummyfct(weights = folds[, i] * weights,
             oobweights = OOBweights[, i])
    
    oobrisk <- papply(1:ncol(folds),
                      function(i) try(dummyfct(weights = folds[, i] * weights,
                                               oobweights = OOBweights[, i]),
                                      silent = TRUE), ...)
  }
  
  ## if any errors occured remove results and issue a warning
  if (any(idx <- sapply(oobrisk, is.character))) {
    
    warning(sum(idx), " fold(s) encountered an error. ",
            "Results are based on ", ncol(folds) - sum(idx),
            " folds only.\n",
            "Original error message(s):\n",
            sapply(oobrisk[idx], function(x) x))
    oobrisk[idx] <- NA
  }
  
  if (!is.null(fun))
    return(oobrisk)
  
  oobrisk <- t(as.data.frame(oobrisk))
  
  #- modification
  # oobrisk <- oobrisk / colSums(OOBweights)
  # oobrisk <- oobrisk
  # -end
  
  colnames(oobrisk) <- grid
  rownames(oobrisk) <- 1:nrow(oobrisk)
  attr(oobrisk, "risk") <- fam_name
  attr(oobrisk, "call") <- call
  attr(oobrisk, "mstop") <- grid
  attr(oobrisk, "type")  <- ifelse(!is.null(attr(folds, "type")),
                                   attr(folds, "type"), "user-defined")
  class(oobrisk) <- "cvrisk"
  oobrisk
}

print.cvrisk <- function(x, ...) {
  cat("\n\t Cross-validated", attr(x, "risk"), "\n\t",
      attr(x, "call"), "\n\n")
  print(colMeans(x, na.rm = TRUE))
  cat("\n\t Optimal number of boosting iterations:", mstop(x), "\n")
  return(invisible(x))
}

plot.cvrisk <- function(x, xlab = "Number of boosting iterations",
                        ylab = attr(x, "risk"),
                        ylim = range(x), main = attr(x, "type"), ...) {

  ## force evaluation of attributes (to not evaluate them in the wrong situation)
  force(ylab)
  force(main)

  mstops <- attr(x, "mstop")

  x <- x[, apply(x, 2, function(y) all(!is.na(y))), drop = FALSE]
  cm <- colMeans(x)
  plot(1:ncol(x), cm, ylab = ylab, ylim = ylim,
       type = "n", lwd = 2, xlab = xlab,
       main = main, axes = FALSE, ...)
  out <- apply(x, 1, function(y) lines(1:ncol(x),y, col = "lightgrey"))
  rm(out)
  ms <- which.min(cm)
  lines(c(ms, ms), c(min(c(0, ylim[1] * ifelse(ylim[1] < 0, 2, 0.5))), cm[ms]),
        lty = 2)
  lines(1:ncol(x), cm, type = "l")
  axis(1, at = 1:ncol(x), labels = mstops)
  axis(2)
  box()
}

mstop.cvrisk <- function(object, ...)
  attr(object, "mstop")[which.min(colSums(object, na.rm = TRUE))]

cv <- function(weights, type = c("bootstrap", "kfold", "subsampling"),
               B = ifelse(type == "kfold", 10, 25),
               prob = 0.5, strata = NULL) {

  type <- match.arg(type)
  n <- length(weights)

  if (is.null(strata)) strata <- gl(1, n)
  if (!is.factor(strata)) stop(sQuote("strata"), " must be a factor")
  folds <- matrix(0, nrow = n, ncol = B)

  ### <FIXME> handling of weights needs careful documentation </FIXME>
  for (s in levels(strata)) {
    indx <- which(strata == s)
    folds[indx,] <- switch(type,
                           "bootstrap" = cvboot(length(indx), B = B, weights[indx]),
                           "kfold" = cvkfold(length(indx), k = B) * weights[indx],
                           "subsampling" = cvsub(length(indx), prob = prob, B = B) * weights[indx])
  }
  attr(folds, "type") <- paste(B, "-fold ", type, sep = "")
  return(folds)
}


cvboot <- function(n, B, weights)
  rmultinom(B, n, weights / sum(weights))

cvkfold <- function(n, k) {
  #if (k > n / 2) stop("k > n/2")
  fl <- floor(n/k)
  folds <- c(rep(c(rep(0, fl), rep(1, n)), k - 1),
             rep(0, n * k - (k - 1) * (fl + n)))
  matrix(folds, nrow = n)[sample(1:n),, drop = FALSE]
}

cvsub <- function(n, prob, B) {
  k <- floor(n * prob)
  indx <- rep(c(0, 1), c(n - k, k))
  replicate(B, sample(indx))[sample(1:n),, drop = FALSE]
}


# #- end -----------------------------------------------------------------------




