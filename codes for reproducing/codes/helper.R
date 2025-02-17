
#-------------------------------------------------------------------------------

load_packages <- function(packages){
  
  # Check if packages are installed
  not_installed <- packages[!(packages %in% installed.packages()[,"Package"])]
  
  # Install packages if they are not already installed
  if (length(not_installed) > 0) {
    install.packages(not_installed)
  }
  
  # Load packages
  suppressPackageStartupMessages({
    lapply(packages, require, character.only = TRUE)
  })
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
    
    # add sel when std = TRUE
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

lollipop <- function(x, 
                     decreasing = T, 
                     topn = 5, 
                     ylim = NULL, 
                     refline = 0, 
                     ylab = "", 
                     xlab = "Observations"){

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
  
  #- points
  points(x, y, cex = 0.8, ylim = ylim, pch = 19)
  
  ord <- order(abs(y), decreasing = T)[1:topn]
  text(x[ord], y[ord], pos = 3, cex = 0.9, labels = ord)
  mtext(xlab, side = 1, line = 2.2)

  mtext(ylab, side = 2, line = 2.2)
  abline(h = 0)
  
}

#- variable selection ----------------------------------------------------------

plot_vs <- function(obj, 
                    method, 
                    name, 
                    xlim, 
                    plot = c("bar", "lol"),
                    center = F, 
                    refline = NULL){
  
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

plot_m1vsm2 <- function(obj, 
                        plot.type = c("blandartman", "corrplot"), 
                        name = ""){
  
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
  
  # use same limits both in x and y axes 
  lim <- c(0, max(c(getlimits(re.dif2)), getlimits(re.dif1)))
  
  if (plot.type == "corrplot"){
    
    corrplot(re.dif1, re.dif2, ylim = lim, xlim = lim, ci.line = F, line = F)
    mtext("M1", side = 1, line = 2.2)
    mtext("M2", side = 2, line = 2.2)
    
    ord1 <- order(re.dif1, decreasing = T)[1:5]
    ord2 <- order(re.dif2, decreasing = T)[1:5]
    ord  <- unique(c(ord1, ord2))
    
    text(re.dif1[ord], re.dif2[ord], labels = ord, pos = 2, cex = 0.9)
  }
  mtext(name, side = 3, line = 0.2, adj = 0)
  
  # add diagonal line 
  lines(c(-100, 100), c(-100, 100), col = "gray")
  
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
  corrplot(re.fix, re.ada, ylim = getlimits(c(re.fix, re.ada)),
           xlim = getlimits(c(re.fix, re.ada)), line = FALSE, 
           ci.line = F)
  lines(c(-100, 100), c(-100, 100))
  mtext("Adaptive", side = 2, line = 2.2)
  mtext("Fixed", side = 1, line = 2.2)
  
  ord1 <- order(re.fix, decreasing = T)[1:3]
  ord2 <- order(re.ada, decreasing = T)[1:3]
  ord  <- unique(c(ord1, ord2))
  
  text(re.fix[ord], re.ada[ord], labels = ord, pos = 2, cex = 1)
  mtext(name, side =3, line = 0.2, adj = 0)
}

# boost path

vs_path <- function(X, y, family, obs.mstop = NULL, nu, method = "M2", nCores = NULL){
  
  X <- scale(X)
  require(mboost)
  require(doParallel)
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (is.null(nCores)) nCores <- detectCores() - 1
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  
  mstops = max(obs.mstop)
  reboost <- foreach(i = 1:(n + 1), 
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
    if (method == "M2")
      scores <- rbind(scores, moddist(mod[-length(mod)]))
    if (method == "M1")
      si <- mod[[length(mod)]]
      scores <- rbind(scores, mapply(ham, sj = mod[-length(mod)], MoreArgs = list(si = si)))
  }
  
  out = list(scores = scores, mstops = obs.mstop)
  class(out) <- "vs_path"
  
  return(out)
}

plot.vs_path <- function(obj, ylab = "M2"){
  
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
  mtext(top2, side = 4, at = at, line = 0, las = 2, cex = 0.7)
  axis(side = 3, at = mstops[top2], tick = T, cex = 0.6,  tcl = 0.3, labels = F)
  mtext(top2, side = 3, at = mstops[top2], line = 0, cex = 0.7)
  
  #- add labels
  mtext("Boosting iterations", line = 2.2, side = 1)
  mtext(ylab, line = 2.2, side = 2)
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
    max1 <- max(c(x, y)) * 1.1
    min1 <- min(c(x, y)) * 0.9
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

cvrisk1 <- function (object, 
                     folds = cv(model.weights(object)),
                     grid = 0:mstop(object), 
                     papply = mclapply,
                     fun = NULL, 
                     mc.preschedule = FALSE,
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

cv <- function(weights, 
               type = c("bootstrap", "kfold", "subsampling"),
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

## plot_predPoint, prediction 

plot_predPoint <- function(obj, 
                           mstops = NULL, 
                           top = 5, 
                           ylim = NULL, 
                           name = "",
                           ref = c("mean", "orig")){

  ## mstops should be a vector of stops
  n <- length(obj)
  
  ## cvm changes by -obs
  cvm <- c()
  for(i in seq_along(mstops)) {
    cvm <- c(cvm, colMeans(obj[[i]])[mstops[i] + 1])
  }
  
  # reference
  if (ref == "mean")
    ref <- mean(cvm[-n])
  if (ref == "orig")
    ref <- cvm[n]
  
  scores <- abs(cvm - ref)[-n]
  scores <- scores / sd(scores)
  
  ## lollipop
  lollipop(scores, topn = top, ylim = ylim, ylab = bquote(Delta ~ "CVM"), 
           xlab = "Observations")
  mtext(name, side = 3, adj = 0, line = 0.2)
  
}


plot_predCorr <- function(obj, 
                          mstops, 
                          ref = c("mean", "orig"),
                          top = 3, 
                          ylim = NULL, 
                          xlim = NULL, 
                          name = ""){
  
  ## mstops should be a vector of stops
  n <- length(obj)
  
  # cvm changes by -obs
  cvm_ada <- c()
  cvm_fix <- c()
  
  for(i in seq_along(mstops)) {
    cvm_ada <- c(cvm_ada, colMeans(obj[[i]])[mstops[i] + 1])
    cvm_fix <- c(cvm_fix, colMeans(obj[[i]])[mstops[n] + 1])
  }
  
  # reference
  if (ref == "mean") {
    ref_ada <- mean(cvm_ada[-n])
    ref_fix <- mean(cvm_fix[-n])
  }
  
  if (ref == "orig"){
    ref_ada <- cvm_ada[n]
    ref_fix <- cvm_fix[n]
  }
  
  scores_ada <- abs(cvm_ada - ref_ada)[-n]
  scores_ada <- scores_ada / sd(scores_ada)
  scores_fix <- abs(cvm_fix - ref_fix)[-n]
  scores_fix <- scores_fix / sd(scores_fix)
  # lollipop
  corrplot(x = scores_ada, 
           y = scores_fix, 
           xlab = "Adaptive", ylab = "Fixed", 
           xlim = xlim, 
           ylim = ylim,
           p.col = "grey60",
           line = F, 
           ci.line = F)
  
  lines(c(-100, 100), c(-100, 100))
  mtext(name, side = 3, adj = 0, line = 0.2)
  
  ord1 <- order(scores_ada, decreasing = T)[seq_len(top)]
  ord2 <- order(scores_fix, decreasing = T)[seq_len(top)]
  
  ord <- unique(c(ord1, ord2))
  text(x = scores_ada[ord], 
       y = scores_fix[ord],
       labels = ord,
       cex = 0.9,
       pos = 4)
  
}



# #- end -----------------------------------------------------------------------




