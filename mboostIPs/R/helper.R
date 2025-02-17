
moddist <- function(mod, std = F, L = 1) {

  #- std, whether p1 need to be standardized
  #- L, norm, 1, manhattan distance; 2, Euclidean distance
  
  n = length(mod)
  allvar <- as.character(unique(unlist(mod)))
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


ham <- function(si, sj){

  if(!is.list(sj)) sj = list(sj)
  ham <- mapply(function(sj){
    s1 <- length(setdiff(si, sj))
    s2 <- length(setdiff(sj, si))
    s1 + s2
  }, sj = sj)
  return(mean(ham))
}


# lollipop function, used in plot_Scores

lollipop <- function(x, 
                     decreasing = T, 
                     topN = 5, 
                     ylim = NULL, 
                     refline = 0, 
                     ylab = "", 
                     xlab = "Observations",
                     ...){
  
  # x, influential scores
  y = x
  x = seq_along(x)
  
  if (is.null(ylim))
    ylim = c(min(y), max(y) * 1.1)
  
  if (is.null(refline))
    refline = 0
  
  plot(x, y, xlab = "", ylab = "", cex = 0.8, ylim = ylim, pch = 19,
       yaxs = "i", ... )
  abline(h = refline, lty = "dashed")
  
  for(i in x){
    lines(c(x[i], x[i]), c(0, y[i]), col = "grey")
  }
  
  #- points
  points(x, y, cex = 0.8, ylim = ylim, pch = 19, ...)
  
  ord <- order(abs(y), decreasing = T)[1:topN]
  text(x[ord], y[ord], pos = 3, cex = 0.9, labels = ord)
  mtext(xlab, side = 1, line = 2.2)
  
  mtext(ylab, side = 2, line = 2.2)
  abline(h = 0)
  
}

## cvrisk modified from mboost

cvrisk_mboost <- function(object, folds = cv(model.weights(object)),
                          grid = 0:mstop(object), papply = mclapply,
                          fun = NULL, mc.preschedule = FALSE,
                          ...) {
  
  papply <- match.fun(papply)
  weights <- stats::model.weights(object)
  if (any(weights == 0))
    warning("zero weights")
  if (is.null(folds)) {
    folds <- rmultinom(25, length(weights), weights/sum(weights))
    attr(folds, "type") <- "25-fold bootstrap"
  } else {
    stopifnot(is.matrix(folds) && nrow(folds) == length(weights))
  }
  fitfct <- object$update
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
  
  ## use case weights as out-of-bag weights (but set inbag to 0)
  OOBweights <- matrix(rep(weights, ncol(folds)), ncol = ncol(folds))
  OOBweights[folds > 0] <- 0
  
  if (identical(papply, mclapply)) {
    oobrisk <- papply(1:ncol(folds),
                      ## here changes,add  * weights
                      function(i) try(dummyfct(weights = folds[, i] * weights,
                                               oobweights = OOBweights[, i]),
                                      silent = TRUE),
                      mc.preschedule = mc.preschedule
                      # ,...
                      )
  } else {
    oobrisk <- papply(1:ncol(folds),
                      ## here changes,add  * weights
                      function(i) try(dummyfct(weights = folds[, i] * weights,
                                               oobweights = OOBweights[, i]),
                                      silent = TRUE),
                      ...)
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
  
  ## here changes
  oobrisk <- oobrisk / colSums(OOBweights)
  colnames(oobrisk) <- grid
  rownames(oobrisk) <- 1:nrow(oobrisk)
  attr(oobrisk, "risk") <- fam_name
  attr(oobrisk, "call") <- call
  attr(oobrisk, "mstop") <- grid
  attr(oobrisk, "type") <- ifelse(!is.null(attr(folds, "type")),
                                  attr(folds, "type"), "user-defined")
  class(oobrisk) <- "cvrisk"
  oobrisk
}



























