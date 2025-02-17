

#' glmboostDrop1 Function
#'
#' @description
#' The `glmboostDrop1` function performs a drop-one analysis for models fitted 
#' with gradient boosting, specifically targeting variable selection and 
#' prediction improvements. It can operate in parallel to speed up computations.
#'
#' @param obj An object of class `mboost` as returned by the `mboost` package 
#' functions. This is the model on which the drop-one analysis will be performed.
#' @param nCores Integer specifying the number of cores to use for parallel execution.
#' If `NULL`, the function will use one less than the total number of available cores.
#' @param fixMstop Optional; an integer specifying using a fixed mstop 
#' for the gradient boosting algorithm. If `NULL`, `mstop` is determined automatically.
#' @param aims A character vector indicating the aim of the analysis: 
#' "variableSelection" to detect IPs for variable seleciton, "prediction" to
#' detect IPs for prediction, or both.
#' @param ref A character string specifying the reference for comparison in 
#' variable selection and prediction. Can be "mean" for leave-one-out models, 
#' or "orig." for the original model. 
#' @param folds the same argument in `cvrisk()`. a weight matrix with number of 
#' rows equal to the number of observations. The number of columns corresponds 
#' to the number of cross-validation runs. Can be computed using function cv 
#' and defaults to 25 bootstrap samples.
#' @param ... Additional arguments passed on to cvrisk. see `mboost` package.
#'
#' @details
#' The function conducts a leave-one-out (LOO) analysis by re-fitting the model 
#' iteratively with one observation removed, and examining the impact on variable 
#' selection or prediction. In the case of variable selection, the function 
#' calculates the difference in selected variables between the leave-one-out 
#' model and the reference model, which can either be the original model or 
#' the mean of all leave-one-out models. If the goal is prediction, the function 
#' compares the mean of cross-validation error between the leave-one-out model 
#' and a reference model. The function also allows for parallel running by 
#' specifying the number of CPU cores.
#'
#' @return A list with three elements:
#' - `vsScore`: A numeric vector of scores indicating the influence 
#'    of each observation on variable selection, if "variableSelection" is 
#'    among the aims. `NULL` otherwise.
#' - `PredScores`: A numeric vector of scores indicating the influence of each 
#'    observation on prediction, if "prediction" is among the aims.
#'     `NULL` otherwise.
#' - `LooObj`: A list of objects resulting from the LOO analysis, including 
#'    details like selected variables, cross-validated risk, and the optimal 
#'    stopping iteration for each left-out observation.
#'
#' @examples
#' 
#' library(mboost)
#' data(golub99)
#' X <- golub99$X
#' X <- scale(X)
#' y <- golub99$y
#' set.seed(1) 
#' foldid <- sample(rep(1:10, length = nrow(X)), nrow(X)) 
#' cv <- sapply(1:max(foldid), function(x) as.numeric(foldid != x))
#' obj <- glmboost(X, y,
#'                 family = Binomial(),
#'                 control = boost_control(mstop = 140,
#'                                         nu = 0.1,
#'                                         risk = "inbag"), 
#'                                         center = FALSE)
#' drop1obj <- glmboostDrop1(obj,
#'                           nCores = 1,
#'                           fixMstop = NULL,
#'                           aims = c("variableSelection", "prediction"),
#'                           ref = "mean",
#'                           folds = cv)
#' plot_Scores(drop1obj)
#' plot_Path(drop1obj, ref = "mean")
#' 
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import mboost
#' @importFrom stats  model.weights rmultinom sd
#' @export glmboostDrop1

glmboostDrop1 <- function(obj, 
                          nCores = 1, 
                          fixMstop = NULL, 
                          aims = c("variableSelection", 
                                   "prediction"),
                          ref = c("mean", "orig."),
                          folds = NULL, 
                          ...){
  
  # argument checking
  stopifnot(inherits(obj, "glmboost"))
  
  if(is.null(folds))
    folds <- cv(model.weights(obj))
  
  # sample size
  n <- length(obj$`(weights)`)
  
  # argument matching
  # aims <- match.arg(aims, c("variableSelection", "prediction"))
  ref  <- match.arg(ref, c("mean", "orig."))
  
  # parameter for parallel running
  ncores <- pmin(parallel::detectCores() - 1, nCores)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  cvboost <- foreach(i = 0:n, .packages = "mboost") %dopar% {
                       
                       # set the deleted case with zero weight
                       wts <- rep(1, n)
                       wts[i] <- 0
                       # combine wts and original weight
                       wts <- wts * obj$`(weights)`
                       
                       # allicate new weights
                       newObj <- obj$update(weights = wts, 
                                            # oobweights = wts,
                                            risk = obj$control$risk)
                       
                       if ((!is.null(fixMstop)) & 
                           identical("variableSelection", aims)) {
                         cvm = NULL
                         optmMstop = fixMstop
                       } else {
                         cvs <- do.call(cvrisk_mboost, 
                                        args = c(list("object" = newObj,
                                                      folds = folds), 
                                                 list(...)))
                         if (is.null(fixMstop))
                           optmMstop <- mstop(cvs)
                         else
                           optmMstop <- fixMstop
                         
                         # output cvm even if needed
                         cvm <- colMeans(cvs)
                       }
                       
                       # output selected variables, even not needed
                       vs <- newObj$xselect()
                       
                       out <- list("obs" = i,
                                   "vs" = vs,
                                   "cvm" = cvm,
                                   "mstop" = optmMstop)
                       
                       return(out) 
                     }
  
  stopCluster(cl)
  
  # IPs in variable selection
  if ("variableSelection" %in% aims ){
    
    drop1vs <- lapply(cvboost, function(x) {
      
      # incase len reached the largest mstop
      len <- ifelse(x[["mstop"]] >= length(x[["vs"]]), 
                    x[["mstop"]], 
                    x[["mstop"]] + 1)
      unique(as.character(x[["vs"]][seq_len(len)]))
    })
    
    # M1 
    if (ref == "orig.") {
      vsScore <- mapply(ham, sj = drop1vs[-1], MoreArgs = list(si = drop1vs[[1]]))
      names(vsScore) <- paste0("obs", 1:n)
    }
    
    # M2 
    if (ref == "mean") {
      vsScore <- moddist(mod = drop1vs[-1], std = F, L = 1)
      names(vsScore) <- paste0("obs", 1:n)
    }
  } else {
    vsScore = NULL
  }
  
  if ("prediction" %in% aims) {
    
    cvms <- lapply(cvboost, function(x)
      x[["cvm"]][x[["mstop"]] + 1]
    )
    
    cvms <- unlist(cvms)
    
    # reference
    if (ref == "mean")
      ref <- mean(cvms[-1])
    if (ref == "orig.")
      ref <- cvms[1]
    
    PredScores <- abs(cvms - ref)[-1]
    PredScores <- PredScores / sd(PredScores)
    # names(PredScores) <- paste0("obs", 1:n)
    
  } else {
    PredScores <- NULL
  }
  
  # output
  out <- list("vsScore" = vsScore, 
              "PredScores" = PredScores,
              "LooObj" = cvboost)
  
  class(out) <- "glmboostDrop1"
  return(out)
  
}





