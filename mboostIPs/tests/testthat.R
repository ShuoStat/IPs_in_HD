
library(testthat)
library(mboostIPs)
library(mboost)
data(golub99)
X <- golub99$X
X <- scale(X)
y <- golub99$y
set.seed(1)
foldid <- sample(rep(1:10, length = nrow(X)), nrow(X))
cv <- sapply(1:max(foldid), function(x) as.numeric(foldid != x))
obj <- mboost::glmboost(X, y,
                family = mboost::Binomial(),
                control = mboost::boost_control(mstop = 140,
                                                nu = 0.1,
                                                risk = "inbag"),
                center = F)

# drop1obj.mean <- glmboostDrop1(obj,
#                                nCores = 8,
#                                fixMstop = NULL,
#                                aims = c("variableSelection", "prediction"),
#                                ref = "mean",
#                                folds = cv)
# 
# drop1obj.orig <- glmboostDrop1(obj,
#                                nCores = 8,
#                                fixMstop = NULL,
#                                aims = c("variableSelection", "prediction"),
#                                ref = "orig.",
#                                folds = cv)
# 
# drop1obj.mean.fixed <- glmboostDrop1(obj,
#                                      nCores = 8,
#                                      fixMstop = 100,
#                                      aims = c("variableSelection", "prediction"),
#                                      ref = "mean",
#                                      folds = cv)
# 
# 
# drop1obj.orig.fixed <- glmboostDrop1(obj,
#                                      nCores = 8,
#                                      fixMstop = 100,
#                                      aims = c("variableSelection", "prediction"),
#                                      ref = "orig.",
#                                      folds = cv)
# 
# #
# save(list = c("drop1obj.mean", "drop1obj.orig",
#               "drop1obj.mean.fixed", "drop1obj.orig.fixed"),
#      file = "./testthat/correctResults.RData")

#

load("./testthat/correctResults.RData")
drop1obj.mean_test <- glmboostDrop1(obj,
                                    nCores = 1,
                                    fixMstop = NULL,
                                    aims = c("variableSelection", "prediction"),
                                    ref = "mean",
                                    folds = cv)

drop1obj.orig_test <- glmboostDrop1(obj,
                                    nCores = 1,
                                    fixMstop = NULL,
                                    aims = c("variableSelection", "prediction"),
                                    ref = "orig.",
                                    folds = cv)

drop1obj.mean.fixed_test <- glmboostDrop1(obj,
                                          nCores = 1,
                                          fixMstop = 100,
                                          aims = c("variableSelection", "prediction"),
                                          ref = "mean",
                                          folds = cv)


drop1obj.orig.fixed_test <- glmboostDrop1(obj,
                                          nCores = 1,
                                          fixMstop = 100,
                                          aims = c("variableSelection", "prediction"),
                                          ref = "orig.",
                                          folds = cv)

test_that("test glmboostDrop1 funciton", {
    expect_equal(drop1obj.mean, drop1obj.mean_test)
    expect_equal(drop1obj.orig, drop1obj.orig_test)
    expect_equal(drop1obj.mean.fixed, drop1obj.mean.fixed_test)
    expect_equal(drop1obj.orig.fixed, drop1obj.orig.fixed_test)
})








