#- loading packages and internal functions

source("./fun.R")
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
#- prediction
#-------------------------------------------------------------------------------

jpeg("./output/predPath.jpeg", width = 9, height = 5, units = "in", res = 300)
layout(matrix(1:6, nrow = 2, byrow = F))
par(mar = c(3, 3, 2, 1), oma = c(1, 1, 1, 1))

for(i in d[c(3, 4, 6)]){
  
  load(paste0("../output/mstop_", i))
  obj <- lapply(mstops$mstops, `[[`, "cvboost")
  plot_predPath(obj, name = gsub(".RData", "", i))
  plot_predScore(obj, ylim = c(0, 13), name = gsub(".RData", "", i))
}

dev.off()

#------------------------------------------------------------------------------
