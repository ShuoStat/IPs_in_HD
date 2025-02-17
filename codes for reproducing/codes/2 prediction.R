#- loading packages and internal functions

source("./helper.R")
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
#- prediction for a single point
#-------------------------------------------------------------------------------

# plot in main text
jpeg("../results/predPoints.jpeg", width = 8, height = 5, 
     units = "in", res = 300)

layout(matrix(1:6, nrow = 2, byrow = F))
par(mar = c(3.5, 3.5, 2, 1), oma = c(1, 1, 1, 1))

for(i in d[c(3, 4, 6)]){
  load(paste0("../output/mstop_", i))
  n <- length(mstops$mstops)
  obj <- lapply(mstops$mstops, `[[`, "cvboost")
  mstop <- unlist(lapply(mstops$mstops, `[[`, "mstop"))
  plot_predPoint(obj, mstops = mstop, top = 5, ylim = c(0, 14), ref = "orig", 
                 name = paste0(gsub(".RData", "", i), "(adaptive)"))
  
  plot_predCorr(obj, mstops = mstop, top = 3, ref = "mean",
                name = paste0(gsub(".RData", "", i)))
  
}

dev.off()


# plot in appendix 
jpeg("../results/predPoints_Appendix.jpeg", width = 8, height = 5, 
     units = "in", res = 300)

layout(matrix(1:6, nrow = 2, byrow = F))
par(mar = c(3, 3, 2, 1), oma = c(1, 1, 1, 1))

for(i in d[-c(3, 4, 6)]){
  
  load(paste0("../output/mstop_", i))
  n <- length(mstops$mstops)
  obj <- lapply(mstops$mstops, `[[`, "cvboost")
  mstop <- unlist(lapply(mstops$mstops, `[[`, "mstop"))
  
  plot_predPoint(obj, mstops = mstop, top = 5, ylim = c(0, 12), ref = "orig",
                 name = paste0(gsub(".RData", "", i), "(adaptive)"))
  
  plot_predCorr(obj, mstops = mstop, top = 3, ref = "mean", name = paste0(gsub(".RData", "", i)))
  
}

dev.off()

#-------------------------------------------------------------------------------
#- prediction for a sets of variables
#-------------------------------------------------------------------------------

# predPath in main text

jpeg("../results/predPath.jpeg", width = 8, height = 5, units = "in", res = 300)
layout(matrix(1:6, nrow = 2, byrow = F))
par(mar = c(3, 3, 2, 1), oma = c(1, 1, 1, 1))

for(i in d[c(3, 4, 6)]){
  
  load(paste0("../output/mstop_", i))
  obj <- lapply(mstops$mstops, `[[`, "cvboost")
  plot_predPath(obj, name = gsub(".RData", "", i))
  plot_predScore(obj, ylim = c(0, 13), name = gsub(".RData", "", i))
}

dev.off()

# predPath in appendix
jpeg("../results/predPath_appendix.jpeg", width = 8, height = 5, units = "in", 
     res = 300)
layout(matrix(1:6, nrow = 2, byrow = F))
par(mar = c(3, 3, 2, 1), oma = c(1, 1, 1, 1))

for(i in d[-c(3, 4, 6)]){
  
  load(paste0("../output/mstop_", i))
  obj <- lapply(mstops$mstops, `[[`, "cvboost")
  plot_predPath(obj, name = gsub(".RData", "", i))
  plot_predScore(obj, ylim = c(0, 13), name = gsub(".RData", "", i))
}

dev.off()


#------------------------------------------------------------------------------



