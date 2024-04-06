#-------------------------------------------------------------------------------
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
                     refline = 0, ylab = ""){

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
  mtext("observations", side = 1, line = 2.2)

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





