 chplot <-
function (formula, data = parent.frame(), chull = TRUE, clevel = 0.95, 
    band.power = 0.2, mar.den = FALSE, descriptives = "mean.sd", 
    dlevel = 0.68, bw = FALSE, ratio = 0.75, plot.points=FALSE, 
    log = "", xlab, ylab, col, lty, legend=TRUE, ...) 
{
    if(is.logical(legend)){
    	lglegend <- TRUE
    	if(legend)nolegend <- FALSE
    	else nolegend <- TRUE
    	}
    else {
    nolegend <- FALSE
    lglegend <- FALSE
    }
   if(lglegend) legend <- list(title=NULL, area.in=NULL,pos=NULL,cex=NULL,bty=NULL,lwd=NULL,col=NULL,lty=NULL,fill=NULL)
 
    legend.title <- legend$title
    area.in <- legend$area.in
    legend.pos <- legend$pos
    legend$title <- legend$area.in <- legend$pos <- NULL

    #if (missing(addpoints)) 
    #    addpoints <- addpoints.control(include=FALSE)
      
    form <- latticeParseFormula(formula, data)
    if (length(names(form$condition)) > 1) 
        stop("Only 1 grouping variable is allowed")
    else if (is.null(form$condition)) {
        groups <- rep(1, length(form$left))
        nogroups <- TRUE
        nolegend <- TRUE
    }
    else {
        groups <- form$condition[[1]]
        if (is.null(legend.title)) 
            legend.title <- names(form$condition)
        nogroups <- FALSE
    }
    if (missing(xlab)) 
        xlabel <- form$right.name
    else xlabel <- xlab
    x <- form$right
    if (missing(ylab)) 
        ylabel <- form$left.name
    else ylabel <- ylab
    y <- form$left
    if (is.null(area.in))  area.in <- chull
    
    par.old <- par(c("mfrow", "mar", "font", "fig", "usr"))
    on.exit(par(par.old))
    require(ellipse)
    require(KernSmooth)
    faktor <- factor(groups)
    nlev <- nlevels(faktor)
    if (bw) {
        if (!missing(lty)) {
            ltyps <- lty
            if (length(ltyps) != nlev) 
                stop("length of lty is not equal to ", nlev)
        }
        else if (nogroups) 
            ltyps <- 1
        else ltyps <- 2:(nlev + 1)
        colrs <- rep(1, nlev)
    }
    else {
        if (!missing(col)) {
            colrs <- col
            if (length(col) != nlev) 
                stop("length of col is not equal to ", nlev)
        }
        else if (nogroups) 
            colrs <- 1
        else colrs <- 2:(nlev + 1)
        if (!missing(lty)) {
            ltyps <- lty
            if (length(ltyps) != nlev) 
                stop("length of lty is not equal to ", nlev)
        }
        else ltyps <- rep(1, nlev)
    }
    
    na.check <- as.logical((!is.na(faktor)) * (!is.na(x)) * (!is.na(y)))
    if (sum(!na.check) > 0) {
        x <- x[na.check]
        y <- y[na.check]
        faktor <- faktor[na.check]
        warning("Missing values excluded listwise.", call. = FALSE)
    }
    my.freq <- function(i, x) {
        x <- x[as.integer(faktor) == i]
        my.hist <- hist(x, plot = FALSE)
        cbind(c(my.hist$breaks[1], my.hist$mids, my.hist$breaks[length(my.hist$breaks)]), 
            c(0, my.hist$counts/diff(my.hist$breaks)/sum(my.hist$counts), 
                0), c(0, diff(my.hist$breaks), 0))
    }
    my.dense <- function(i, x, is.log) {
        x <- x[as.integer(faktor) == i]
        cbind(density(x)$x, density(x)$y)
    }
    maxi <- function(matrika, i, m = "max") {
        m.fun <- get(m)
        m.fun(matrika[, i])
    }
    area.fun <- function(x) {
        y <- c(x[, 2], x[1, 2])
        x <- c(x[, 1], x[1, 1])
        i <- 2:length(x)
        return(0.5 * sum(x[i] * y[i - 1] - x[i - 1] * y[i]))
    }
    my.polygon <- function(d, data, faktor) {
        points <- (data[as.integer(faktor) == d, ])[chull(data[as.integer(faktor) == 
            d, 1:2]), 1:2]
        polygon(points, lty = ltyps[d], border = colrs[d])
        area.fun(points)/length(data[as.integer(faktor) == d])
    }
    my.line.list <- function(d, data, faktor, clevel, bandwidth) {
        est <- bkde2D(data[as.integer(faktor) == d, ], bandwidth = nrow(data[as.integer(faktor) == 
            d, ])^(-bandwidth))
        line.list <- contourLines(est$x1, est$x2, est$fhat, nlevels = 1, 
            levels = 1 - clevel)
        line.list <- cbind(line.list[[1]][[2]], line.list[[1]][[3]])
        lines(line.list, lty = ltyps[d], col = colrs[d])
        area.fun(line.list)/length(data[as.integer(faktor) == d, ])
    }
    my.lines.mean <- function(d, x, y, faktor, sd, dlevel) {
        msd <- abs(qnorm((1 - dlevel)/2))
        podx <- x[as.integer(faktor) == d]
        pody <- y[as.integer(faktor) == d]
        if (sd) {
            lines(c(mean(podx) - msd * sqrt(var(podx)), mean(podx) + 
                msd * sqrt(var(podx))), c(mean(pody), mean(pody)), 
                lty = ltyps[d], col = colrs[d])
            lines(c(mean(podx), mean(podx)), c(mean(pody) - msd * 
                sqrt(var(pody)), mean(pody) + msd * sqrt(var(pody))), 
                lty = ltyps[d], col = colrs[d])
        }
        else {
            lines(c(mean(podx) - msd * sqrt(var(podx)/length(podx)), 
                mean(podx) + msd * sqrt(var(podx)/length(podx))), 
                c(mean(pody), mean(pody)), lty = ltyps[d], col = colrs[d])
            lines(c(mean(podx), mean(podx)), c(mean(pody) - msd * 
                sqrt(var(pody)/length(pody)), mean(pody) + msd * 
                sqrt(var(pody)/length(pody))), lty = ltyps[d], 
                col = colrs[d])
        }
    }
    my.lines.median <- function(d, x, y, faktor) {
        podx <- x[as.integer(faktor) == d]
        pody <- y[as.integer(faktor) == d]
        lines(c(quantile(podx, 0.25), quantile(podx, 0.75)), 
            c(median(pody), median(pody)), lty = ltyps[d], col = colrs[d])
        lines(c(median(podx), median(podx)), c(quantile(pody, 
            0.25), quantile(pody, 0.75)), lty = ltyps[d], col = colrs[d])
    }
    my.ellipse <- function(d, x, y, faktor, level) {
        podx <- x[as.integer(faktor) == d]
        pody <- y[as.integer(faktor) == d]
        lines(ellipse(cor(cbind(podx, pody)), centre = c(mean(podx), 
        mean(pody)), scale = c(sqrt(var(podx)/length(podx)), sqrt(var(pody)/length(pody))), 
        level = level), lty = ltyps[d], col = colrs[d])
    }
    par(fig = c(0, 1, 0, 1))
    par(mar = c(5, 4, 2, 2))
    if (ratio < 1) {
        frame()
        par(fig = c(0, ratio, 0, ratio), mar = c(5, 4, 0, 0))
    }
    kode <- 1:nlevels(faktor)
        
    args <- list(x = x, y = y, xlab = xlabel, ylab = ylabel, log = log,...)
     
    if(bw){
    	if(is.null(args$pch))args$pch <- as.integer(faktor)+1
    }
    else{
    	if(is.null(args$col))args$col <- as.integer(faktor)+1 
    }
    if(!plot.points) args$type <- "n"
       
    do.call("plot", args)
  
    usr <- par("usr")
    is.xlog <- par("xlog")
    is.ylog <- par("ylog")
    if(is.xlog){
	na.check <- x>0
	if (sum(!na.check) > 0) {
	    x <- x[na.check]
	    y <- y[na.check]
	    faktor <- faktor[na.check]
	    warning("Cases with x<=0 excluded listwise.", call. = FALSE)
	}
    }
    if(is.ylog){
    	na.check <- y>0
    	if (sum(!na.check) > 0) {
    	    x <- x[na.check]
    	    y <- y[na.check]
    	    faktor <- faktor[na.check]
    	    warning("Cases with y<=0 excluded listwise.", call. = FALSE)
    	}
    }
    if (chull) 
        area <- unlist(lapply(kode, my.polygon, cbind(x, y), 
            faktor))
    else area <- unlist(lapply(kode, my.line.list, cbind(x, y), 
        faktor, clevel, band.power))
    areaout <- area
    if (area.in) 
        area <- format(area, digits = 2)
    if (descriptives == "mean.sd") 
        lapply(kode, my.lines.mean, x, y, faktor, sd = TRUE, 
            dlevel)
    else if (descriptives == "mean.se") 
        lapply(kode, my.lines.mean, x, y, faktor, sd = FALSE, 
            dlevel)
    else if (descriptives == "median") 
        lapply(kode, my.lines.median, x, y, faktor)
    else if (descriptives == "ellipse") 
        lapply(kode, my.ellipse, x, y, faktor, dlevel)
    if (mar.den) 
        lister <- get("my.dense")
    else lister <- get("my.freq")
    if (is.xlog) 
        x <- log(x)
    if (ratio < 1) {
        par(fig = c(0, ratio, ratio, 1), mar = c(0, 4, 2, 0), 
            new = TRUE)
        lista <- lapply(1:nlevels(faktor), lister, x)
        maxbreak <- 1
        if (!mar.den) 
            maxbreak <- max(unlist(lapply(lista, maxi, 3)))
        height <- max(unlist(lapply(lista, maxi, 2))) * maxbreak
        args2 <- list(x=c(min(x), max(x)),y=c(0, height),axes=FALSE,
        	 xlab = "", ylab = "", type = "n")
        if (is.xlog) 
            args2$log <- "x"
        if(!is.null(args$xlim))args2$xlim <- args$xlim
        do.call("plot",args2)
        for (i in 1:length(lista)) {
            points(lista[[i]][, 1], lista[[i]][, 2] * maxbreak, 
                type = "l", lty = ltyps[i], col = colrs[i])
        }
        usru <- par("usr")
        ticks <- round(height/3, 2)
        ticks2 <- round(height * 2/3, 2)
        height1 <- round(height, 2)
        axis(2, at = c(0, ticks, ticks2, height1), labels = c("", 
            ticks, ticks2, height1), cex.axis = 0.75)
        if (is.ylog) 
            y <- log(y)
        par(fig = c(ratio, 1, 0, ratio), mar = c(5, 0, 0, 2), 
            new = TRUE)
        lista <- lapply(1:nlevels(faktor), lister, y)
        maxbreak <- 1
        if (!mar.den) 
            maxbreak <- max(unlist(lapply(lista, maxi, 3)))
        height <- max(unlist(lapply(lista, maxi, 2))) * maxbreak
        args2 <- list(c(0, height), c(min(y), max(y)), axes = FALSE, 
                xlab = "", ylab = "", type = "n")
	if (is.xlog) 
	      args2$log <- "y"
	if(!is.null(args$ylim))args2$ylim <- args$ylim
        do.call("plot",args2)
        for (i in 1:length(lista)) {
            points(lista[[i]][, 2] * maxbreak, lista[[i]][, 1], 
                type = "l", lty = ltyps[i], col = colrs[i])
        }
        usrr <- par("usr")
        ticks <- round(height/3, 2)
        ticks2 <- round(height * 2/3, 2)
        height1 <- round(height, 2)
        axis(1, at = c(0, ticks, ticks2, height1), labels = c("", 
            ticks, ticks2, height1), cex.axis = 0.75)
        par(fig = c(ratio, 1, ratio, 1), mar = c(0, 0, 2, 2), 
            new = TRUE)
        plot(1, 1, col = 0, axes = FALSE, xlim = c(0, 1), ylim = c(0, 
            1))
    }
    max.len <- max(nchar(levels(faktor)))
    if (max.len < 10) 
        max.len1 <- -max.len * 1.1
    else max.len1 <- format(-max.len * 1.01, nsmall = 2)
    sprintf.format <- paste("%", max.len1, "s", sep = "")
    cum <- function(line, sprintf.format) {
        temp2 <- paste("(", line[1], ")", sep = "")
        temp1 <- sprintf(sprintf.format, line[2])
        paste(temp1, temp2, sep = " ")
    }
    
    if (!nolegend){
    	
        if (area.in) 
            par(font = 10)
        out <- 2
        extra <- 1
        if (!is.null(legend.pos)) {
            out <- as.numeric(legend.pos == "out")
        }
        if (area.in) 
            legend.text <- apply(cbind(area, levels(faktor)), 
                1, cum, sprintf.format)
        else legend.text <- levels(faktor)
        l.args <- list(legend=legend.text)
        if (ratio > 0.75 & out != 1 | out == 0) {
            par(fig = c(0, ratio, 0, ratio), mar = c(5, 4, 0, 
                0))
            xy <- locator(1)
            if (!is.null(legend$cex)) 
                cex.size <- legend$cex
            else cex.size <- 1
            if (is.null(legend$bty)) l.args$bty <- "n"
	    if (is.null(legend$lwd)) l.args$lwd <- 2*cex.size
	    if (is.null(legend$col)) l.args$col <- colrs
            if (is.null(legend$lty)) l.args$lty <- ltyps
            if (!is.null(legend$fill)) {
            	l.args$fill <- colrs
            	legend$fill <- NULL
            	l.args$lty <- NULL
            	l.args$lwd <- NULL
            }
        }
        else {
            wz <- list(x = -0.05, y = 0.95)
            xy <- list(x = 0, y = 0.9)
           
           
                if (area.in) 
                  cex.size <- (2/3)^((ratio - 0.65) * 10) * 12/(max.len + 
                    6) * (0.93)
                else cex.size <- (2/3)^((ratio - 0.65) * 10) * 
                  10/max.len * (0.93)
                if (cex.size > 1) 
                  cex.size <- 1
           
             if (!is.null(legend$cex)) 
                cex.size <- legend$cex*cex.size
            l.args$cex <- cex.size
            text(wz, legend.title, pos = 4, cex = cex.size +  0.1)
            if (is.null(legend$bty)) l.args$bty <- "n"
            if (is.null(legend$lwd)) l.args$lwd <- 2*cex.size
            if (is.null(legend$col)) l.args$col <- colrs
            if (is.null(legend$lty)) l.args$lty <- ltyps
            if (!is.null(legend$fill)) {
		l.args$fill <- colrs
		legend$fill <- NULL
		l.args$lty <- NULL
		l.args$lwd <- NULL
            }
        }
        l.args$x <-xy
        if(!lglegend) l.args[names(legend)] <- legend
        do.call("legend",l.args)
    }
    par(fig = c(0, 1, 0, 1), mar = c(5, 4, 2, 2), new = TRUE)
    par(usr = c(0, 1, 0, 1))
    plot(1, 1, type="n", axes = FALSE, xlab = "", ylab = "")
    if (ratio == 1) 
        usru <- usrr <- NA
    out <- list(area=areaout,usrc = usr, usru = c(usr[1:2], usru[3:4]), usrr = c(usrr[1:2], 
        usr[3:4]), ratio = ratio, is.xlog = is.xlog, is.ylog = is.ylog)
}
