chplot <-function(formula,data=parent.frame(),chull=TRUE,clevel=0.95,band.power=.2,mar.den=FALSE,descriptives="mean.sd",dlevel=.68,bw=FALSE,ratio=.75,plot.points=FALSE,log="",xlab,ylab,col,lty,legend,...){
        
    if(missing(legend))legend <- legend.control(...)
    
    form <- latticeParseFormula(formula, data)	
    
    if(length(names(form$condition))>1)stop("Only 1 grouping variable is allowed")
    else if(is.null(form$condition)){
    	groups <- rep(1,length(form$left))
    	nogroups <- TRUE
    	legend$include <- FALSE
    }
    else{    
	    groups <- form$condition[[1]]	
	    if(is.na(legend$title)) legend.title <- names(form$condition)
	    else legend.title <- legend$title
	    nogroups <- FALSE
    }
    
    
        
    if(missing(xlab)) xlabel <- form$right.name
    else xlabel<-xlab
    x <- form$right
    
    if(missing(ylab)) ylabel<- form$left.name
    else ylabel<-ylab
    y <- form$left
    
    
    if(is.na(legend$area.in))area.in<-chull
    else area.in <- legend$area.in
        
    par.old<-par(c("mfrow","mar","font","fig","usr"))
    on.exit(par(par.old))
    

    require(ellipse)
    require(KernSmooth)
    
 
    faktor<-factor(groups)
    nlev <- nlevels(faktor)
    
    
    if(bw){
    	if(!missing(lty)){
    		ltyps <- lty
    		if(length(ltyps)!=nlev)stop("length of lty is not equal to ",nlev)
    	}
    	else if(nogroups) ltyps <- 1
    	else ltyps <- 2:(nlev+1)
    	colrs <- rep(1,nlev)
    }
    else{
    	if(!missing(col)){
    		colrs <- col
    		if(length(col)!=nlev)stop("length of col is not equal to ",nlev)
    	}
    	else if(nogroups) colrs <- 1
    	else colrs <- 2:(nlev+1)
    	if(!missing(lty)){
    		ltyps <- lty
    		if(length(ltyps)!=nlev)stop("length of lty is not equal to ",nlev)
    	}
    	else ltyps <- rep(1,nlev)
    }
    

    na.check<-as.logical((!is.na(faktor))*(!is.na(x))*(!is.na(y)))
    if(sum(!na.check)>0) {
	x<-x[na.check]
	y<-y[na.check]
	faktor<-faktor[na.check]
	warning("Missing values excluded listwise.",call.=FALSE)
    }


    my.freq <- function(i,x){
	x<-x[as.integer(faktor)==i]
	my.hist<-hist(x,plot=FALSE)
	cbind(c(my.hist$breaks[1],my.hist$mids,my.hist$breaks[length(my.hist$breaks)]),c(0,my.hist$counts/diff(my.hist$breaks)/sum(my.hist$counts),0),c(0,diff(my.hist$breaks),0))
    }

    my.dense <- function(i,x,is.log){
	x<-x[as.integer(faktor)==i]
	cbind(density(x)$x,density(x)$y)
    }

    maxi <- function(matrika,i,m="max"){
	m.fun<-get(m)
	m.fun(matrika[,i])
    }

    area.fun <- function(x){
        y<-c(x[,2],x[1,2])
	x<-c(x[,1],x[1,1])
        i <- 2:length(x)
        return(0.5 * sum(x[i] * y[i - 1] - x[i - 1] * y[i]))
    }
	
   	
    my.polygon <- function(d,data,faktor,area.in){
	points<-(data[as.integer(faktor)==d,])[chull(data[as.integer(faktor)==d,1:2]),1:2]
	polygon(points,lty=ltyps[d],border=colrs[d])
	if(area.in)area.fun(points)/length(data[as.integer(faktor)==d])
    }

    my.line.list <- function(d,data,faktor,clevel,bandwidth,area.in){
	est<-bkde2D(data[as.integer(faktor)==d,],bandwidth=nrow(data[as.integer(faktor)==d,])^(-bandwidth))
	line.list<-contourLines(est$x1,est$x2,est$fhat,nlevels=1,levels=1-clevel)
	line.list<-cbind(line.list[[1]][[2]],line.list[[1]][[3]])
	lines(line.list,lty=ltyps[d],col=colrs[d])
	if(area.in)area.fun(line.list)/length(data[as.integer(faktor)==d,])
    }    


    my.lines.mean <- function(d,x,y,faktor,sd,dlevel){
        msd<-abs(qnorm((1-dlevel)/2))
	podx<-x[as.integer(faktor)==d]
	pody<-y[as.integer(faktor)==d]
	if(sd){
	    lines(c(mean(podx)-msd*sqrt(var(podx)),mean(podx)+msd*sqrt(var(podx))),c(mean(pody),mean(pody)),lty=ltyps[d],col=colrs[d])
	    lines(c(mean(podx),mean(podx)),c(mean(pody)-msd*sqrt(var(pody)),mean(pody)+msd*sqrt(var(pody))),lty=ltyps[d],col=colrs[d])
	}
	else{
	    lines(c(mean(podx)-msd*sqrt(var(podx)/length(podx)),mean(podx)+msd*sqrt(var(podx)/length(podx))),c(mean(pody),mean(pody)),lty=ltyps[d],col=colrs[d])
	    lines(c(mean(podx),mean(podx)),c(mean(pody)-msd*sqrt(var(pody)/length(pody)),mean(pody)+msd*sqrt(var(pody)/length(pody))),lty=ltyps[d],col=colrs[d])
	}
    }

    my.lines.median <- function(d,x,y,faktor){
	podx<-x[as.integer(faktor)==d]
	pody<-y[as.integer(faktor)==d]
	lines(c(quantile(podx,.25),quantile(podx,.75)),c(median(pody),median(pody)),lty=ltyps[d],col=colrs[d])
	lines(c(median(podx),median(podx)),c(quantile(pody,.25),quantile(pody,.75)),lty=ltyps[d],col=colrs[d])
    }



    my.ellipse <- function(d,x,y,faktor,level){
	podx<-x[as.integer(faktor)==d]
	pody<-y[as.integer(faktor)==d]
	lines(ellipse(cov(cbind(podx,pody)),centre=c(mean(podx),mean(pody)),scale=c(sqrt( 1/length(podx) ),sqrt( 1/length(pody) )),level=level),lty=ltyps[d],col=colrs[d])
    }



    
    
    par(fig=c(0,1,0,1))
    par(mar=c(5,4,2,2))
    if(ratio<1){			
    frame()
    par(fig=c(0,ratio,0,ratio),mar=c(5,4,0,0))
    }					
    kode<-1:nlevels(faktor)

    if(!plot.points) plot(x,y,type="n",xlab=xlabel,ylab=ylabel,log=log)
    else plot(x,y,col=as.integer(faktor)+1,xlab=xlabel,ylab=ylabel,log=log)
    
	    
    usr<-par("usr")
    is.xlog<-par("xlog")
    is.ylog<-par("ylog")

    

    if(chull) area<-unlist(lapply(kode,my.polygon,cbind(x,y),faktor,area.in))
    else area<-unlist(lapply(kode,my.line.list,cbind(x,y),faktor,clevel,band.power,area.in))
    

    if(area.in) area<-format(area,digits=2)


    if(descriptives=="mean.sd")	lapply(kode,my.lines.mean,x,y,faktor,sd=TRUE,dlevel)
	
    else if(descriptives=="mean.se") lapply(kode,my.lines.mean,x,y,faktor,sd=FALSE,dlevel)
	
    else if(descriptives=="median") lapply(kode,my.lines.median,x,y,faktor)
	
    else if(descriptives=="ellipse") lapply(kode,my.ellipse,x,y,faktor,dlevel)
    

    if(mar.den)lister<-get("my.dense")
    else lister<-get("my.freq")


    if(is.xlog)x<-log(x)
    if(ratio<1){				
    par(fig=c(0,ratio,ratio,1),mar=c(0,4,2,0),new=TRUE)
    lista<-lapply(1:nlevels(faktor),lister,x)
    maxbreak<-1
    if(!mar.den)maxbreak<-max(unlist(lapply(lista,maxi,3)))
    height<-max(unlist(lapply(lista,maxi,2)))*maxbreak
    if(is.xlog) plot(c(min(x),max(x)),c(0,height),axes=FALSE,xlab="",ylab="",type="n",log="x")
    else plot(c(min(x),max(x)),c(0,height),axes=FALSE,xlab="",ylab="",type="n")

    
    for(i in 1:length(lista)){
	points(lista[[i]][,1],lista[[i]][,2]*maxbreak,type="l",lty=ltyps[i],col=colrs[i])
    }

    usru<-par("usr")
    ticks <- round(height/3,2)
    ticks2 <- round(height * 2/3,2)
    height1 <- round(height,2)
    axis(2, at = c(0,ticks,ticks2, height1), labels = c("",ticks, ticks2, height1), cex.axis = 0.75)
        
    if(is.ylog)y<-log(y)
    par(fig=c(ratio,1,0,ratio),mar=c(5,0,0,2),new=TRUE)
    lista<-lapply(1:nlevels(faktor),lister,y)
    maxbreak<-1
    if(!mar.den)maxbreak<-max(unlist(lapply(lista,maxi,3)))
    height<-max(unlist(lapply(lista,maxi,2)))*maxbreak
    if(is.ylog) plot(c(0,height),c(min(y),max(y)),axes=FALSE,xlab="",ylab="",type="n",log="y")
    else plot(c(0,height),c(min(y),max(y)),axes=FALSE,xlab="",ylab="",type="n") 	

    for(i in 1:length(lista)){
	points(lista[[i]][,2]*maxbreak,lista[[i]][,1],type="l",lty=ltyps[i],col=colrs[i])

    }


    usrr <- par("usr")
    ticks <- round(height/3,2)
        ticks2 <- round(height * 2/3,2)
    height1 <- round(height,2)
    axis(1, at = c(0, ticks, ticks2, height1), labels = c("",ticks, ticks2, height1), cex.axis = 0.75)
	

    par(fig=c(ratio,1,ratio,1),mar=c(0,0,2,2),new=TRUE)
    plot(1,1,col=0,axes=FALSE,xlim=c(0,1),ylim=c(0,1))
    } 							
    max.len<- max(nchar(levels(faktor)))
    if(max.len<10) max.len1<--max.len*1.1
    else max.len1<-format(-max.len*1.01,nsmall=2)

    sprintf.format<-paste("%",max.len1,"s",sep="")

    cum<-function (line,sprintf.format) {
	temp2<-paste("(",line[1],")",sep="")
	temp1<-sprintf(sprintf.format, line[2])
	paste(temp1,temp2,sep=" ")
    }
   
    if(legend$include){
	if(area.in)par(font=10)
	out<-2
	extra<-1
	if(!is.na(legend$pos)){out<-as.numeric(legend$pos=="out")}
	if(area.in)legend.text<-apply(cbind(area,levels(faktor)),1,cum,sprintf.format)
	else legend.text<-levels(faktor)
	if(ratio>0.75&out!=1|out==0){
	    par(fig=c(0,ratio,0,ratio),mar=c(5,4,0,0))
	    xy<-locator(1)
	    if(!is.na(legend$cex))cex.size<-legend$cex
	    else cex.size<-1
	    if(!is.na(legend$bty))btype<-legend$bty
	    else btype<-"o"
	}
	else{
	    wz<-list(x=-.05,y=.95)
	    xy<-list(x=0,y=.9)
	    if(!is.na(legend$cex))cex.size<-legend$cex
	    else{
		if(area.in)cex.size<-(2/3)^((ratio-.65)*10)*12/(max.len+6)*(.93)
		else cex.size<-(2/3)^((ratio-.65)*10)*10/max.len*(.93)
		if(cex.size>1)cex.size<-1
	    }
	    text(wz,legend.title,pos=4,cex=cex.size+.1)	
	    if(!is.na(legend$bty))btype<-legend$bty
	    else btype<-"n"
	}
	legend(xy,legend=legend.text, lty=ltyps,col=colrs,ncol=1,bty=btype,cex=cex.size,lwd=2*cex.size)
    }
    par(fig=c(0,1,0,1),mar=c(5,4,2,2),new=TRUE)
    par(usr=c(0,1,0,1))
    plot(1,1,col=0,axes=FALSE,xlab="",ylab="")
    if(ratio==1)usru<-usrr<-NA
    out <- list(usrc=usr,usru=c(usr[1:2],usru[3:4]),usrr=c(usrr[1:2],usr[3:4]),ratio=ratio,is.xlog=is.xlog,is.ylog=is.ylog)
}
