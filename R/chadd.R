chadd <- function(param,pos=1,add.fun,...){
	
	add.fun<-match.fun(add.fun)
	if(param$ratio==1&pos!=1)stop("when ratio is 1, pos can only be 1")
	if(pos==1){
		if(param$ratio==1)par(mar=c(5,4,2,2))
		else par(mar=c(5,4,0,0))
		par(xlog=param$is.xlog,ylog=param$is.ylog,fig=c(0,param$ratio,0,param$ratio),usr=param$usrc)
	}
	else if(pos==2){par(xlog=param$is.xlog);par(fig=c(0,param$ratio,param$ratio,1),mar=c(0,4,2,0),usr=param$usru)}
	else if(pos==4){par(ylog=param$is.ylog);par(fig=c(param$ratio,1,0,param$ratio),mar=c(5,0,0,2),usr=param$usrr)}
	else if(pos==3)par(fig=c(param$ratio,1,param$ratio,1),mar=c(0,0,4,2))
	add.fun(...)
	par(fig=c(0,1,0,1),mar=c(5,4,2,2),new=TRUE)
		par(usr=c(0,1,0,1))
	plot(1,1,col=0,axes=FALSE,xlab="",ylab="")
}


