legend.control <- function(include=TRUE,area.in,pos,cex,bty,title){
	if(missing(area.in))area.in<-NA
	else if(!is.logical(area.in))
		stop("value of legend 'area.in' must be logical")
	if(missing(cex))cex<-NA
	else if (!is.numeric(cex) || cex <= 0) 
		stop("value of legend cex must be > 0")
	if(missing(pos))pos<-NA
	else if(pos!="in"&pos!="out")
		stop("value of legend 'pos' must be either 'in' or 'out'")
	if(missing(bty))bty<-NA
	if(missing(title))title<-NA
	list(include=include,area.in=area.in,pos=pos,cex=cex,bty=bty,title=title)
}
