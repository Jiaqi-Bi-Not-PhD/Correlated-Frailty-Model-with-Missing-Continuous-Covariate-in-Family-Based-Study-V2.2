inv2.surv <- function(val, base.dist, parms){
	out<-try(uniroot(surv.dist, lower=0,upper=100000, base.dist=base.dist, parms=parms, alpha=val[2], xbeta=val[1], res=val[3])$root)
  if(is.null(attr(out,"class"))) return(out)
	else print(c(parms, val))
	}
