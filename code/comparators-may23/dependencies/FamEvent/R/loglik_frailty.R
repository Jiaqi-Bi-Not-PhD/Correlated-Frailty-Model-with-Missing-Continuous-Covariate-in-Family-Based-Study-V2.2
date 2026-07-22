loglik_frailty<- function(X, Y, theta, cuts=NULL, nbase, data, design, base.dist, frailty.dist, agemin, vec=FALSE)
{
  
if(!design %in% c("pop", "pop+"))  stop("Frailty model is only available for POP or POP+ design.")

    
if(base.dist=="lognormal") bparms <- c(theta[1], exp(theta[2]))
else bparms <- exp(theta[1:nbase])

nX <- dim(X)[2]
xbeta <- c(X%*%theta[(nbase+1):(nbase+nX)])
kappa <- exp(theta[length(theta)])

time0 <- Y[,1] - agemin
cuts0 <- cuts - agemin
status <- Y[,2]

ip <- data$proband == 1
ip_fam <- aggregate(data$proband, list(data$famID.byuser), sum)[,2] # indicates if family has proband or not
i_ap <- ifelse(data$proband == 1 & data$time < data$currentage, 1, 0)[ip] ### indicates if proband is affected ###

wt <- data$weight
wt_fam <- wt[!duplicated(data$famID.byuser)]

bhaz <- hazards(base.dist, time0, bparms, cuts=cuts0)
bcumhaz <- cumhaz(base.dist, time0, bparms, cuts=cuts0)
H <- bcumhaz*exp(xbeta)
logh <- log(bhaz) + xbeta
loglik <-  wt * (status*logh)
df <- data$df[!duplicated(data$famID.byuser)]
s <- aggregate(H, list(data$famID.byuser), sum)[,2]
logdL <- wt_fam*log( dlaplace(frailty.dist, g=s, d=df, k=kappa) )


# Ascertainment correction (design = pop, pop+)
# families with no proband, logasc=0
# families with affected proband, logasc = log(1-P(T>a))
# families with unaffected proband, logasc = log(P(T>a))

  cagep <- data$currentage[ip]-agemin
  xbeta.p <- xbeta[ip]
  bcumhaz.p <- cumhaz(base.dist, cagep, bparms, cuts=cuts0)
  laplace.p <- laplace(frailty.dist, bcumhaz.p*exp(xbeta.p), kappa)
  logasc.p <- ifelse(i_ap==1, log(1-laplace.p), log(laplace.p))
  logasc <- wt_fam*ifelse(ip_fam==0, 0, logasc.p)
  #ip_fam[ip_fam!=0] <- i_ap
  
  logasc[logasc == -Inf] <- 0
  sloglik <- sum(loglik[loglik!=-Inf], na.rm=T) + sum(logdL[logdL!=-Inf], na.rm=T) - sum(logasc[logasc!=-Inf], na.rm=T)
  loglik[ip] <- loglik[ip] + logdL[ip_fam==1] - logasc[ip_fam==1]
  
  #print(c(theta, -sloglik))
  #print(c(theta, -sloglik, sum(logdL), sum(logasc)))
  if(vec) return(-loglik)
  else  return(-sloglik)
  
}

