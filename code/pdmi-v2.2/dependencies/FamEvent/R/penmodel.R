# frailty.dist option 
penmodel <- function(formula, cluster="famID", gvar="mgene", parms, cuts=NULL, data, design="pop", base.dist="Weibull", frailty.dist="none", agemin=NULL, robust=FALSE){
  
  if(!frailty.dist%in%c("none", "gamma", "lognormal")) stop("frailty.dist should be one of \"none\", \"gamma\", or \"lognormal\".")
  nfp <- ifelse(frailty.dist=="none", 0, 1)
  
#  if(any(is.na(data[, gvar]))) cat("Individuals with missing genetic information were removed.") #stop("data include missing genetic information, use penmodelEM function.")
  options(na.action='na.omit')
  
  if(is.null(agemin)){
    agemin <- attr(data, "agemin")
    if(is.null(agemin)) {
          agemin <- 0
          warning("agemin = 0 was used or assign agemin to attr(data, \"agemin\").")
    }
  }
  
  
  if(sum(data$time <=  agemin, na.rm = TRUE) > 0) cat("Individuals with time < agemin (", agemin,") were removed.\n")
  
  data <- data[data$time >=  agemin, ]
  data$famID.byuser <- data[, cluster]
  
  x.names <- attr(terms(formula), "term.labels")
  i.missing <- apply(data[, x.names], 1, anyNA) # indicates if any x variable is missing
  newdata <- data[!i.missing,]

  if( sum(i.missing) >0 ) cat("Individuals with missing covariates were removed.\n")
  m <- model.frame(formula, data) # missing data were removed
  
  Terms <- attr(m, "terms")
  Y <- model.extract(m, "response")
  
  if (!inherits(Y, "Surv")) stop("Response must be a survival object.")
  type <- attr(Y, "type")
  if (type == "counting") stop("start-stop type Surv objects are not supported.")
  if (type == "mright" || type == "mcounting") stop("multi-state survival is not supported.")

  if(base.dist=="piecewise"){
    if(is.null(cuts)) stop("The cuts should be specified")
    if(any(cuts > max(Y[,1]) | cuts < min(Y[,1]))) stop("Some value(s) of the cuts are beyond the range.")
  }


  X <- model.matrix(Terms, m)

  nvar <- ncol(X)-1
  var.names <- colnames(X)[-1]
  if(nvar==1) X <- matrix(X[,-1])
  else X <- X[,-1]

  #number of parameters for baseline
  nbase <- ifelse(base.dist=="logBurr", 3, ifelse(base.dist=="piecewise", length(cuts)+1, 2)) 

  #nk <- ifelse(is.null(frailty.dist),0,1)
  colnames(X) <- var.names
  vbeta <- parms[(nbase+1):(nbase+nvar)]
  kappa <- parms[(nbase+nvar+1)]
#  if(length(vbeta) != nvar) stop("The size of parms is incorrect.")
  if(length(parms) != (nvar+nbase+nfp) ) stop("The size of parms is incorrect.")
  
  if(is.null(newdata$weight)) newdata$weight <- 1
 
  if(base.dist=="lognormal"){
  	if(parms[2] <= 0) stop("parms[2] has to be > 0")
  	else parms[1] <- exp(parms[1])
  } 
  else if(any(parms[1:nbase]<=0)) stop("All baseline parameters should be > 0")
  
  if(frailty.dist=="none"){
    est1 <- optim(c(log(parms[1:nbase]), vbeta), loglik_ind, X=X, Y=Y, cuts=cuts, 
                  nbase=nbase, data=newdata, design=design, base.dist=base.dist, 
                  agemin=agemin, control = list(maxit = 50000), hessian=TRUE)
  }
  else{ # frailty model
    if(is.null(newdata$df)){
      df <- aggregate(Y[,2], list(newdata$famID), sum)[,2]
      fsize <- aggregate(Y[,2], list(newdata$famID), length)[,2]
      newdata$df <- rep(df, fsize)
    }
    
    est1 <- optim(c(log(parms[1:nbase]), vbeta, log(kappa)), loglik_frailty, X=X, Y=Y, 
                  cuts=cuts, nbase=nbase, data=newdata, design=design, base.dist=base.dist, 
                  frailty.dist=frailty.dist, agemin=agemin, vec=FALSE,
                  control = list(maxit = 50000), hessian=TRUE)
  }
    
  logLik <- -est1$value
    EST <- est1$par
    H <- est1$hessian
    Var <- try(solve(H), TRUE)
    ## Convergence check
    Conv <- est1$convergence

  if(Conv != 0) stop(paste("Model did not converge.\n Convergence code "), Conv)
  if(!is.null(attr(Var,"class"))) stop("Model did not converge.\n  Try again with different initial values.")
  else{ 
    bparms.name <- c("log.lambda","log.rho", "log.eta")
    if(base.dist=="lognormal") bparms.name[1] <- "lambda" 
    else if(base.dist=="piecewise") bparms.name <- paste0("log.q", 1:nbase,")")
    
    parms.cov <- Var
    parms.se <- sqrt(diag(parms.cov))
    parms.cov.robust <- parms.se.robust <- NULL
    if(frailty.dist=="none") parms.name <- c(bparms.name[1:nbase], colnames(X))
    else parms.name <- c(bparms.name[1:nbase], colnames(X), "log.kappa")
    
    names(EST) <- names(parms.se)  <- rownames(parms.cov) <- colnames(parms.cov) <- parms.name
    
    if(robust){
      if(frailty.dist=="none")  grad <- jacobian(loglik_ind, est1$par, X=X, Y=Y, cuts=cuts, nbase=nbase, data=newdata, design=design, base.dist=base.dist, agemin=agemin, vec=TRUE)
      else grad <- jacobian(loglik_frailty, est1$par, X=X, Y=Y, cuts=cuts, nbase=nbase, data=newdata, design=design, base.dist=base.dist, frailty.dist=frailty.dist, agemin=agemin, vec=TRUE)
      Jscore <- t(grad)%*%grad
 		 parms.cov.robust <- Var%*%(Jscore)%*%Var
 		 parms.se.robust <- sqrt(diag(parms.cov.robust))
 		 rownames(parms.cov.robust) <- colnames(parms.cov.robust) <- parms.name
    }
  }
    
    aic = 2*length(EST) - 2*logLik
    
    out <- list(estimates = EST, varcov = parms.cov, varcov.robust = parms.cov.robust, se = parms.se, se.robust = parms.se.robust,logLik = logLik, AIC=aic)  
    class(out) <- "penmodel"
    attr(out, "design") <- design
    attr(out, "base.dist") <- base.dist
    attr(out, "frailty.dist") <- frailty.dist
    attr(out, "agemin") <- agemin
    attr(out, "cuts") <- cuts
    attr(out, "nbase") <- nbase
    attr(out, "data") <- newdata
    attr(out, "robust") <- robust
    attr(out, "formula") <- formula
    attr(out, "X") <- X
    attr(out, "Y") <- Y
    invisible(out)
    
}#end
  