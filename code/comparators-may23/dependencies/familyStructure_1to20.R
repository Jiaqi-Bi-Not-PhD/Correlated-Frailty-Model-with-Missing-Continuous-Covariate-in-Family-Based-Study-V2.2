familyStructure <- function (
    i, cumind = cumind, m.carrier = 1, variation = "none", interaction = FALSE, 
    add.x = FALSE, x.dist = NULL, x.parms = NULL, depend = NULL, 
    base.dist = "Weibull", frailty.dist = NULL, base.parms = c(0.016, 3), 
    vbeta = c(1, 1), allelefreq = 0.02, dominant.m = TRUE, dominant.s = TRUE, 
    mrate = 0, probandage = c(45, 2), agemin = 20, agemax = 100,
    ## turn on to align to female-only size distribution (fast; no pools)
    align_female_only = TRUE,
    femsize_vals  = NULL,
    femsize_probs = NULL
) {
  
  ## ---------------- helpers ----------------
  .default_fem_table <- function() {
    vals  <- c(1:15, 19, 20)
    probs <- c(0.113131313, 0.121212121, 0.147474747, 0.129292929,
               0.145454545, 0.115151515, 0.072727273, 0.042424242,
               0.026262626, 0.024242424, 0.020202020, 0.014141414,
               0.014141414, 0.006060606, 0.004040404, 0.002020202,
               0.002020202)
    setNames(probs, vals)
  }
  .sanitize_fem_table <- function(fem_vals, fem_probs, min_fem = 3L, max_fem = 31L) {
    if (is.null(fem_vals) || is.null(fem_probs)) {
      tab <- .default_fem_table()
      fem_vals  <- as.integer(names(tab))
      fem_probs <- as.numeric(tab)
    }
    keep <- (fem_vals >= min_fem) & (fem_vals <= max_fem) & (fem_probs > 0)
    if (!any(keep)) stop("No feasible female-only sizes in provided/default table.")
    fem_vals  <- as.integer(fem_vals[keep])
    fem_probs <- as.numeric(fem_probs[keep]); fem_probs <- fem_probs / sum(fem_probs)
    list(vals = fem_vals, probs = fem_probs)
  }
  .conv <- function(p, q) {
    out <- numeric(length(p) + length(q) - 1L)
    for (a in seq_along(p)) for (b in seq_along(q)) out[a+b-1L] <- out[a+b-1L] + p[a]*q[b]
    out
  }
  
  ## original discrete laws (already "multinomial"/categorical in your code)
  p_sec <- c(`2`=0.4991, `3`=0.2720, `4`=0.1482, `5`=0.0807)  # SecNum
  p_kid <- c(`2`=0.4991, `3`=0.2720, `4`=0.1482, `5`=0.0807)  # ThiNum
  
  ## precompute sum-of-children PMFs for s = 2..5 (support 2s..5s)
  pmf_children_sum <- function(s) {
    pmf <- p_kid
    if (s == 1L) {
      supp <- 2:5
      return(list(support = supp, pmf = setNames(as.numeric(pmf), as.character(supp))))
    }
    cur <- pmf
    for (t in 2:s) cur <- .conv(cur, p_kid)
    supp <- (2*s):(5*s)
    list(support = supp, pmf = setNames(as.numeric(cur), as.character(supp)))
  }
  
  ## sample structure conditional on Ftarget (female-only family size)
  .sample_structure_given_F <- function(Ftarget) {
    # Step 1: choose s in {2,3,4,5}
    s_vals <- 2:5
    w_s <- numeric(length(s_vals))
    for (idx in seq_along(s_vals)) {
      s <- s_vals[idx]
      Nf <- Ftarget - 1L - s  # number of female children required
      if (Nf < 0L) { w_s[idx] <- 0; next }
      pcs <- pmf_children_sum(s)
      K <- pcs$support
      pK <- pcs$pmf[as.character(K)]
      # weight ∝ p_sec(s) * Σ_K pK * choose(K, Nf) * 0.5^K  (only K >= Nf)
      ok <- (K >= Nf)
      if (!any(ok)) { w_s[idx] <- 0; next }
      w_s[idx] <- as.numeric(p_sec[as.character(s)]) *
        sum( pK[ok] * choose(K[ok], Nf) * (0.5^K[ok]) )
    }
    if (sum(w_s) <= 0) stop("Ftarget=", Ftarget, " is infeasible under the original structure.")
    w_s <- w_s / sum(w_s)
    s <- sample(s_vals, 1L, prob = w_s)
    
    # Step 2: choose total children K given s and Nf
    Nf <- Ftarget - 1L - s
    pcs <- pmf_children_sum(s)
    K   <- pcs$support
    pK  <- pcs$pmf[as.character(K)]
    ok <- (K >= Nf)
    wK <- numeric(length(K)); wK[ok] <- pK[ok] * choose(K[ok], Nf) * (0.5^K[ok])
    wK <- wK / sum(wK)
    K_total <- sample(K, 1L, prob = wK)
    
    # Step 3: sample composition (ThiNum_1..ThiNum_s) given K_total
    # build DP for sums: pmf_sum[j][r] = prob sum of j couples equals r
    pmf_sum <- vector("list", s)
    support <- vector("list", s)
    pmf_sum[[1]] <- as.numeric(p_kid); names(pmf_sum[[1]]) <- as.character(2:5)
    support[[1]] <- 2:5
    if (s >= 2) {
      for (j in 2:s) {
        prev <- pmf_sum[[j-1]]; prev_supp <- support[[j-1]]
        cur <- numeric(length(prev) + length(p_kid) - 1L)
        names(cur) <- as.character((min(prev_supp)+2):(max(prev_supp)+5))
        # convolve prev with p_kid
        for (r in prev_supp) for (k in 2:5) {
          idx <- as.character(as.integer(r)+k)
          cur[idx] <- cur[idx] + prev[as.character(r)] * p_kid[as.character(k)]
        }
        pmf_sum[[j]] <- cur
        support[[j]] <- as.integer(names(cur))
      }
    }
    
    # backward sampling of ThiVec
    ThiVec <- integer(s)
    r <- K_total
    for (j in s:2) {
      probs <- numeric(4); names(probs) <- as.character(2:5)
      for (k in 2:5) {
        rem <- r - k
        probs[as.character(k)] <-
          if (rem %in% support[[j-1]]) p_kid[as.character(k)] * pmf_sum[[j-1]][as.character(rem)] else 0
      }
      if (sum(probs) <= 0) stop("Composition sampling failed; no mass for j=", j, " r=", r)
      probs <- probs / sum(probs)
      kj <- as.integer(sample(2:5, 1L, prob = probs))
      ThiVec[j] <- kj
      r <- r - kj
    }
    ThiVec[1] <- r  # remaining
    
    # Step 4: sample per-couple female-child counts via sequential hypergeometric
    m_fem <- integer(s)
    Nf_left <- Nf
    K_left  <- K_total
    for (j in 1:(s-1)) {
      Kj <- ThiVec[j]
      # rhyper draws #whites among k draws:
      # whites=Kj, blacks=(K_left - Kj), draws=Nf_left
      xj <- if (Nf_left > 0) rhyper(1, m = Kj, n = K_left - Kj, k = Nf_left) else 0
      xj <- as.integer(xj)
      if (xj > Kj) xj <- Kj
      m_fem[j] <- xj
      Nf_left <- Nf_left - xj
      K_left  <- K_left  - Kj
    }
    m_fem[s] <- Nf_left
    
    # Step 5: construct child genders per couple (0=female, 1=male) in random order
    child_sex <- vector("list", s)
    for (j in 1:s) {
      Kj <- ThiVec[j]; fj <- m_fem[j]
      g <- c(rep(0L, fj), rep(1L, Kj - fj))
      if (Kj > 1L) g <- sample(g, Kj, replace = FALSE)
      child_sex[[j]] <- g
    }
    
    # second-gen index person's sex per couple: Bern(0.5), then force proband female
    tmpgender <- sample(c(1L,0L), s, replace=TRUE, prob=c(0.5,0.5))
    tmpgender[1] <- 0L
    
    list(SecNum = s, ThiVec = ThiVec, tmpgender = tmpgender, child_sex = child_sex)
  }
  
  ## ---------------- ORIGINAL generator driven either by (a) conditional structure or (b) free draws ----------------
  build_family <- function(struct = NULL) {
    tmpdata <- numeric()
    indID <- c(cumind+1, cumind+2)
    motherID <- c(0,0); fatherID <- c(0,0)
    gender <- c(1,0)       # 1=male, 0=female (your original convention)
    proband <- c(0,0)
    generation <- c(1,1)
    relation  <- c(4,4)
    
    # draw or use given SecNum/tmpgender/ThiVec/child_sex
    if (is.null(struct)) {
      SecNum <- as.integer(sample(2:5, 1, replace=TRUE, prob=as.numeric(p_sec)))
      tmpgender <- sample(c(1L,0L), SecNum, replace=TRUE, prob=c(0.5,0.5))
      tmpgender[1] <- 0L
      ThiVec <- vapply(seq_len(SecNum), function(.) as.integer(sample(2:5,1,prob=as.numeric(p_kid))), 1L)
      child_sex <- lapply(ThiVec, function(n) if (n>0L) sample(c(1L,0L), n, replace=TRUE, prob=c(0.5,0.5)) else integer(0))
    } else {
      SecNum    <- struct$SecNum
      tmpgender <- struct$tmpgender
      ThiVec    <- struct$ThiVec
      child_sex <- struct$child_sex
    }
    
    if (SecNum > 0) {
      NumMem <- 2*SecNum + cumind + 2
      for (j in 1:SecNum) {
        if (j==1) {
          proband  <- c(proband, c(1,0))
          relation <- c(relation, c(1,6))
          tmpgender[j] <- 0L
        } else {
          proband  <- c(proband, c(0,0))
          relation <- c(relation, c(2,7))
        }
        indID    <- c(indID, 2*j+cumind+1, 2*j+cumind+2)
        fatherID <- c(fatherID, c(cumind+1, 0))
        motherID <- c(motherID, c(cumind+2, 0))
        if (tmpgender[j]==0L) gender <- c(gender, c(0,1)) else gender <- c(gender, c(1,0))
        generation <- c(generation, c(2,0))
        
        ThiNum <- ThiVec[j]
        if (ThiNum > 0L) {
          for (k in 1:ThiNum) {
            proband <- c(proband, 0)
            indID   <- c(indID, NumMem + k)
            if (j==1) relation <- c(relation, 3) else relation <- c(relation, 5)
            if (gender[indID==(2*j+cumind+1)]==0L) {
              fatherID <- c(fatherID, 2*j+cumind+2)
              motherID <- c(motherID, 2*j+cumind+1)
            } else {
              fatherID <- c(fatherID, 2*j+cumind+1)
              motherID <- c(motherID, 2*j+cumind+2)
            }
            # use pre-specified or free child sex
            child_g <- if (is.null(struct)) sample(c(1L,0L), 1, replace=TRUE, prob=c(0.5,0.5)) else child_sex[[j]][k]
            gender     <- c(gender, child_g)
            generation <- c(generation, 3)
          }
        }
        NumMem <- NumMem + ThiNum
      } # j
      
      nIndi <- length(indID)
      famID <- rep(i, nIndi)
      ageonset <- rep(0, nIndi); censorage <- rep(0, nIndi); status <- rep(0, nIndi)
      affected <- rep(0, nIndi); disgene.m <- rep(0, nIndi); disgene.s <- rep(0, nIndi)
      ParentsG.m <- rep(0, nIndi)
      pos <- 1:nIndi
      
      ## --- frailty / kinship (unchanged) ---
      if (is.null(frailty.dist)) alpha <- 0
      else if (frailty.dist=="lognormal") alpha <- rnorm(1, mean=0, sd=1/sqrt(depend))
      else if (frailty.dist=="gamma")     alpha <- log(rgamma(1, shape=depend, scale=1/depend))
      else stop("unrecognized frailty distribution")
      
      if (variation == 'kinship'){
        sex <- ifelse(gender == 1, 'male', 'female')
        mped <- pedigree(indID, fatherID, motherID, sex = sex, famid = famID)
        kmat <- as.matrix(kinship(mped))
        alpha <- MASS::mvrnorm(1, mu = rep(0,nrow(kmat)), Sigma = 2*kmat/depend)
      }
      
      ## --- ages (unchanged) ---
      prob.age <- rtruncnorm(1, a=agemin, b=agemax, mean=probandage[1], sd=probandage[2])
      genepos <- pos[generation==1]
      censorage[genepos] <- rnorm(length(genepos), mean=prob.age+20, sd=1.5)
      min.page1 <- min(censorage[genepos])
      genepos <- pos[generation==0]
      if (length(genepos) > 0L) censorage[genepos] <- rnorm(length(genepos), mean=probandage[1], sd=probandage[2])
      genepos <- pos[generation==2]
      for (jj in 1:length(genepos)) {
        if (jj==1) censorage[genepos[1]] <- prob.age
        else censorage[genepos[jj]] <- rtruncnorm(1, a=agemin, b=min.page1-14, mean=prob.age[1]-jj-1, sd=1.5)
        sonpos <- pos[fatherID==indID[genepos[jj]] | motherID==indID[genepos[jj]]]
        if (length(sonpos) > 0L) {
          min.page2 <- min(censorage[indID==fatherID[sonpos[1]] | indID==motherID[sonpos[1]]])
          for (k in 1:length(sonpos)) censorage[sonpos[k]] <- rnorm(1, mean=min.page2-20-k, sd=1.5)
        }
      }
      
      ## --- extra covariate ---
      if (add.x){
        if (x.dist =="normal") newx <- rnorm(nIndi, mean=x.parms[1], sd=x.parms[2])
        else if (x.dist =="binomial") newx <- rbinom(nIndi, size=x.parms[1], prob=x.parms[2])
        else if (x.dist == "mvnormal") {
          if (variation %in% c("frailty", "kinship")) {
            sex <- ifelse(gender == 1, 'male', 'female')
            mped <- pedigree(indID, fatherID, motherID, sex = sex, famid = famID)
            kmat <- as.matrix(kinship(mped))
            Sigma <- 2*(x.parms[2]^2)*kmat
            newx  <- MASS::mvrnorm(1, mu=rep(x.parms[1], nrow(kmat)), Sigma=Sigma)
          } else {
            stop('x.dist = "mvnormal" requires variation = "frailty" or "kinship".')
          }
        }
      }
      
      ## --- genetics (unchanged) ---
      if (length(allelefreq)==1) allelefreq <- c(allelefreq,0)
      AAq <- allelefreq^2; Aaq <- 2*allelefreq*(1-allelefreq); aaq <- (1-allelefreq)^2
      G <- cbind(AAq, Aaq, aaq)
      
      prob.age <- censorage[proband==1]
      prob.sex <- gender[proband==1]
      if (add.x) prob.x <- newx[proband==1]
      prob.alpha <- if (variation == "kinship") alpha[proband==1] else alpha
      
      if (add.x) pGene <- fgeneZX(base.dist, frailty.dist, depend=depend, affage=prob.age-agemin,
                                  affsex=prob.sex, affx=prob.x, interaction=interaction, variation=variation,
                                  base.parms=base.parms, vbeta=vbeta, alpha=prob.alpha, pg=c(0,0), m.carrier=m.carrier,
                                  dominant.m=dominant.m, aq=allelefreq)
      else        pGene <- fgeneZ(base.dist, frailty.dist, depend=depend, affage=prob.age-agemin,
                                  affsex=prob.sex, interaction=interaction, variation=variation,
                                  base.parms=base.parms, vbeta=vbeta, alpha=prob.alpha, pg=c(0,0), m.carrier=m.carrier,
                                  dominant.m=dominant.m, aq=allelefreq)
      
      if (variation=="secondgene") ngene <- c(1,2) else ngene <- 1
      for (g in ngene) {
        if (m.carrier==1 & g==1) {
          gg <- if (dominant.m) 1:2 else 1
          prob.G <- sample(gg, 1, replace=TRUE, prob=pGene[1, gg])
        } else prob.G <- sample(c(1,2,3), 1, replace=TRUE, prob=pGene[g,])
        G1 <- parents.g(prob.G, g=g, allelefreq=allelefreq)
        nsibs <- sum(proband==0 & generation==2)
        sib.G <- kids.g(nsibs, G1)
        if (g==1) {
          disgene.m[generation==1] <- G1
          disgene.m[proband==1] <- prob.G
          disgene.m[proband==0 & generation==2] <- sib.G
        } else {
          disgene.s[generation==1] <- G1
          disgene.s[proband==1] <- prob.G
          disgene.s[proband==0 & generation==2] <- sib.G
        }
      }
      disgene.m[generation==0] <- sample(c(1,2,3), sum(generation==0), replace=TRUE, prob=G[1,])
      if (variation=="secondgene") {
        disgene.s[generation==0] <- sample(c(1,2,3), sum(generation==0), replace=TRUE, prob=G[2,])
      }
      for (ii in indID[generation==3]) {
        m.g <- disgene.m[indID==motherID[indID==ii]]
        f.g <- disgene.m[indID==fatherID[indID==ii]]
        disgene.m[indID==ii] <- kids.g(1, c(m.g, f.g))
        if (variation=="secondgene") {
          disgene.s[indID==ii] <- kids.g(1, c(disgene.s[indID==motherID[indID==ii]],
                                              disgene.s[indID==fatherID[indID==ii]]))
        }
      }
      
      if (dominant.m) majorgene <- ifelse(disgene.m==3, 0, 1) else majorgene <- ifelse(disgene.m==1, 1, 0)
      if (variation=="secondgene") {
        if (dominant.s) secondgene <- ifelse(disgene.s==3, 0, 1) else secondgene <- ifelse(disgene.s==1, 1, 0)
      } else secondgene <- rep(0, nIndi)
      
      x <- cbind(gender, majorgene)
      if (interaction)               x <- cbind(gender, majorgene, gender*majorgene)
      if (variation == "secondgene") x <- cbind(x, secondgene)
      if (add.x)                     x <- cbind(x, newx)
      xbeta <- c(x %*% vbeta)
      
      if (variation == "frailty" | variation == "kinship") {
        uni <- runif(nIndi, 0, 1)
        ageonset <- apply(cbind(xbeta, alpha, uni), 1, inv2.surv, base.dist=base.dist, parms=base.parms)
      } else {
        genepos <- pos[proband==1]; xvbeta <- xbeta[genepos]
        affage  <- censorage[genepos]; affage.min <- ifelse(affage>agemin, affage-agemin, 0)
        uni <- runif(length(genepos), 0, 1)
        ageonset[genepos] <- apply(cbind(xvbeta, affage.min, uni), 1, inv.survp, base.dist=base.dist, parms=base.parms, alpha=alpha)
        
        genepos <- pos[proband==0]; xvbeta <- xbeta[genepos]
        uni <- runif(length(genepos), 0, 1)
        ageonset[genepos] <- apply(cbind(xvbeta, uni), 1, inv.surv, base.dist=base.dist, parms=base.parms, alpha=alpha)
      }
      
      ageonset <- ageonset + agemin
      currentage <- ifelse(censorage > agemax, agemax, censorage)
      time <- pmin(currentage, ageonset)
      status <- ifelse(currentage >= ageonset, 1, 0)
      
      mgene <- majorgene
      mm <- round((length(mgene)-1) * mrate)
      if (mm > 0L) mgene[is.element(indID, sample(indID[!proband], mm))] <- NA
      
      fsize <- length(famID)
      naff  <- sum(status)
      
      if (add.x) {
        tmpdata <- cbind(famID, indID, gender, motherID, fatherID, proband, generation,
                         majorgene=disgene.m, secondgene=disgene.s, ageonset, currentage,
                         time, status, mgene, newx, relation, fsize, naff)
      } else {
        tmpdata <- cbind(famID, indID, gender, motherID, fatherID, proband, generation,
                         majorgene=disgene.m, secondgene=disgene.s, ageonset, currentage,
                         time, status, mgene, relation, fsize, naff)
      }
      return(tmpdata)
    } # SecNum>0
    return(tmpdata)
  } # build_family
  
  ## ---------------- choose path ----------------
  if (!align_female_only) {
    return(build_family(NULL))  # original behavior
  } else {
    # sanitize table; note max feasible female count is 31 (s=5, all 5 children female)
    st <- .sanitize_fem_table(femsize_vals, femsize_probs, min_fem = 3L, max_fem = 31L)
    Ftarget <- sample(st$vals, 1L, prob = st$probs)
    struct  <- .sample_structure_given_F(Ftarget)
    return(build_family(struct))
  }
}


