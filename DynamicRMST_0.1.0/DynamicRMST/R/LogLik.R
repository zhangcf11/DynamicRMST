LogLik.RMSTGH <-
  function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    betas <- thetas$betas
    sigma <- exp(thetas$log.sigma)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    delta<-exp(thetas$log.delta)
    D <- thetas$D
    D <- if (diag.D) exp(D) else JM:::chol.transf(D)
    # linear predictors
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    Y <- eta.yxT + Ztime.b
    eta.tw <- if (!is.null(W)) as.vector(W %*% gammas) else rep(0, n)
    eta.t <- eta.tw + alpha * Y
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id); dimnames(log.p.yb) <- NULL
    mu.T<-eta.t
    log.p.Tb<-dnorm(p1,mu.T,delta,TRUE)
    log.p.b<-rep(JM:::dmvnorm(b,rep(0,ncol(Z)),D,TRUE),each=n)
    p.ytb <- exp(log.p.yb + log.p.Tb + log.p.b)
    dimnames(p.ytb) <- NULL
    p.yt <- c(p.ytb %*% wGH)
    p.byt <- p.ytb / p.yt
    log.p.yt <- log(p.yt)
    - sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
  }
