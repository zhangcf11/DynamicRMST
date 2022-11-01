Score.RMSTGH <-
  function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    betas <- thetas$betas
    sigma <- exp(thetas$log.sigma)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    delta<-exp(thetas$log.delta)
    D <- thetas$D
    D <- if (diag.D) exp(D) else JM:::chol.transf(D)
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
    post.b <-  p.byt %*% (b * wGH)
    post.vb <- if (ncol(Z) == 1) {
      c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
    } else {
      (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
    }
    Zb <- if (ncol(Z) == 1) post.b[id] else rowSums(Z * post.b[id, ], na.rm = TRUE)
    mu <- y - eta.yx
    tr.tZZvarb <- sum(ZtZ * post.vb, na.rm = TRUE)
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    sc2<-numeric(ncol(X))
    for (i in 1:ncol(X)){
      sc2[i]<--sum(alpha*Xtime[,i]/delta^2*p1-alpha*Xtime[,i]/delta^2*c((p.byt*eta.t)  %*% wGH),na.rm=TRUE)
    }
    score.y <- c(c(sc1 + sc2),
                 - sigma * (- N / sigma + drop(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / sigma^3))
    sc.gammas<-numeric(ncol(W))
    for (i in 1:ncol(W))
    {
      sc.gammas[i]<- - sum(W[,i]/delta^2 * p1 -  W[,i]/delta^2*c((p.byt*eta.t)%*% wGH), na.rm = TRUE)
    }
    sc.alpha <- - sum((p.byt*(p1-eta.tw-alpha*Y)*Y/delta^2)%*% wGH, na.rm = TRUE)
    sc.delta<--delta*(-n/delta+sum((p.byt *((p1-eta.t)^2/delta^3))%*% wGH,na.rm = TRUE))
    score.t <- c(sc.gammas, sc.alpha,sc.delta)
    score.b <- if (diag.D) {
      svD <- 1 / D
      svD2 <- svD^2
      cS.postVB <- colSums(as.matrix(post.vb), na.rm = TRUE)
      dim(cS.postVB) <- c(ncz, ncz)
      D * 0.5 * (n * svD - diag(cS.postVB) * svD2 - colSums(as.matrix(post.b^2), na.rm = TRUE) * svD2)
    } else {
      svD <- solve(D)
      dD <- JM:::deriv.D(D)
      ndD <- length(dD)
      D1 <- sapply(dD, function (x) sum(svD * x))
      D2 <- t(sapply(dD, function (x) c(svD %*% x %*% svD)))
      cS.postVB <- colSums(as.matrix(post.vb), na.rm = TRUE)
      out <- numeric(ndD)
      for (i in seq_along(dD)) {
        D.mat <- D2[i, ]
        dim(D.mat) <- c(ncz, ncz)
        out[i] <- sum(D2[i, ] * cS.postVB, na.rm = TRUE) + sum((post.b %*% D.mat) * post.b, na.rm = TRUE)
      }
      J <- JM:::jacobian2(attr(D, "L"), ncz)
      drop(0.5 * (n * D1 - out) %*% J)
    }
    c(score.y, score.t, score.b)
  }
