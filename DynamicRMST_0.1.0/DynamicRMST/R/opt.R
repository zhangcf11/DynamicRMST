opt.survRMST <-
  function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    delta<-exp(thetas$log.delta)
    eta.tw <- if (!is.null(W)) as.vector(W %*% gammas) else 0
    eta.t <-  eta.tw + alpha * Y
    mu.T<-eta.t
    log.p.Tb<-dnorm(p1,mu.T,delta,TRUE)
    p.bytn <- p.byt * log.p.Tb
    -sum(p.bytn %*% wGH, na.rm = TRUE)
  }
