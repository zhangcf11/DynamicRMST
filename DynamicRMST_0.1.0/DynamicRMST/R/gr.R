
gr.longRMST <-
  function (betas) {
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    Y <- eta.yxT + Ztime.b
    eta.t <- eta.tw + alpha * Y
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    sc2<-numeric(ncol(X))
    for (i in 1:ncol(X)){
      sc2[i]<--sum(alpha*Xtime[,i]/delta^2*p1-alpha*Xtime[,i]/delta^2*c((p.byt*eta.t)  %*% wGH),na.rm=TRUE)
    }
    c(sc1 + sc2)
  }
gr.survRMST <-
  function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    delta<-exp(thetas$log.delta)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.tw <- if (!is.null(W)) as.vector(W %*% gammas) else rep(0, n)
    Y <- eta.yxT + Ztime.b
    eta.t <- eta.tw + alpha * Y
    sc.gammas<-numeric(ncol(W))
    for (i in 1:ncol(W))
    {
      sc.gammas[i]<- - sum(W[,i]/delta^2 * p1 - W[,i]/delta^2*c((p.byt*eta.t)%*% wGH), na.rm = TRUE)
    }
    sc.alpha <- - sum((p.byt*(p1-eta.t)*Y/delta^2)%*% wGH, na.rm = TRUE)
    sc.delta<--delta*(-n/delta+sum((p.byt *((p1-eta.t)^2/delta^3))%*% wGH,na.rm = TRUE))
    c(sc.gammas, sc.alpha,sc.delta)
  }
