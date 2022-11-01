log.posterior.b <-
  function (b, y, time, ii) {
    id.i <- id %in% ii
    idT.i <- idT %in% ii
    X.i <- X[id.i, , drop = FALSE]
    Z.i <- Z[id.i, , drop = FALSE]
    data.id[[timeVar]]<- time
    mfX.id<-model.frame(TermsX,data.id)
    mfZ.id<-model.frame(TermsZ,data.id)
    Xtime.i <- model.matrix(formYx, mfX.id)[idT.i, , drop = FALSE]
    Ztime.i <- model.matrix(formYz, mfZ.id)[idT.i, , drop = FALSE]
    mu.y <- as.vector(X.i %*% betas.new) + rowSums(Z.i * rep(b, each = nrow(Z.i)))
    logNorm <- dnorm(y[id.i], mu.y, sigma.new, TRUE)
    log.p.yb <- sum(logNorm)
    log.p.b <- JM:::dmvnorm(b, rep(0, ncol(Z)), D.new, TRUE)
    eta.tw <- if (!is.null(W)) {
      if (!LongFormat)
        as.vector(W[ii, , drop = FALSE] %*% gammas.new)
      else
        as.vector(W[idT.i, , drop = FALSE] %*% gammas.new)
    } else 0
    Y<-as.vector( Xtime.i %*% betas.new+rowSums(Ztime.i*rep(b,each=nrow(Ztime.i))))
    mu.T <- eta.tw+alpha.new*Y
    Y1<-tail(y[id.i],1)
    time1<-eta.tw+alpha.new*Y1
    log.p.Tb<-dnorm(time1,mu.T,delta.new,TRUE)
    log.p.yb + log.p.Tb + log.p.b
  }
