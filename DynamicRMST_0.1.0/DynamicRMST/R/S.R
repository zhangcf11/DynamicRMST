S.b <-
  function (t, b, ii) {
    idT.i <- idT %in% ii
    eta.tw <- if (!is.null(W))
      as.vector(W[ii, , drop = FALSE] %*% gammas.new)
    else 0
    data.id[[timeVar]]<- t
    mfX.id<-model.frame(TermsX,data.id)
    mfZ.id<-model.frame(TermsZ,data.id)
    Xtime.i <- model.matrix(formYx, mfX.id)[idT.i, , drop = FALSE]
    Ztime.i <- model.matrix(formYz, mfZ.id)[idT.i, , drop = FALSE]
    Y<-as.vector(Xtime.i %*% betas.new+rowSums(Ztime.i*rep(b,each=nrow(Ztime.i))))
    mu.T <- eta.tw+alpha.new*Y
    mu.T
  }
