# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Title
#'
#' @param lmeObjects Longitudinal submodel
#' @param survObjects Survival submodel
#' @param Time The survival time of subject
#' @param status The status of subject
#' @param timeVar The observation time of longitudinal time-dependent covarites
#'
#' @return \item{coefficients} {a list with the estimated coefficients}
#' @export
#'
#' @examples
#' library(MASS)
#' library(pseudo)
#' library(geepack)
#' library(JM)
#' library("lattice")
#' data(pbc2)
#' data(pbc2.id)
#' head(pbc2)
#' head(pbc2.id)

#'data<-data.frame(id=pbc2$id,time=pbc2$years,status=pbc2$status,
#'                year=pbc2$year,drug=pbc2$drug,sex=pbc2$sex,
#'                age=pbc2$age,serBilir=pbc2$serBilir)
#'data.id<-data.frame(id=pbc2.id$id,time=pbc2.id$years,
#'                    status=pbc2.id$status,drug=pbc2.id $drug,
#'                    sex= pbc2.id$sex,age= pbc2.id$age,serBilir=
#'                     pbc2.id$serBilir)

#'data$status<-as.numeric(data$status=="dead")
#'data$drug<-as.numeric(data$drug=="D-penicil")
#'data$sex<-as.numeric(data$sex=="female")
#'data.id$status<-as.numeric(data.id$status=="dead")
#'data.id$drug<-as.numeric(data.id$drug=="D-penicil")
#'data.id$sex<-as.numeric(data.id$sex=="female")

#' fit_nonlinear3 <- lme(log(serBilir) ~ year+ age + sex + drug,
#'                      random =~year|id,
#'                       data = data)
#' p<-pseudomean(data.id$time,data.id$status)
#' fitSURV<-geeglm(p~drug+sex+age,data=data.id,id= id)
#' object<-RMST_Joint(fit_nonlinear3,fitSURV,data.id$time,data.id$status,"year")
RMST_Joint<-function(lmeObjects,survObjects,Time,status,timeVar)
{
  cl <- match.call()
  idT<-lmeObjects$groups[[1]]
  id<-match(idT,unique(idT))
  b<-data.matrix(ranef(lmeObjects))
  D<-lapply(pdMatrix(lmeObjects$modelStruct$reStruct),"*",lmeObjects$sigma^2)[[1]]
  TermsX<-lmeObjects$terms
  data <- lmeObjects$data[all.vars(TermsX)]#
  data <- data[complete.cases(data), ]
  mfX<-model.frame(TermsX,data=data)
  y<-model.response(mfX,"numeric")
  formYx<-formula(lmeObjects)
  X<-model.matrix(formYx,mfX)#fixed effect
  formYz<-formula(lmeObjects$modelStruct$reStruct[[1]])
  mfz<-model.frame(terms(formYz),data=data)
  TermsZ<-attr(mfz,"terms")
  Z<-model.matrix(formYz,mfz)#random effect
  data.id <- data[!duplicated(id), ]#
  data.id[[timeVar]]<-Time
  Xtime<-model.matrix(formula(lmeObjects),model.frame(TermsX,data=data.id))
  Ztime<-model.matrix(formYz,model.frame(terms(formYz),data=data.id))
  long<-c(X %*% fixef(lmeObjects))+rowSums(Z*b[id,])
  long.id<-tapply(long, id, tail,1)
  ######Pseudo observtion RMST regreesion
  W<-survObjects$geese$X
  WW<-W[,-1]
  p1<-survObjects$y
  fit<-geeglm(p1~WW+long.id,id=c(1:length(p1)))
  gammas<-head(fit$coefficients,-1)
  alpha<-fit$coefficients["long.id"]
  betas<-fixef(lmeObjects)
  sigma<-lmeObjects$sigma
  delta<-sqrt(var(fit$residuals))
  initial_values<-c(list(betas=betas,sigma=sigma,D=D,gammas=gammas,alpha=alpha,delta=delta))
  x<-list(X=X,Z=Z,W=W,idT=idT,Xtime=Xtime,Ztime=Ztime)
  yy<-list(y=y,p1=p1,Time=Time,d=status)
  n<-length(p1)
  N<-length(y)
  ni <- as.vector(tapply(id, id, length))
  XtX <- crossprod(X)
  ZtZ <- lapply(split(Z, id), function (x) crossprod(matrix(x, ncol = ncol(Z))))
  names(ZtZ) <- NULL
  ZtZ <- matrix(unlist(ZtZ), n, ncol(Z) * ncol(Z), TRUE)
  ncz=ncol(Z)
  ##Gauss-Hermite quadrature rule components
  GHk<-if(ncz<3 && nrow(Z)<2000) 15 else 9
  GH <- JM:::gauher(GHk)#11points
  b <- as.matrix(expand.grid(rep(list(GH$x), ncol(Z))))
  k <- nrow(b)
  wGH <- as.matrix(expand.grid(rep(list(GH$w), ncol(Z))))
  wGH <- 2^(ncol(Z)/2) * apply(wGH, 1, prod) * exp(rowSums(b * b))
  diag.D <- !is.matrix(D)
  if (!diag.D) dimnames(D) <- NULL else names(D) <- NULL
  inv_chol_D<-solve(chol(solve(D)))
  det_inv_chol_D<-det(inv_chol_D)
  b <- sqrt(2) * t(inv_chol_D %*% t(b))
  wGH <- wGH * det_inv_chol_D
  dimnames(b) <- NULL
  b2 <- if (ncol(Z) == 1) b * b else t(apply(b, 1, function (x) x %o% x))
  Ztb <- Z %*% t(b)
  Ztime.b <- Ztime %*% t(b)
  environment(opt.survRMST) <- environment(gr.survRMST) <- environment(gr.longRMST) <- environment()
  environment(LogLik.RMSTGH) <- environment(Score.RMSTGH) <- environment()
  ##EM iterations
  iter=50
  Y.mat<-matrix(0,iter+1,ncol(X)+1)
  T.mat<-matrix(0,iter+1,ncol(W)+2)
  B.mat<-matrix(0,iter+1,ncol(Z)*ncol(Z))
  lgLik<-numeric(iter+1)
  conv<-TRUE
  verbose=FALSE
  parscale=NULL
  optimizer="optim"
  iter.qN=350
  only.EM=FALSE
  for(it in 1:iter){
    #save parameter values in matrix
    Y.mat[it, ] <- c(betas, sigma)
    T.mat[it, ] <- c(gammas, alpha,delta)
    B.mat[it,] <- D
    # linear predictors
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    Y <- eta.yxT + Ztime.b
    eta.tw <- if (!is.null(W)) as.vector(W %*% gammas) else rep(0, n)
    eta.t <- eta.tw + alpha * Y
    #E-step
    mu.y<-eta.yx+Ztb
    logNorm<-dnorm(y,mu.y,sigma,TRUE)
    log.p.yb<-rowsum(logNorm,id,reorder=FALSE)
    dimnames(log.p.yb)<-NULL
    mu.T<-eta.t
    log.p.Tb<-dnorm(p1,mu.T,delta,TRUE)
    dimnames(log.p.Tb)<-NULL
    log.p.b<-rep(JM:::dmvnorm(b,rep(0,ncol(Z)),D,TRUE),each=n)
    p.ytb<-exp(log.p.yb+log.p.Tb+log.p.b)
    p.yt<-c(p.ytb %*% wGH)
    p.byt<-p.ytb/p.yt
    post.b<-p.byt %*% (b*wGH)
    post.vb<-if (ncol(Z)==1){
      c(p.byt %*% (b2*wGH))-c(post.b*post.b)
    }else{
      p.byt %*%(b2*wGH)-t(apply(post.b,1,function(x) x %o% x))
    }
    ##compute log-likelihood
    log.p.yt<-log(p.yt)
    lgLik[it]<-sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
    # print results if verbose
    if (verbose) {
      cat("\n\niter:", it, "\n")
      cat("log-likelihood:", lgLik[it], "\n")
      cat("betas:", round(betas, 4), "\n")
      cat("sigma:", round(sigma, 4), "\n")
      if (!is.null(WW))
        cat("gammas:", round(gammas, 4), "\n")
      cat("alpha:", round(alpha, 4), "\n")
      cat("D:", if (!diag.D) round(D[lower.tri(D, TRUE)], 4) else round(D, 4), "\n")
    }
    # check convergence
    if (it > 5) {
      if (lgLik[it] > lgLik[it - 1]) {
        thets1 <- c(Y.mat[it - 1, ], T.mat[it - 1, ], B.mat[it - 1, ])
        thets2 <- c(Y.mat[it, ], T.mat[it, ], B.mat[it, ])
        check1 <- max(abs(thets2 - thets1) / (abs(thets1) + 1e-03)) < 1e-04
        check2 <- (lgLik[it] - lgLik[it - 1]) < 1e-10 * (abs(lgLik[it - 1]) + 1e-10)
        if (check1 || check2) {
          conv <- FALSE
          if (verbose)
            cat("\n\nconverged!\n")
          break}
      }
    }
    if (iter == 0) break
    #M-step
    Zb<-rowSums(Z*post.b[id,],na.rm=TRUE)
    mu<-y-eta.yx
    tr.tZZvarb<-sum(ZtZ*post.vb,na.rm=TRUE)
    sigman<-sqrt(c(crossprod(mu,mu-2*Zb)+crossprod(Zb)+tr.tZZvarb)/N)
    Dn<-matrix(colMeans(p.byt %*% (b2*wGH),na.rm=TRUE),ncol(Z),ncol(Z))
    Dn<-if (diag.D) diag(Dn) else 0.5 * (Dn + t(Dn))
    Hbetas<-JM:::nearPD(JM:::fd.vec(betas,gr.longRMST))
    scbetas<-gr.longRMST(betas)
    betasn<-betas-c(solve(Hbetas,scbetas))

    list.thetas <- list(gammas = gammas, alpha = alpha,log.delta=log(delta))
    list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
    thetas <- unlist(as.relistable(list.thetas))
    optz.surv <- optim(thetas, opt.survRMST, gr.survRMST, method = "BFGS",
                       control = list(maxit = if (it < 5) 20 else 5,
                                      parscale = if (it < 5) rep(0.01, length(thetas)) else rep(0.1, length(thetas))))
    thetasn <- relist(optz.surv$par, skeleton = list.thetas)

    # update parameter values

    betas <- betasn
    sigma <- sigman
    D <- Dn
    gammas <- thetasn$gammas
    alpha <- thetasn$alpha
    delta<-exp(thetasn$log.delta)

  }
  list.thetas <- list(betas = betas, log.sigma = log(sigma), gammas = gammas, alpha = alpha, log.delta=log(delta), D = if (diag.D) log(D) else JM:::chol.transf(D))

  list.thetas <- list.thetas[!sapply(list.thetas, is.null)]

  thetas <- unlist(as.relistable(list.thetas))
  lgLik <- - LogLik.RMSTGH(thetas)
  # if not converged, start quasi-Newton iterations
  if (conv && !only.EM) {
    if (is.null(parscale))
      parscale <- rep(0.01, length(thetas))
    if (verbose)
      cat("\n\nquasi-Newton iterations start.\n\n")
    out <- if (optimizer == "optim") {
      optim(thetas, LogLik.RMSTGH, Score.RMSTGH, method = "BFGS",
            control = list(maxit =iter.qN, parscale = parscale,
                           trace = 10 * verbose))
    } else {
      nlminb(thetas, LogLik.RMSTGH,Score.RMSTGH, scale =parscale,
             control = list(iter.max = iter.qN, trace = 1 * verbose))
    }
    if ((conv <- out$convergence) == 0 || - out[[2]] > lgLik) {
      lgLik <- - out[[2]]
      thetas <- relist(out$par, skeleton = list.thetas)
      betas <- thetas$betas
      sigma <- exp(thetas$log.sigma)
      gammas <- thetas$gammas
      alpha <- thetas$alpha
      delta<- exp(thetas$log.delta)
      D <- thetas$D
      D <- if (diag.D) exp(D) else JM:::chol.transf(D)
      it <- it + if (optimizer == "optim") out$counts[1] else out$iterations
      # compute posterior moments for thetas after quasi-Newton
      eta.yx <- as.vector(X %*% betas)
      eta.tw <- as.vector(W %*% gammas)
      Y <- as.vector(Xtime %*% betas) + Ztime.b
      eta.t <- eta.tw + alpha * Y
      mu.y <- eta.yx + Ztb
      logNorm <- dnorm(y, mu.y, sigma, TRUE)
      log.p.yb <- rowsum(logNorm, id)
      mu.T<-eta.t
      log.p.Tb<-dnorm(p1,mu.T,delta,TRUE)
      log.p.b<-rep(JM:::dmvnorm(b,rep(0,ncol(Z)),D,TRUE),each=n)
      p.ytb <- exp(log.p.yb + log.p.Tb + log.p.b)
      p.yt <- c(p.ytb %*% wGH)
      p.byt <- p.ytb / p.yt
      post.b <- p.byt %*% (b * wGH)
      post.vb <- if (ncol(Z) == 1) {
        c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
      } else {
        (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
      }
      Zb <- if (ncol(Z) == 1) post.b[id] else rowSums(Z * post.b[id, ], na.rm = TRUE)
    }
  }
  Score <- Score.RMSTGH(unlist(thetas))
  # calculate Hessian matrix
  if (verbose) cat("\ncalculating Hessian...\n")
  numeriDeriv = "fd"
  Hessian <- if (numeriDeriv == "fd") {
    JM:::fd.vec(unlist(thetas), Score.RMSTGH, eps = 1e-06)
  } else {
    JM:::cd.vec(unlist(thetas), Score.RMSTGH, eps = 1e-06)
  }
  names(betas) <- names(initial_values$betas)
  if (!diag.D) dimnames(D) <- dimnames(initial_values$D) else names(D) <- names(initial_values$D)
  names(gammas) <- colnames(W)
  nams <- c(paste("Y.", c(names(betas), "sigma"), sep = ""), paste("T.", c(names(gammas), "alpha","delta"), sep = ""),
            paste("B.", if (!diag.D) paste("D", seq(1, ncz * (ncz + 1) / 2), sep = "") else names(D), sep = ""))
  dimnames(Hessian) <- list(nams, nams)
  colnames(post.b) <- colnames(Z)
  out<-list(coefficients = list(betas = betas, sigma = sigma, gammas = gammas, alpha = alpha,
                                delta=delta, D = as.matrix(D)), Score = Score, Hessian = Hessian,
            logLik = lgLik, EB = list(post.b = post.b, post.vb = post.vb,Zb = if (iter == 0) rowSums(Z * post.b[id, ], na.rm = TRUE) else Zb,Ztimeb = rowSums(Ztime * post.b)),  iters = it, convergence = conv,n = n, N = N, ni = ni, id = id)
  H <- out$Hessian
  if (any(is.na(H) | !is.finite(H))) {
    warning("infinite or missing values in Hessian at convergence.\n")
  } else {
    ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    if (!all(ev >= -1e-06 * abs(ev[1])))
      warning("Hessian matrix at convergence is not positive definite.\n")
  }
  out$coefficients <- out$coefficients[!sapply(out$coefficients, is.null)]
  out$x <- x
  out$yy <- yy
  out$times <- data[[timeVar]]
  out$data <- data
  out$data.id <- data.id
  out$termsYx <- TermsX
  out$termsYz <- TermsZ
  out$termsT <- survObjects$terms
  out$formYx <- formYx
  out$formYz <- formYz
  out$formT <- survObjects$geese$formula
  out$timeVar <- timeVar
  out$assignY <- attr(lmeObjects$fixDF, "assign")[-1]
  out$assignT <- survObjects$assign
  out$call <- cl
  class(out) <- "jointModel_RMST"
  out}
summary.jointModel_RMST<-function(object){
  VarCov <- vcov(object)
  betas <- object$coefficients$betas
  indY <- seq(1, length(betas))
  sds <- sqrt(diag(VarCov[indY, indY]))
  coefsY <- cbind("Value" = betas, "Std.Err" = sds, "z-value" = betas / sds, "Lower_CI"=betas-1.96*sds,"Upper_CI"=betas+1.96*sds,"p-value" = 2 * pnorm(abs(betas / sds), lower.tail = FALSE))
  gammas <- c(object$coefficients$gammas,
              "Assoct" = as.vector(object$coefficients$alpha),
              "Assoct.s" = as.vector(object$coefficients$Dalpha))
  indT <- grep("T.", colnames(VarCov), fixed = TRUE)
  indT<-head(indT,-1)
  jj <- grep("Assoct[!^\\.s]", names(gammas))
  ii <- setdiff(grep("Assoct", names(gammas)), jj)
  if (length(ii) > 1) {
    nn <- names(object$coefficients$alpha)
    names(gammas)[ii] <- if (length(nn) == 1) "Assoct" else {
      if (nn[1] == "")
        c("Assoct", paste("Assoct", nn[-1], sep = ":"))
      else
        paste("Assoct", nn, sep = ":")
    }
  }
  if (length(jj) > 1) {
    nn <- names(object$coefficients$Dalpha)
    names(gammas)[jj] <- if (length(nn) == 1) "Assoct.s" else {
      if (nn[1] == "")
        c("Assoct.s", paste("Assoct.s", nn[-1], sep = ":"))
      else
        paste("Assoct.s", nn, sep = ":")
    }
  }
  sds <- if (length(indT) > 1) sqrt(diag(VarCov[indT, indT])) else sqrt(VarCov[indT, indT])
  if (!is.null(object$scaleWB))
    sds <- c(sds, NA)
  coefsT <- cbind("Value" = gammas, "Std.Err" = sds, "z-value" = gammas / sds,"Lower_CI"=gammas-1.96*sds,"Upper_CI"=gammas+1.96*sds,"p-value" = 2 * pnorm(abs(gammas / sds), lower.tail = FALSE))
  out1<-list("CoefTable-Long" = coefsY, "CoefTable-Event" = coefsT, "sigma"=object$coefficients$sigma, "delta"=object$coefficients$delta, "D"=object$coefficients$D)
  class(out1) <- "summary.jointModel_RMST"
  out1
}

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
vcov.jointModel_RMST <-
  function (object, ...) {
    out <- try(solve(object$Hessian), silent = TRUE)
    vmat <- if (!inherits(out, "try-error"))
      structure(out, dimnames = dimnames(object$Hessian))
    else
      structure(ginv(object$Hessian), dimnames = dimnames(object$Hessian))
    (vmat + t(vmat)) / 2
  }


#' Title
#'individual prediction
#' @param object the result of dynamic RMST prediction model
#' @param newdata a data frame that contain the longitudinal and  covariate information for which prediction of dynamic RMST value is requried
#' @param idVar the name of the variable in newdata that identififies the different subjects
#' @param simulate logical; if TRUE, a Monte Carlo approach is used to estimate dynamic RMST value.If FALSE, a joint likelihood function is used to estimate dynamic RMST
#' @param survTimes a numeric vector of times for which prediction dynamic RMST value are to be computed
#' @param last.time this refers to the time the subject is known to survive. If last.time=NULL,then this is automatically taken as the last time each subject provided a longitudinal measurement. If last.time is a numeric vector, then it is assumed to contain this last time point for each subject
#' @param M integer denoting how many Monte Carlo samples to use
#' @param CI.levels a numeric vector of length two that specififies which quantiles to use for the calculation of confifidence interval for the predicted probabilities
#' @param scale a numeric scalar that controls the acceptance rate of the Metropolis-Hastings algorithm
#'
#' @return \item{summaries} {a list with elements numeric matrices with numeric summaries of the predicted dynamic RMST for each subject}
#'
#' @export
#'
#' @examples
#' library(MASS)
#' library(pseudo)
#' library(geepack)
#' library(JM)
#' library("lattice")
#' data(pbc2)
#' data(pbc2.id)
#' head(pbc2)
#' head(pbc2.id)

#'data<-data.frame(id=pbc2$id,time=pbc2$years,status=pbc2$status,
#'                year=pbc2$year,drug=pbc2$drug,sex=pbc2$sex,
#'                age=pbc2$age,serBilir=pbc2$serBilir)
#'data.id<-data.frame(id=pbc2.id$id,time=pbc2.id$years,
#'                    status=pbc2.id$status,drug=pbc2.id $drug,
#'                    sex= pbc2.id$sex,age= pbc2.id$age,serBilir=
#'                     pbc2.id$serBilir)

#'data$status<-as.numeric(data$status=="dead")
#'data$drug<-as.numeric(data$drug=="D-penicil")
#'data$sex<-as.numeric(data$sex=="female")
#'data.id$status<-as.numeric(data.id$status=="dead")
#'data.id$drug<-as.numeric(data.id$drug=="D-penicil")
#'data.id$sex<-as.numeric(data.id$sex=="female")

#' fit_nonlinear3 <- lme(log(serBilir) ~ year+ age + sex + drug,
#'                      random =~year|id,
#'                       data = data)
#' p<-pseudomean(data.id$time,data.id$status)
#' fitSURV<-geeglm(p~drug+sex+age,data=data.id,id= id)
#' object<-RMST_Joint(fit_nonlinear3,fitSURV,data.id$time,data.id$status,"year")
#' newdata <- data[data$id ==5, ]
#' survRes<-survJM(object=object,newdata=newdata, idVar = "id", simulate = TRUE,survTimes = NULL, last.time = NULL,  M = 200, CI.levels = c(0.025, 0.975), scale = 1.6)
#' survRes$summaries[[1]]
survJM<-function(object,newdata,idVar, simulate = TRUE, survTimes = NULL, last.time = NULL, M = 200, CI.levels = c(0.025, 0.975), scale = 1.6){
  if (is.null(survTimes) || !is.numeric(survTimes))
    survTimes <- seq(min((object$yy$Time)),
                     max((object$yy$Time)) + 0.1, length.out = 35)
  timeVar <- object$timeVar
  TermsX <- object$termsYx
  TermsZ <- object$termsYz
  LongFormat <- FALSE
  mfX <- model.frame(TermsX, data = newdata)
  mfZ <- model.frame(TermsZ, data = newdata)
  formYx <- reformulate(attr(delete.response(TermsX), "term.labels"))
  formYz <- object$formYz
  na.ind <- as.vector(attr(mfX, "na.action"))
  na.ind <- if (is.null(na.ind)) {
    rep(TRUE, nrow(newdata))
  } else {
    !seq_len(nrow(newdata)) %in% na.ind
  }
  id <- as.numeric(unclass(newdata[[idVar]]))
  id <- id. <- match(id, unique(id))
  id <- id[na.ind]
  y <- model.response(mfX)
  X <- model.matrix(formYx, mfX)
  Z <- model.matrix(formYz, mfZ)[na.ind, , drop = FALSE]
  TermsT <- object$termsT
  data.id <-newdata[tapply(row.names(newdata), id, tail, n = 1L),]
  idT <- data.id[[idVar]]
  idT <- match(idT, unique(idT))
  mfT <- model.frame(delete.response(TermsT), data = data.id)
  formT <- if (!is.null(kk <- attr(TermsT, "specials")$strata)) {
    strt <- eval(attr(TermsT, "variables"), data.id)[[kk]]
    tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
    reformulate(attr(tt, "term.labels"))
  } else if (!is.null(kk <- attr(TermsT, "specials")$cluster)) {
    tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
    reformulate(attr(tt, "term.labels"))
  } else {
    tt <- attr(delete.response(TermsT), "term.labels")
    if (length(tt)) reformulate(tt) else reformulate("1")
  }
  W <- model.matrix(formT, mfT)
  obs.times <- split(newdata[[timeVar]][na.ind], id)
  last.time <- if (is.null(last.time)) {
    tapply(newdata[[timeVar]], id., tail, n = 1)
  } else if (is.character(last.time) && length(last.time) == 1) {
    tapply(newdata[[last.time]], id., tail, n = 1)
  } else if (is.numeric(last.time) && length(last.time) == nrow(data.id)) {
    last.time
  } else {
    stop("\nnot appropriate value for 'last.time' argument.")
  }

  times.to.pred <- lapply(last.time, function (t) survTimes[survTimes > t])
  n <- object$n
  n.tp <- length(last.time)
  ncx <- ncol(X)
  ncz <- ncol(Z)
  ncww <- ncol(W)
  betas <- object$coefficients[['betas']]
  sigma <- object$coefficients[['sigma']]
  D <- object$coefficients[['D']]
  diag.D <- ncol(D) == 1 & nrow(D) > 1
  D <- if (diag.D) diag(c(D)) else D
  gammas <- object$coefficients[['gammas']]
  alpha <- object$coefficients[['alpha']]
  delta<-object$coefficients[['delta']]
  list.thetas <- list(betas = betas, log.sigma = log(sigma), gammas = gammas, alpha = alpha,
                      delta =log( delta),   D = if (diag.D) log(diag(D)) else JM:::chol.transf(D))
  list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
  thetas <- unlist(as.relistable(list.thetas))
  Var.thetas <- vcov(object)
  environment(log.posterior.b) <- environment(S.b) <- environment()


  # calculate the Empirical Bayes estimates and their (scaled) variance
  modes.b <- matrix(0, n.tp, ncz)
  Vars.b <- vector("list", n.tp)
  for (i in seq_len(n.tp)) {
    betas.new <- betas
    sigma.new <- sigma
    D.new <- D
    gammas.new <- gammas
    alpha.new <- alpha
    delta.new <- delta
    ff <- function (b, y, tt, i) -log.posterior.b(b, y, time=tt, ii = i)
    opt <- try(optim(rep(0, ncz), ff, y = y, tt=last.time[i], i = i,
                     method = "BFGS", hessian = TRUE), TRUE)
    if (inherits(opt, "try-error")) {
      gg <- function (b, y, tt, i) JM:::cd(b, ff, y = y, tt=tt, i = i)
      opt <- optim(rep(0, ncz), ff, gg, y = y, tt=last.time[i],
                   i = i, method = "BFGS", hessian = TRUE)
    }
    modes.b[i, ] <- opt$par
    Vars.b[[i]] <- scale * solve(opt$hessian)
  }
  if (!simulate) {
    res <- vector("list", n.tp)
    for (i in seq_len(n.tp)) {
      S.last <- S.b(last.time[i], modes.b[i, ], i)
      S.pred <- numeric(length(times.to.pred[[i]]))
      for (l in seq_along(S.pred))
        S.pred[l] <- S.b(times.to.pred[[i]][l], modes.b[i, ], i)
      res[[i]] <- cbind(times = c(last.time[i],times.to.pred[[i]][l]), predSurv = c(S.last-last.time[i],S.pred-times.to.pred[[i]][l]))
      rownames(res[[i]]) <- seq_along(c(S.last,S.pred))
    }
  } else {
    out <- vector("list", M)
    success.rate <- matrix(FALSE, M, n.tp)
    b.old <- b.new <- modes.b
    if (n.tp == 1)
      dim(b.old) <- dim(b.new) <- c(1, ncz)
    for (m in 1:M) {
      # Step 1: simulate new parameter value
      set.seed(0104+m)
      thetas.new <- mvrnorm(1, thetas, Var.thetas)
      thetas.new <- relist(thetas.new, skeleton = list.thetas)
      betas.new <- thetas.new$betas
      sigma.new <- exp(thetas.new$log.sigma)
      gammas.new <- thetas.new$gammas
      alpha.new <- thetas.new$alpha
      delta.new<-exp(thetas.new$delta)
      D.new <- thetas.new$D
      D.new <- if (diag.D) exp(D.new) else JM:::chol.transf(D.new)
      SS <- vector("list", n.tp)
      for (i in seq_len(n.tp)) {
        # Step 2: simulate new random effects values
        proposed.b <- JM:::rmvt(1, modes.b[i, ], Vars.b[[i]], 4)
        dmvt.old <- JM:::dmvt(b.old[i, ], modes.b[i, ], Vars.b[[i]], 4, TRUE)
        dmvt.proposed <- JM:::dmvt(proposed.b, modes.b[i, ], Vars.b[[i]], 4, TRUE)
        a <- min(exp(log.posterior.b(proposed.b, y, last.time[[i]], ii = i) + dmvt.old -
                       log.posterior.b(b.old[i, ], y, last.time[[i]], ii = i) - dmvt.proposed), 1)
        ind <- runif(1) <= a
        success.rate[m, i] <- ind
        if (!is.na(ind) && ind)
          b.new[i, ] <- proposed.b
        # Step 3: compute dynamic RMST value

        S.last <- S.b(last.time[i], b.new[i, ], i)
        S.pred <- numeric(length(times.to.pred[[i]]))
        for (l in seq_along(S.pred)){
          S.pred[l] <- S.b(times.to.pred[[i]][l], b.new[i, ], i)
          SS[[i]]<-c(S.last-last.time[i],S.pred-times.to.pred[[i]])
        }

      }
      b.old <- b.new
      out[[m]] <- SS
    }
    res <- vector("list", n.tp)
    for (i in seq_len(n.tp)) {
      rr <- sapply(out, "[[", i)
      if (!is.matrix(rr))
        rr <- rbind(rr)
      res[[i]] <- cbind(
        times =c(last.time[i],times.to.pred[[i]]),
        "Mean" = rowMeans(rr, na.rm = TRUE),
        "Median" = apply(rr, 1, median, na.rm = TRUE),
        "Lower" = apply(rr, 1, quantile, probs = CI.levels[1], na.rm = TRUE),
        "Upper" = apply(rr, 1, quantile, probs = CI.levels[2], na.rm = TRUE)
      )
      rownames(res[[i]]) <- as.character(seq_len(NROW(res[[i]])))
    }
  }
  y <- split(y, id)
  newdata. <- do.call(rbind, mapply(function (d, t) {
    d. <- rbind(d, d[nrow(d), ])
    d.[[timeVar]][nrow(d.)] <- t
    d.
  }, split(newdata, id.), last.time, SIMPLIFY = FALSE))
  id. <- as.numeric(unclass(newdata.[[idVar]]))
  id. <- match(id., unique(id.))
  mfX. <- model.frame(delete.response(TermsX), data = newdata.)
  mfZ. <- model.frame(TermsZ, data = newdata.)
  X. <- model.matrix(formYx, mfX.)
  Z. <- model.matrix(formYz, mfZ.)
  fitted.y <- split(c(X. %*% betas) + rowSums(Z. * modes.b[id., , drop = FALSE]), id.)
  names(res) <-names(y) <- names(last.time) <- names(obs.times) <- unique(unclass(newdata[[idVar]]))
  res <- list(summaries = res, survTimes = survTimes, last.time = last.time,
              obs.times = obs.times, y =y,
              fitted.times = split(newdata.[[timeVar]], factor(newdata.[[idVar]])),
              fitted.y = fitted.y, ry = range(object$y$y, na.rm = TRUE))
  if (simulate) {
    res$full.results <- out
    res$success.rate <- success.rate
    rm(list = ".Random.seed", envir = globalenv())
  }
  class(res) <- "survJM"
  res
}

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



#' Title
#'
#' @param object the dynamic RMST prediction model fitted on the training set
#' @param newdata the test set
#' @param time the name of the variable in data that indicates survival time
#' @param status the name of the variable in data that indicates survival status
#' @param Tstart numeric scalar denoting the time point up to which
#' @param Thoriz longitudinal information is to be used to derive predictions
#' @param idVar the name of the variable in newdata that identififies the different subjects
#' @param simulate this parameter is the same as the parameter in the survJM
#' @param M this parameter is the same as the parameter in the survJM function
#'
#' @return \item{cindex} {a numeric scale denoting the model discrimination}
#' @export
#'
#' @examples
#' library(MASS)
#' library(pseudo)
#' library(geepack)
#' library(JM)
#' library("lattice")
#' data(pbc2)
#' data(pbc2.id)
#' head(pbc2)
#' head(pbc2.id)

#'data<-data.frame(id=pbc2$id,time=pbc2$years,status=pbc2$status,
#'                year=pbc2$year,drug=pbc2$drug,sex=pbc2$sex,
#'                age=pbc2$age,serBilir=pbc2$serBilir)
#'data.id<-data.frame(id=pbc2.id$id,time=pbc2.id$years,
#'                    status=pbc2.id$status,drug=pbc2.id $drug,
#'                    sex= pbc2.id$sex,age= pbc2.id$age,serBilir=
#'                     pbc2.id$serBilir)

#'data$status<-as.numeric(data$status=="dead")
#'data$drug<-as.numeric(data$drug=="D-penicil")
#'data$sex<-as.numeric(data$sex=="female")
#'data.id$status<-as.numeric(data.id$status=="dead")
#'data.id$drug<-as.numeric(data.id$drug=="D-penicil")
#'data.id$sex<-as.numeric(data.id$sex=="female")

#' fit_nonlinear3 <- lme(log(serBilir) ~ year+ age + sex + drug,
#'                      random =~year|id,
#'                       data = data)
#' p<-pseudomean(data.id$time,data.id$status)
#' fitSURV<-geeglm(p~drug+sex+age,data=data.id,id= id)
#' object<-RMST_Joint(fit_nonlinear3,fitSURV,data.id$time,data.id$status,"year")
#' DC<-Dynamic_cindex(object,data,"time","status",5,7)
Dynamic_cindex<-function(object,newdata,time,status,Tstart,Thoriz=NULL,idVar = "id",simulate=FALSE,M=100)
{

  Thoriz <- Thoriz + 1e-07
  id <- newdata[[idVar]]
  id <- match(id, unique(id))
  TermsT<-object$termsT
  Time <-newdata[[time]]
  timeVar <- object$timeVar
  ordTime <- order(Time)
  newdata2 <- newdata[ordTime, ]
  newdata2 <- newdata2[Time[ordTime] > Tstart, ]
  newdata2 <- newdata2[newdata2[[timeVar]] <= Tstart, ]
  pi.u.t <- survJM(object, newdata = newdata2, idVar = idVar, survTimes = Thoriz,
                   simulate = simulate, M = M)
  pi.u.t <- sapply(pi.u.t$summaries, "[", 2,2)
  # find comparable subjects
  id <- newdata2[[idVar]]
  Time <- newdata2[!duplicated(id),][[time]]
  event <- newdata2[!duplicated(id),][[status]]
  names(Time) <- names(event) <- as.character(unique(id))
  if (any(dupl <- duplicated(Time))) {
    Time[dupl] <- Time[dupl] + 1e-07
  }
  if (!all(names(pi.u.t) == names(Time)))
    stop("...")
  total<-con<-0
  wh<-which(event==1)
  nt<-length(Time)
  if (max(wh)==nt) {
    for (m in head(wh,-1)){
      for(n in (m+1):nt){
        if(Time[n]>Time[m]){
          total<-total+1
          if(pi.u.t[n]>pi.u.t[m])
            con<-con+1
          if(pi.u.t[m]==pi.u.t[n])
            con<-con+0.5
        }
      }
    }
  }
  else {
    for (m in wh){
      for(n in (m+1):nt){
        if(Time[n]>Time[m]){
          total<-total+1
          if(pi.u.t[n]>pi.u.t[m])
            con<-con+1
          if(pi.u.t[m]==pi.u.t[n])
            con<-con+0.5
        }
      }
    }
  }

  cindex<-con/total
  out <- list(cindex = cindex, Tstart = Tstart, Thoriz = Thoriz, nr = length(unique(id)),
              classObject = class(object), nameObject = deparse(substitute(object)))
  class(out) <- "Dynamic_cindex"
  out

}

