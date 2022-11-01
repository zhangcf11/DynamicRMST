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
#' @return
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
