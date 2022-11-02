print.RMST_Joint<-function(object){
  if (!inherits(object,"RMST_Joint"))
    stop("Object must be of class 'RMST_Joint")
  cat ("\nCall:\n",paste(deparse(object$call),sep="\n",collapse="\n"),"\n\n",sep="")
  cat("Variance Components:\n")
  D<-object$coefficients$D
  ncz<-nrow(D)
  diag.D<-ncz !=ncol(D)
  sds<-if(diag.D) sqrt(D) else sqrt(diag(D))
  if(ncz>1){
    if(diag.D){
      dat<-round(c(D),4)
      names(dat)<-rownames(D)
    }
    else{
      corrs<-cov2cor(D)
      corrs[upper.tri(corrs,TRUE)]<-0
      mat<-round(cbind(sds,corrs[,-ncz]),4)
      mat <- apply(mat, 2, sprintf, fmt = "% .4f")
      mat[mat == mat[1, 2]] <- ""
      mat[1, -1] <- abbreviate(colnames(mat)[-1], 6)
      colnames(mat) <- c(colnames(mat)[1], rep("", ncz - 1))
      dat <- data.frame(mat, check.rows = FALSE, check.names = FALSE)
      names(dat) <- c("StdDev", "Corr", if (ncz > 2) rep(" ", ncz - 2) else NULL)
      row.names(dat) <- dimnames(D)[[1]]
    }
  }else{
    dat <- data.frame("StdDev" = sds, row.names = rownames(D),
                      check.rows = FALSE, check.names = FALSE)
  }
  lis <- list("Random-Effects" = dat, "Residual" = c("Longitudinal" = object$coefficients$sigma,
                                                     "Pseudo-observation RMST" = object$coefficients$delta))
  print(lapply(lis, function (x) if (!is.numeric(x)) x else round(x,  4)))
  cat("\nCoefficients:\n")
  gammas <- c(object$coefficients$gammas,"Assoct" = as.vector(object$coefficients$alpha))
  print(lapply(list("Longitudinal Process" = object$coefficients$betas, "Pseudo-observation RMST Process" = gammas),
               round, digits = 4))
  cat("\nlog-Lik:", object$coefficients$logLik)
  cat("\n\n")
  invisible(object)
}
