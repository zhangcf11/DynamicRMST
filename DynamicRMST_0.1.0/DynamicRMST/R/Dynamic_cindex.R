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
#' DC<-Dynamic_cindex(object,data,"time","status",5,7)
Dynamic_cindex<-function(object,newdata,time,status,Tstart,Thoriz=NULL,idVar = "id",simulate=FALSE,M=100)
{

  Thoriz <- Thoriz + 1e-07
  id <- newdata[[idVar]]
  id <- match(id, unique(id))
  TermsT<-object$termsT
  #SurvT <- model.response(model.frame(TermsT, newdata))
  #Time <- SurvT[, 1]
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
  #SurvT <- model.response(model.frame(TermsT, newdata2))
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
