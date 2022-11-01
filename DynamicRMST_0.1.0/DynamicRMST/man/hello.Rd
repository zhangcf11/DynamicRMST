\name{DynamicRMST}
\alias{DynamicRMST}
\title{Dynamic RMST prediction model}
\usage{
RMST_Joint(lmeObjects,survObjects,Time,status,timeVar)
survJM(object,newdata,idVar, simulate = TRUE, survTimes = NULL, last.time = NULL, M = 200, CI.levels = c(0.025, 0.975), scale = 1.6)
Dynamic_cindex(object,newdata,time,status,Tstart,Thoriz=NULL,idVar = "id",Dt=NULL,simulate=FALSE,M=100)
}
\description{
This function called RMST_Joint is about building  dynamic RMST prediction model;
This function called survJM is about individual prediction;
This function called Dynamic_cindex is about cindex calculation
}
\examples{
library(MASS)
library(pseudo)
library(geepack)
library(JM)
library("lattice")
data(pbc2)
data(pbc2.id)
data<-data.frame(id=pbc2$id,time=pbc2$years,status=pbc2$status,
year=pbc2$year,drug=pbc2$drug,sex=pbc2$sex,
age=pbc2$age,serBilir=pbc2$serBilir)
data.id<-data.frame(id=pbc2.id$id,time=pbc2.id$years,
status=pbc2.id$status,drug=pbc2.id $drug,
sex= pbc2.id$sex,age= pbc2.id$age,serBilir=
pbc2.id$serBilir)

data$status<-as.numeric(data$status=="dead")
data$drug<-as.numeric(data$drug=="D-penicil")
data$sex<-as.numeric(data$sex=="female")
data.id$status<-as.numeric(data.id$status=="dead")
data.id$drug<-as.numeric(data.id$drug=="D-penicil")
data.id$sex<-as.numeric(data.id$sex=="female")
#longitudinal submodel
fit_nonlinear3 <- lme(log(serBilir) ~ year+ age + sex + drug,
                      random =~year|id,
                       data = data)
#pseudo-observation estimation
p<-pseudomean(data.id$time,data.id$status)
#survival submodel
fitSURV<-geeglm(p~drug+sex+age,data=data.id,id= id)
#Dynamic RMST prediction model
object<-RMST_Joint(fit_nonlinear3,fitSURV,data.id$time,data.id$status,"year")
newdata <- data[data$id ==5, ]
survRes<-survJM(object=object,newdata=newdata, idVar = "id", simulate = TRUE,survTimes = NULL, last.time = NULL,  M = 200, CI.levels = c(0.025, 0.975), scale = 1.6)
survRes$summaries[[1]]
DC<-Dynamic_cindex(object,data,"time","status",5,7)
}


