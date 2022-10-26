library(data.table)
library(ggplot2)
rm(list=ls())

demoDat <- data.table(read.csv("C:/Users/mm890045/Documents/UCL/Documents/Pharmacology course/2015-2016/case studies/pk models references/population/Main_Population.csv"))[,c("BSA","CRCL","BMI","LBM") := round(.SD,2), .SDcols=c("BSA","CRCL","BMI","LBM")]


# Drug name: Rivaroxaban
# Usual dose is 20 mg QD
# Model structure - PO 1COMP PK
# Covariates: AGE and SCR on CL, AGE and LBM on V
# ETA on CL and V

### functions ###

popFun <- function(n=20) {
  popDat <- demoDat[][sample(1:length(ID),n)]
  popDat <- cbind(popDat,data.table(OM1=rnorm(n),OM2=rnorm(n)))
  popDat[order(ID)]
}

parFun <- function(simPop,input){
  cbind(simPop,simPop[,simPar(as.list(.SD),input),by="ID"])
}


simPar <- function(indDat,input){
  with(c(indDat,input,par), {
    
    # Clearance
    CLi = (CL*(10**pCL))*exp(OM1*sqrt(ETA1))*(1 -0.0069 * (AGE - 65) - 0.2690 * (SCR - 1.05))
    # Central volume of distribution
    V2i = (V2*(10**pV2))*exp(OM2*sqrt(ETA2))*(1 - 0.00486 * (AGE - 65)  + 0.0082 * (LBM - 56.62))
    # Absorption rate
    KAi = (KA*(10**pKA))
    # Dose
    A1 = DOSE*1000

    list(CL=CLi,V2=V2i,KA=KAi,A1=A1,AUC=A1/CLi,CMAX=(A1/V2i),TMAX=0)
  })
}

par <- list(CL=7.16,V2=68.69,KA=1.23,ETA1=0.16,ETA2=0.08,EPS=0.23)

simInd <- function(indDat,input,time= exp(seq(0,log(49),len=24))-1){
  with(c(indDat,input), {
    
    # Calculate rates
    k<-CL/V2
    l1<-CL/V2
    C1<-1/V2

    nest <- function(time,dno){
      DV <- A1*C1*(KA/(KA-l1))*(exp(-l1*time)-exp(-KA*time))
      if(dno>=ADDL | sum(time>II)==0) return(DV)
      DV[time>II] <- DV[time>II]+nest(time[time>II]-II,dno+1)
      DV
    }
    
    list(TIME=time,DV=nest(time,1))
  })
}


#### input from web ###
# Dose is in mg
input <- list(N=100,ADDL=2,II=24,DOSE=20,pCL=0,pV2=0,pKA=0,tps=12,COR1="CL",COR2="WT")


### to caculate ##

times <- do.call(c,lapply(1:(input$ADDL+1),function(n)(exp(seq(0,log(input$II/3+1),len=input$tps))-1)*3+(n-1)*input$II))
popDat <- popFun(input$N)
simPop <- parFun(popDat,input)
simDat <- simPop[,simInd(as.list(.SD),input,times),by="ID"]


### example plots sent to web ###

ggplot(simDat,aes(TIME,DV))+geom_line(aes(group=ID,col=ID),alpha=1/4)+geom_line(data=simDat[,list(TIME=median(TIME),DV=median(DV)),by="TIME"],col="red",size=1)+theme(legend.position="none")

ggplot(simPop,aes(WT,CL))+geom_point()+labs(y=input$COR1,x=input$COR2)

