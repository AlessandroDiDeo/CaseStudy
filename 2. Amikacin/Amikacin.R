library(data.table)
library(ggplot2)
rm(list=ls())

demoDat <- data.table(read.csv("C:/Users/mm890045/Documents/UCL/Documents/Pharmacology course/2015-2016/case studies/pk models references/population/Main_Population.csv"))[,c("BSA","CRCL","BMI","LBM") := round(.SD,2), .SDcols=c("BSA","CRCL","BMI","LBM")]


# Drug name: Amikacin
# Usual dose is 20 mg/kg QD
# Model structure - IV 2COMP PK
# Covariates: CRCL on CL and WT on V2
# ETA on all parameters

### functions ###

popFun <- function(n=20) {
  popDat <- demoDat[][sample(1:length(ID),n)]
  popDat <- cbind(popDat,data.table(OM1=rnorm(n),OM2=rnorm(n),OM3=rnorm(n),OM4=rnorm(n)))
  popDat[order(ID)]
}

parFun <- function(simPop,input){
  cbind(simPop,simPop[,simPar(as.list(.SD),input),by="ID"])
}


simPar <- function(indDat,input){
  with(c(indDat,input,par), {
    
    # Clearance
    CLi = (CL*(10**pCL))*exp(OM1*sqrt(ETA1))*(CRCL/82)**(-0.7)
    # Central volume of distribution
    V2i = (V2*(10**pV2))*exp(OM2*sqrt(ETA2))*(WT/78)**(-0.9)
    # Peripheral volume of distribution
    V3i = (V2*(10**pV3))*exp(OM3*sqrt(ETA3))
    # Inter-comparmental clearance
    Qi = (Q*(10**pQ))*exp(OM4*sqrt(ETA4))
    # Dose
    A1 = DOSE*WT

    list(CL=CLi,V2=V2i,V3=V3i,Q=Qi,D1=D1,A1=A1,AUC=A1/CLi,CMAX=(A1/V2i),TMAX=0)
  })
}

par <- list(CL=4.3,V2=15.9,V3=21.4,Q=12.2,D1=0.5,ETA1=0.3,ETA2=0.2,ETA3=0.3,ETA4=0.5,EPS=0.1)

simInd <- function(indDat,input,time= exp(seq(0,log(49),len=24))-1){
  with(c(indDat,input), {
    # Calculate rates
    k20<-CL/V2
    k23<-Q/V2
    k32<-Q/V3
    
    a<-k20+k23+k32
    
    l1<-(a+sqrt((a^2)-(4*k20*k32)))/2
    l2<-(a-sqrt((a^2)-(4*k20*k32)))/2
    C1<-((k32-l1)/(l2-l1))/V2
    C2<-((k32-l2)/(l1-l2))/V2
    R <- A1/D1
    
    nest <- function(time,dno){
      DV <- R*(C1/l1*(1-exp(-l1*ifelse(time<=D1,time,D1)))*exp(-l1*(ifelse(time<=D1,0,time-D1)))+C2/l2*(1-exp(-l2*ifelse(time<=D1,time,D1)))*exp(-l2*(ifelse(time<=D1,0,time-D1))))
      if(dno>=ADDL | sum(time>II)==0) return(DV)
      DV[time>II] <- DV[time>II]+nest(time[time>II]-II,dno+1)
      DV
    }
    
    list(TIME=time,DV=nest(time,1))
  })
}


#### input from web ###
# Dose is in mg/kg
input <- list(N=100,ADDL=3,II=24,DOSE=20,pCL=0,pV2=0,pV3=0,pQ=0,pD1=0,tps=12,COR1="CL",COR2="WT")


### to caculate ##

times <- do.call(c,lapply(1:(input$ADDL+1),function(n)(exp(seq(0,log(input$II/3+1),len=input$tps))-1)*3+(n-1)*input$II))
popDat <- popFun(input$N)
simPop <- parFun(popDat,input)
simDat <- simPop[,simInd(as.list(.SD),input,times),by="ID"]

### example plots sent to web ###

ggplot(simDat,aes(TIME,DV))+geom_line(aes(group=ID,col=ID),alpha=1/4)+geom_line(data=simDat[,list(TIME=median(TIME),DV=median(DV)),by="TIME"],col="red",size=1)+theme(legend.position="none")

ggplot(simPop,aes(WT,CL))+geom_point()+labs(y=input$COR1,x=input$COR2)

