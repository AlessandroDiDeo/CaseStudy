library(shiny)
library(ggplot2)
library(data.table)
library(DT)
library(rootSolve)
library(plyr)

labDat <- data.table(read.csv("labs.csv",head=T,stringsAsFactors=F))[DRUG=="Drug_X",c("VAR","LAB"),with=F]
demoDat <- data.table(read.csv("Main_Population.csv"))[,c("BSA","CRCL","BMI","LBM") := round(.SD,2), .SDcols=c("BSA","CRCL","BMI","LBM")][AGE<=65 & AGE>=18 & SCR<4]

popFun <- function(n=20, MIC) {
  popDat <- demoDat[][sample(1:length(ID),n)]
  popDat <- cbind(popDat,data.table(OM1=rnorm(n),OM2=rnorm(n)))
  ID_MIC <- expand.grid(ID = popDat$ID, MICi = as.numeric(strsplit(MIC,",")[[1]]))
  popDat <- merge(popDat, ID_MIC)
  popDat$ID_MIC <- paste0(popDat$ID,'_',popDat$MICi)
  popDat[order(ID_MIC)]
}

parFun <- function(simPop,input){
  merge(simPop,simPop[,simPar(as.list(.SD),input), by=c("ID_MIC")], by=c("ID_MIC"))
}


simPar <- function(indDat,input){
  with(c(indDat,input,par), {
    
    # Clearance
    CLi = (CL*((1/100)*pCL))*exp(OM1*sqrt(ETA1))*(1 -0.0069 * (AGE - 65) - 0.2690 * (SCR - 1.05))
    # Central volume of distribution
    V2i = (V2*((1/100)*pV2))*exp(OM2*sqrt(ETA2))*(1 - 0.00486 * (AGE - 65)  + 0.0082 * (LBM - 56.62))
    # Absorption rate
    KAi = (KA*((1/100)*pKA))
    # Dose
    A1 = DOSE1
    
    k<-CLi/V2i
    l1<-CLi/V2i
    C1<-1/V2i
    
    nest <- function(time,dno=1){
      DV <- A1*C1*(KA/(KA-l1))*(exp(-l1*time)-exp(-KA*time))
      if(dno>=ADDL | sum(time>II)==0) return(DV)
      DV[time>II] <- DV[time>II]+nest(time[time>II]-II,dno+1)
      DV
    }
    
    TMIC <- function(MIC) {
      root <- uniroot.all(function(x)(nest(x))-MIC,c(0,II))
      if (length(root)==1) { return(II - root[1])
      } else if (length(root)==0) { return(0)
      } else {
        return(diff(root))
      }
    }
    
    list(CL=CLi,V2=V2i,KA=KAi,A1=A1,AUC=A1/CLi,CMAX=A1/V2i, TMIC=(TMIC(MICi)/II)*100, TOTAL=0)
  })
}

par <- list(CL=0.4,V2=7,KA=0.8,ETA1=0.16,ETA2=0.08,EPS=0.23, pCL=100, pV2=100, pKA=100)

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

shinyServer(function(input, output) {
  
  popDat <- reactive(popFun(input$N, input$MIC))
  simPop <- reactive(parFun(popDat(),reactiveValuesToList(input)))
  
  output$PKPD <- renderPlot({
    Y <- input$COR1
    X <- 'Var'
    
    layer1 <- geom_boxplot(aes(x=factor(""))) 
    
    labY <- c(labDat$LAB,Y)[match(Y,c(labDat$VAR,Y))]
    
    if(input$COR1 == 'TMIC') {
      pCur <- ggplot(simPop(),aes(factor(MICi), TMIC, group = factor(MICi))) + geom_boxplot()
      pCur+theme_bw()+labs(y='T>MIC (%)',x='MIC (mg/L)',title="Mycobactesin",col="")
    } else if (input$COR1 == 'PTA') {
      thres <- as.numeric(input$thres)
      PTA <- ddply(simPop(), .(MICi), here(summarise), PTA = sum(TMIC>thres)/length(TMIC)*100)
      pCur <- ggplot(PTA,aes(factor(MICi), PTA, group = 1)) + geom_line() + geom_point() + geom_hline(yintercept = 90, col = 'green', linetype = 2)
      pCur+theme_bw()+labs(y='Probability of taget attainment (%)',x='MIC (mg/L)',title="Mycobactesin",col="")
    }
    
  })
  
  output$AUC <- renderPlot({
    Y <- 'AUC'
    X <- 'Var'
    
    layer1 <- geom_boxplot(aes(x=factor(""))) 
    
    labY <- c(labDat$LAB,Y)[match(Y,c(labDat$VAR,Y))]
    
    pCur <- ggplot(simPop(),aes_string(X,Y))+layer1
    pCur+theme_bw()+labs(y=labY,x='',title="",col="")
    
  })
  
  output$Cmax <- renderPlot({
    Y <- 'CMAX'
    X <- 'Var'
    
    layer1 <- geom_boxplot(aes(x=factor(""))) 
    
    labY <- c(labDat$LAB,Y)[match(Y,c(labDat$VAR,Y))]
    
    pCur <- ggplot(simPop(),aes_string(X,Y))+layer1+geom_hline(yintercept = 7.5, col = 'red', linetype = 2)
    pCur+theme_bw()+labs(y=labY,x='',title="",col="")
    
  })
  
  output$pkPlot <- renderPlot({
    tps <- 24
    times <- seq(0, input$II*input$ADDL, 0.5)
    simDat <- simPop()[,simInd(as.list(.SD),reactiveValuesToList(input),times),by=c("ID_MIC")]
    
    log <- scale_y_continuous()
    if (input$log) log <- scale_y_log10()
    
    return(ggplot(simDat[,list(TIME2=median(TIME),DV=median(DV),U=quantile(DV,0.95),L=quantile(DV,0.05),U2=quantile(DV,0.75),L2=quantile(DV,0.25)),by=c("TIME")],aes(TIME2,DV))+geom_ribbon(aes(ymin=L,ymax=U),alpha=1/4)+geom_ribbon(aes(ymin=L2,ymax=U2),alpha=1/4)+geom_line(aes(TIME,DV),linewidth=1)+geom_hline(yintercept = 7.5, col = 'red', linetype = 2)+theme_bw()+coord_cartesian(xlim=c(0,input$II*(input$ADDL)))+log+labs(x="Time (h)" , y="Concentration (mg/L)", title="Mycobactesin", col="", fill="")+theme(legend.position="bottom"))
  })
  
  output$popSum <- renderDT({
    popSum <- popDat()[,2:11,with=F]
    popSum <- data.table(Median=apply(popSum,2,median),
                         Mean=apply(popSum,2,mean),
                         SD=apply(popSum,2,sd),
                         "5th percentile"=apply(popSum,2,quantile,0.05),
                         "95th percentile"=apply(popSum,2,quantile,0.95))
    cbind(Var=names(popDat())[2:11],popSum[,round(.SD,2), .SDcols=1:5])
  },options=list(paging=F,ordering=F,searching=F,info=F))#,include.rownames=F,comment=T)
})