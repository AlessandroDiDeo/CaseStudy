library(shiny)
library(ggplot2)
library(data.table)

labDat <- data.table(read.csv("labs.csv",head=T,stringsAsFactors=F))[DRUG=="DOCETAXEL",c("VAR","LAB"),with=F]
demoDat <- data.table(read.csv("Main_Population.csv"))[,c("BSA","CRCL","BMI","LBM") := round(.SD,2), .SDcols=c("BSA","CRCL","BMI","LBM")][AGE<=65 & AGE>=18]

  
popFun <- function(n=20) {
  popDat <- demoDat[][sample(1:length(ID),n)]
  GENi <- rbinom(n,1,0.66)
  popDat <- cbind(popDat,data.table(OM1=rnorm(n),OM2=rnorm(n),OM3=rnorm(n),OM4=rnorm(n),GEN=GENi, C1236T=rbinom(n,1,0.5)*GENi))
  popDat[order(ID)]
}

parFun <- function(simPop,input){
  merge(simPop,simPop[,simPar(as.list(.SD),input),by="ID"])
}


simPar <- function(indDat,input){
  with(c(indDat,input,par), {
    
    # Clearance
    CLi = (CL*((1/100)*pCL))*exp(OM1*sqrt(ETA1))*(0.719**((C1236T)*GEN))*(1.05**(1-C1236T*GEN))
    
    # Central volume of distribution - 1
    V1i = (V1*((1/100)*pV1))*exp(OM2*sqrt(ETA2))
    # Peripheral volume of distribution - 2
    V2i = (V2*((1/100)*pV2))
    # Peripheral volume of distribution - 3
    V3i = (V3*((1/100)*pV3))*exp(OM3*sqrt(ETA3))
    # Inter-comparmental clearance 12
    Q2i = (Q2*((1/100)*1))
    # Inter-comparmental clearance 13
    Q3i = (Q3*((1/100)*1))*exp(OM4*sqrt(ETA4))
    # Dose
    # A1 = DOSE*BSA
    A1 = c(DOSE,DOSE1,DOSE2)[GEN+C1236T+1]*BSA
    
    
    k10 <- CLi/V1i
    k12 <- Q2i/V1i
    k13 <- Q3i/V1i
    k21 <- Q2i/V2i
    k31 <- Q3i/V3i
    
    a0 <- k10*k21*k31
    a1 <- k10*k31+k21*k31+k21*k31+k10*k21+k31*k12
    a2 <- k10+k12+k13+k21+k31
    p <- a1-(a2^2)/3
    q <- (2*(a2^3))/27-(a1*a2/3)+a0
    r1 <- sqrt(-((p^3)/27))
    phi <- acos((-q/(2*r1)))/3
    r2 <- 2*exp(log(r1)/3)
    
    l1 <- -1*(cos(phi)*r2-a2/3)
    l2 <- -1*(cos(phi+(2*pi)/3)*r2+a2/3)
    l3 <- -1*(cos(phi+(4*pi)/3)*r2-a2/3)
    
    C1 <- (((k21-l1)*(k31-l1))/((l1-l2)*(l1-l3)))/V1i
    C2 <- (((k21-l2)*(k31-l2))/((l2-l1)*(l2-l3)))/V1i
    C3 <- (((k21-l3)*(k31-l3))/((l3-l2)*(l3-l1)))/V1i
    
    R <- A1/pD1
    
    list(CL=CLi,V1=V1i,V2=V2i,V3=V3i,Q2=Q2i,Q3=Q3i,D1=pD1,A1=A1,AUC=A1/CLi,CMAX=R*((C1/l1)*(1-exp(-l1*pD1))+(C2/l2)*(1-exp(-l2*pD1))+(C3/l3)*(1-exp(-l3*pD1))),CSS=(A1/CLi)/II*24,TMAX=pD1,TOTAL=0)
  })
}

par <- list(CL=54.2,V1=11,V2=14.3,V3=315,D1=1,Q2=6.57,Q3=14,ETA1=0.13,ETA2=0.09,ETA3=0.31,ETA4=0.01,EPS=0.14, pV2=100, pV3=100, pQ2=100, pQ3=100, pCL=100, pV1=100)

simInd <- function(indDat,input,time= exp(seq(0,log(49),len=24))-1){
  with(c(indDat,input), {
    
    k10 <- CL/V1
    k12 <- Q2/V1
    k13 <- Q3/V1
    k21 <- Q2/V2
    k31 <- Q3/V3
    
    a0 <- k10*k21*k31
    a1 <- k10*k31+k21*k31+k21*k31+k10*k21+k31*k12
    a2 <- k10+k12+k13+k21+k31
    p <- a1-(a2^2)/3
    q <- (2*(a2^3))/27-(a1*a2/3)+a0
    r1 <- sqrt(-((p^3)/27))
    phi <- acos((-q/(2*r1)))/3
    r2 <- 2*exp(log(r1)/3)
    
    l1 <- -1*(cos(phi)*r2-a2/3)
    l2 <- -1*(cos(phi+(2*pi)/3)*r2+a2/3)
    l3 <- -1*(cos(phi+(4*pi)/3)*r2-a2/3)
    
    C1 <- (((k21-l1)*(k31-l1))/((l1-l2)*(l1-l3)))/V1
    C2 <- (((k21-l2)*(k31-l2))/((l2-l1)*(l2-l3)))/V1
    C3 <- (((k21-l3)*(k31-l3))/((l3-l2)*(l3-l1)))/V1
    
    R <- A1/D1
    
    nest <- function(time,dno){
      DV <- R*((C1/l1)*(1-exp(-l1*ifelse(time<=D1,time,D1)))*exp(-l1*ifelse(time<=D1,0,time-D1))+(C2/l2)*(1-exp(-l2*ifelse(time<=D1,time,D1)))*exp(-l2*ifelse(time<=D1,0,time-D1))+(C3/l3)*(1-exp(-l3*ifelse(time<=D1,time,D1)))*exp(-l3*ifelse(time<=D1,0,time-D1)))
      if(dno>=ADDL | sum(time>II*24)==0) return(DV)
      DV[time>II*24] <- DV[time>II*24]+nest(time[time>II*24]-II*24,dno+1)
      DV
    }
    
    list(TIME=time,DV=nest(time,1))
  })
}


shinyServer(function(input, output) {
 
   popDat <- reactive(popFun(input$N))
   simPop <- reactive(parFun(popDat(),reactiveValuesToList(input)))

   output$corPlot <- renderPlot({
     Y <- input$COR1
     X <- input$COR2

    layer1 <- geom_point()
    if(X=="SEX") { layer1 <- geom_boxplot(aes_string(x='c("Female","Male")[factor(SEX)]'))   
    } else if (X=="TOTAL") {
      layer1 <- geom_boxplot(aes(x=factor(""))) 
    } else if (X=="GEN") {
      layer1 <- geom_boxplot(aes_string(x='c("Wild type","C1236T Heterozygous","C1236T Homozygous")[GEN+C1236T+1]'))
    } 
    
    labX <- c(labDat$LAB,X)[match(X,c(labDat$VAR,X))]
    labY <- c(labDat$LAB,Y)[match(Y,c(labDat$VAR,Y))]
    
    pCur <- switch(input$STRAT,
                 SEX = ggplot(simPop(),aes_string(X,Y,col='c("Female","Male")[factor(SEX)]'))+layer1+labs(y=labY,x=labX,title="Docetaxel"),
                 AGE = ggplot(simPop(),aes_string(X,Y,col='c("<65","65+")[1+(AGE>=65)]'))+layer1+labs(y=labY,x=labX,title="Docetaxel"),
                 GEN = ggplot(simPop(),aes_string(X,Y,col='c("Wild type","C1236T Heterozygous","C1236T Homozygous")[GEN+C1236T+1]'))+layer1+labs(y=labY,x=labX,title="Docetaxel"),
                 pCur <- ggplot(simPop(),aes_string(X,Y))+layer1
    )

     if(input$lm) pCur <- pCur + geom_smooth(method="lm",se=FALSE)
     pCur+theme_bw()+labs(y=labY,x=labX,title="Docetaxel",col="")
   })
   
   output$lm_info <- renderPrint({
     Y <- input$COR1
     X <- input$COR2
     
     corPop <- simPop()
     corPop$GRP <- switch(input$STRAT,
                            SEX = c("Female","Male")[factor(corPop$SEX)],
                            AGE =c("<65","65+")[1+(corPop$AGE>=65)],
                          GEN =c("Wild type","C1236T Heterozygous","C1236T Homozygous")[corPop$GEN+corPop$C1236T+1],
                            "All")
     
     invisible(corPop[,cat(paste0("(Pearson) Correlation Coefficient (",GRP,"): ",round(cor(get(Y),get(X)),2)),"\n"),by=GRP])
     
     if (!input$lm) return(invisible(NULL))

     invisible(corPop[,cat("\n**********",GRP,"**********",gsub("‘|’","\'",gsub("get\\(X\\)",paste0(X,paste0(rep(" ",6-nchar(X)),collapse="")),gsub("get\\(Y\\)",Y,capture.output(summary(lm(get(Y)~get(X),data=.SD)))))),sep="\n"),by=GRP])
     

     })
   
  output$pkPlot <- renderPlot({
    times <- c(do.call(c,lapply(1:(input$ADDL),function(n)(exp(seq(0,log(input$II*24/4+1),len=input$tps))-1)*4+(n-1)*input$II*24)),(exp(seq(log(input$II*24/4+1),log((input$II*24*2)/4+1),log(input$II*24)/input$tps))-1)*4+(input$ADDL-1)*input$II*24,input$II*24*(input$ADDL+1))
    simDat <- simPop()[,simInd(as.list(.SD),reactiveValuesToList(input),times),by="ID"]
    
    log <- scale_y_continuous()
    enh <- NULL
    if (input$log) log <- scale_y_log10()
    if (input$enh) enh <- geom_point()
    if (input$sum) return(ggplot(simDat[,list(DV=median(DV),U=quantile(DV,0.95, na.rm=T),L=quantile(DV,0.0, na.rm=T),U2=quantile(DV,0.75, na.rm=T),L2=quantile(DV,0.25, na.rm=T)),by="TIME"],aes(TIME/24,DV))+geom_ribbon(aes(ymin=L,ymax=U),alpha=1/4,fill="#d9230f")+geom_ribbon(aes(ymin=L2,ymax=U2),alpha=1/4,fill="#d9230f")+geom_line(aes(TIME/24,DV),col="red",size=1)+theme_bw()+coord_cartesian(xlim=c(0,input$II*(input$ADDL+1)))+log+enh+labs(x="Time (d)" , y="Concentration (mg/l)", title="Docetaxel"))
    return(ggplot(simDat,aes(TIME/24,DV))+geom_line(aes(group=ID,col=ID),alpha=1/4)+geom_line(data=simDat[,list(DV=median(DV,na.rm=T)),by="TIME"],col="red",size=1)+theme_bw()+theme(legend.position="none")+coord_cartesian(xlim=c(0,input$II*(input$ADDL+1)))+log+enh+labs(x="Time (d)" , y="Concentration (mg/l)", title="Docetaxel"))#+scale_colour_manual(values=c("black", "grey50", "grey30", "grey70", "#d9230f")))
  })
  
    output$popSum <- renderDataTable({
      popSum <- popDat()[,2:11,with=F]
      popSum <- data.table(Median=apply(popSum,2,median),
                           Mean=apply(popSum,2,mean),
                           SD=apply(popSum,2,sd),
                           "5th percentile"=apply(popSum,2,quantile,0.05),
                           "95th percentile"=apply(popSum,2,quantile,0.95))
      cbind(Var=names(popDat())[2:11],popSum[,round(.SD,2), .SDcols=1:5])
    },options=list(paging=F,ordering=F,searching=F,info=F))#,include.rownames=F,comment=T)
  output$dlDat <- downloadHandler(
    filename = 'popDat.csv',
    content =        function(file) {
      write.csv(simPop(), file)
    }
  )
})