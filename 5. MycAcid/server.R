library(shiny)
library(ggplot2)
library(data.table)

labDat <- data.table(read.csv("labs.csv",head=T,stringsAsFactors=F))[DRUG=="MYCACID",c("VAR","LAB"),with=F]
demoDat <- data.table(read.csv("Main_Population.csv"))[,c("BSA","CRCL","BMI","LBM") := round(.SD,2), .SDcols=c("BSA","CRCL","BMI","LBM")][AGE<=65 & AGE>=18]

  
popFun <- function(n=20) { 
  popDat <- demoDat[][sample(1:length(ID),n)]
  popDat <- cbind(popDat,data.table(OM1=rnorm(n),OM2=rnorm(n),OM3=rnorm(n),OM4=rnorm(n),OM5=rnorm(n)),CYC=rbinom(n,1,0.5))
  popDat[order(ID)]
}

parFun <- function(simPop,input){
  merge(simPop,simPop[,simPar(as.list(.SD),input),by="ID"])
}


ss_fun <- function(A1, KA, l1, l2, C1, C2, II, time) {
  (A1*exp((-KA)*time-l1*time-l2*time)*KA*(C2*exp(II*l2+KA*time+l1*time)*KA-C2*exp(II*KA+II*l2+KA*time+l1*time)*KA-C2*exp(II*l1+II*l2+KA*time+l1*time)*KA+C2*exp(II*KA+II*l1+II*l2+KA*time+l1*time)*KA+C1*exp(II*l1+KA*time+l2*time)*KA-C1*exp(II*KA+II*l1+KA*time+l2*time)*KA-C1*exp(II*l1+II*l2+KA*time+l2*time)*KA+C1*exp(II*KA+II*l1+II*l2+KA*time+l2*time)*KA-C1*exp(II*KA+l1*time+l2*time)*KA-C2*exp(II*KA+l1*time+l2*time)*KA+C1*exp(II*KA+II*l1+l1*time+l2*time)*KA+C2*exp(II*KA+II*l1+l1*time+l2*time)*KA+C1*exp(II*KA+II*l2+l1*time+l2*time)*KA+C2*exp(II*KA+II*l2+l1*time+l2*time)*KA-C1*exp(II*KA+II*l1+II*l2+l1*time+l2*time)*KA-C2*exp(II*KA+II*l1+II*l2+l1*time+l2*time)*KA-C2*exp(II*l2+KA*time+l1*time)*l1+C2*exp(II*KA+II*l2+KA*time+l1*time)*l1+C2*exp(II*l1+II*l2+KA*time+l1*time)*l1-C2*exp(II*KA+II*l1+II*l2+KA*time+l1*time)*l1+C2*exp(II*KA+l1*time+l2*time)*l1-C2*exp(II*KA+II*l1+l1*time+l2*time)*l1-C2*exp(II*KA+II*l2+l1*time+l2*time)*l1+C2*exp(II*KA+II*l1+II*l2+l1*time+l2*time)*l1-C1*exp(II*l1+KA*time+l2*time)*l2+C1*exp(II*KA+II*l1+KA*time+l2*time)*l2+C1*exp(II*l1+II*l2+KA*time+l2*time)*l2-C1*exp(II*KA+II*l1+II*l2+KA*time+l2*time)*l2+C1*exp(II*KA+l1*time+l2*time)*l2-C1*exp(II*KA+II*l1+l1*time+l2*time)*l2-C1*exp(II*KA+II*l2+l1*time+l2*time)*l2+C1*exp(II*KA+II*l1+II*l2+l1*time+l2*time)*l2))/((-1+exp(II*KA))*(-1+exp(II*l1))*(-1+exp(II*l2))*(KA-l1)*(KA-l2))
}


simPar <- function(indDat,input){
  with(c(indDat,input,par), {
    
    # Clearance
    CLi = (CL*((1/100)*pCL))*exp(OM1*sqrt(ETA1))*((CRCL/48)**(-0.12))*(1.1**SEX)*((ALB/30)**(-1.07))*(1+(WT*dCSA*CYC/450)**0.31)
    # Central volume of distribution
    V2i = (V2*((1/100)*pV2))*exp(OM2*sqrt(ETA2))*((CRCL/48)**(-0.62))*((ALB/30)**(-1.13))
    # Peripheral volume of distribution
    V3i = (V2*((1/100)*pV3))*exp(OM3*sqrt(ETA3))
    # Inter-comparmental clearance
    Qi = (Q*((1/100)*pQ))*exp(OM4*sqrt(ETA4))
    # Absorption rate
    KAi = (KA*((1/100)*pKA))*exp(OM5*sqrt(ETA5))
    
    k20<-CLi/V2i
    k23<-Qi/V2i
    k32<-Qi/V3i
    
    a<-k20+k23+k32
    
    l1<-(a+sqrt((a^2)-(4*k20*k32)))/2
    l2<-(a-sqrt((a^2)-(4*k20*k32)))/2
    C1<-((k32-l1)/(l2-l1))/V2i
    C2<-((k32-l2)/(l1-l2))/V2i
    
    # Dose
    # A1 = DOSE
    A1 = DOSE*(1-CYC)+DOSE2*CYC
    
    CTROUGH <- ss_fun(A1, KA, l1, l2, C1, C2, II, 0)
    
    TMAX=(optimize(function(time) A1*(C1*KAi/(KAi-l1)*(exp(-l1*time)-exp(-KAi*time))+(C2*KAi/(KAi-l2)*(exp(-l2*time)-exp(-KAi*time)))), interval=c(0, 300), maximum=TRUE)$maximum)
    CMAX= (function(time) A1*(C1*KAi/(KAi-l1)*(exp(-l1*time)-exp(-KAi*time))+(C2*KAi/(KAi-l2)*(exp(-l2*time)-exp(-KAi*time)))))(TMAX)

    list(CL=CLi,V2=V2i,V3=V3i,Q=Qi,KA=KAi,A1=A1,AUC=A1/CLi*(24/II)*0.5,CMAX=CMAX,CTROUGH=CTROUGH,TMAX=TMAX,TOTAL=0,CSA=WT*dCSA*CYC)
  })
}

par <- list(CL=33,V2=91,V3=237,Q=35,KA=4.1,ETA1=0.03,ETA2=0.5,ETA3=0.4,ETA4=0.3,ETA5=0.6,EPS=0.13, pV3=100, pQ=100, pCL=100, pV2=100, pKA=100)

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
  
    nest <- function(time,dno){
      DV <- A1*(C1*KA/(KA-l1)*(exp(-l1*time)-exp(-KA*time))+(C2*KA/(KA-l2)*(exp(-l2*time)-exp(-KA*time))))
      if(dno>=ADDL | sum(time>II)==0) return(DV)
      DV[time>II] <- DV[time>II]+nest(time[time>II]-II,dno+1)
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
    } else if (X=="CSA") {
      layer1 <- geom_boxplot(aes_string(x='c("w/o Cyclosporin","w/ Cyclosporin")[CYC+1]')) 
      # layer1 <- geom_boxplot(aes_string(x="factor(seq(0,20000,100)[cut(CSA,breaks=c(-Inf,seq(0,20000,100),Inf),labels=F)])")) 
      
      } 
    
    labX <- c(labDat$LAB,X)[match(X,c(labDat$VAR,X))]
    labY <- c(labDat$LAB,Y)[match(Y,c(labDat$VAR,Y))]

    pCur <- switch(input$STRAT,
                 SEX = ggplot(simPop(),aes_string(X,Y,col='c("Female","Male")[factor(SEX)]'))+layer1+labs(y=labY,x=labX,title="Mycophenolic acid"),
                 AGE = ggplot(simPop(),aes_string(X,Y,col='c("<65","65+")[1+(AGE>=65)]'))+layer1+labs(y=labY,x=labX,title="Mycophenolic acid"),
                 CSA = ggplot(simPop(),aes_string(X,Y,col='c("w/o Cyclosporin","w/ Cyclosporin")[CYC + 1]'))+layer1+labs(y=labY,x=labX,title="Mycophenolic acid"),
                 # CSA = ggplot(simPop(),aes_string(X,Y,col="factor(seq(0,20000,100)[cut(CSA,breaks=c(-Inf,seq(0,20000,100),Inf),labels=F)])"))+layer1+labs(y=labY,x=labX,title="Mycophenolic acid"),
                 pCur <- ggplot(simPop(),aes_string(X,Y))+layer1
    )

     if(input$lm) pCur <- pCur + geom_smooth(method="lm",se=FALSE)
     pCur+theme_bw()+labs(y=labY,x=labX,title="Mycophenolic acid",col="")
   })
   
   output$lm_info <- renderPrint({
     Y <- input$COR1
     X <- input$COR2
     
     corPop <- simPop()
     corPop$GRP <- switch(input$STRAT,
                            SEX = c("Female","Male")[factor(corPop$SEX)],
                            AGE =c("<65","65+")[1+(corPop$AGE>=65)],
                            "All")
     
     invisible(corPop[,cat(paste0("(Pearson) Correlation Coefficient (",GRP,"): ",round(cor(get(Y),get(X)),2)),"\n"),by=GRP])
     
     if (!input$lm) return(invisible(NULL))

     invisible(corPop[,cat("\n**********",GRP,"**********",gsub("‘|’","\'",gsub("get\\(X\\)",paste0(X,paste0(rep(" ",6-nchar(X)),collapse="")),gsub("get\\(Y\\)",Y,capture.output(summary(lm(get(Y)~get(X),data=.SD)))))),sep="\n"),by=GRP])
     

     })
   
  output$pkPlot <- renderPlot({
    times <- c(do.call(c,lapply(1:(input$ADDL),function(n)(exp(seq(0,log(input$II/4+1),len=input$tps))-1)*4+(n-1)*input$II)),(exp(seq(log(input$II/4+1),log((input$II*2)/4+1),log(input$II)/input$tps))-1)*4+(input$ADDL-1)*input$II,input$II*(input$ADDL+1))
    simDat <- simPop()[,simInd(as.list(.SD),reactiveValuesToList(input),times),by="ID"]
    
    log <- scale_y_continuous()
    enh <- NULL
    if (input$log) log <- scale_y_log10()
    if (input$enh) enh <- geom_point()
    if (input$sum) return(ggplot(simDat[,list(DV=median(DV),U=quantile(DV,0.95),L=quantile(DV,0.05),U2=quantile(DV,0.75),L2=quantile(DV,0.25)),by="TIME"],aes(TIME,DV))+geom_ribbon(aes(ymin=L,ymax=U),alpha=1/4,fill="#d9230f")+geom_ribbon(aes(ymin=L2,ymax=U2),alpha=1/4,fill="#d9230f")+geom_line(aes(TIME,DV),col="red",size=1)+theme_bw()+coord_cartesian(xlim=c(0,input$II*(input$ADDL+1)))+log+enh+labs(x="Time (h)" , y="Concentration (mg/l)", title="Mycophenolic acid"))
    return(ggplot(simDat,aes(TIME,DV))+geom_line(aes(group=ID,col=ID),alpha=1/4)+geom_line(data=simDat[,list(DV=median(DV)),by="TIME"],col="red",size=1)+theme_bw()+theme(legend.position="none")+coord_cartesian(xlim=c(0,input$II*(input$ADDL+1)))+log+enh+labs(x="Time (h)" , y="Concentration (mg/l)", title="Mycophenolic acid"))#+scale_colour_manual(values=c("black", "grey50", "grey30", "grey70", "#d9230f")))
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