library(shiny)
library(ggplot2)
library(data.table)

labDat <- data.table(read.csv("labs.csv",head=T,stringsAsFactors=F))[DRUG=="RIVAROXABAN",c("VAR","LAB"),with=F]
demoDat <- data.table(read.csv("Main_Population.csv"))[,c("BSA","CRCL","BMI","LBM") := round(.SD,2), .SDcols=c("BSA","CRCL","BMI","LBM")][AGE<=65 & AGE>=18 & SCR<4]

  
popFun <- function(n=20) {
  popDat <- demoDat[][sample(1:length(ID),n)]
  popDat <- cbind(popDat,data.table(OM1=rnorm(n),OM2=rnorm(n),RTV=rbinom(n,1,0.5)))
  popDat$POP_LAB <- c("HIV negative", "HIV positive")[popDat$RTV+1]
  popDat[order(ID)]
}

parFun <- function(simPop,input){
  merge(simPop,simPop[,simPar(as.list(.SD),input),by="ID"])
}


simPar <- function(indDat,input){
  with(c(indDat,input,par), {
    
    # Clearance
    CLi = (CL*((1/100)*pCL))*exp(OM1*sqrt(ETA1))*(1 -0.0069 * (AGE - 65) - 0.2690 * (SCR - 1.05)) * (0.67**RTV) 
    # Central volume of distribution
    V2i = (V2*((1/100)*pV2))*exp(OM2*sqrt(ETA2))*(1 - 0.00486 * (AGE - 65)  + 0.0082 * (LBM - 56.62))
    # Absorption rate
    KAi = (KA*((1/100)*pKA))
    # Dose
    A1 = c(DOSE1,DOSE2)[RTV+1]*1000
    
    k<-CLi/V2i
    l1<-CLi/V2i
    C1<-1/V2i
    list(CL=CLi,V2=V2i,KA=KAi,A1=A1,AUC=A1/CLi,CMAX=(A1/V2i),TMAX=(optimize(function(time) A1*C1*(KAi/(KAi-l1))*(exp(-l1*time)-exp(-KAi*time)), interval=c(0, 300), maximum=TRUE)$maximum),TOTAL=0)
  })
}

par <- list(CL=7.16,V2=68.69,KA=1.23,ETA1=0.16,ETA2=0.08,EPS=0.23, pCL=100, pV2=100, pKA=100)

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
    
    list(TIME=time,DV=nest(time,1), POP_LAB=POP_LAB)
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
     } else if (X=="RTV") {
       layer1 <- geom_boxplot(aes_string(x='c("w/o RTV","w/ RTV")[1+RTV]'))  
     } else if (X=="TOTAL") {
       print(X)
       layer1 <- geom_boxplot(aes(x=factor(""))) 
     } 
    
    labX <- c(labDat$LAB,X)[match(X,c(labDat$VAR,X))]
    labY <- c(labDat$LAB,Y)[match(Y,c(labDat$VAR,Y))]
    
    pCur <- switch(input$STRAT,
                 SEX = ggplot(simPop(),aes_string(X,Y,col='c("Female","Male")[factor(SEX)]'))+layer1+labs(y=labY,x=labX,title="Rivaroxaban"),
                 AGE = ggplot(simPop(),aes_string(X,Y,col='c("<65","65+")[1+(AGE>=65)]'))+layer1+labs(y=labY,x=labX,title="Rivaroxaban"),
                 RTV = ggplot(simPop(),aes_string(X,Y,col='c("w/o RTV","w/ RTV")[1+RTV]'))+layer1+labs(y=labY,x=labX,title="Rivaroxaban"),
                 pCur <- ggplot(simPop(),aes_string(X,Y))+layer1
    )

     if(input$lm) pCur <- pCur + geom_smooth(method="lm",se=FALSE)
     pCur+theme_bw()+labs(y=labY,x=labX,title="Rivaroxaban",col="")
   })
   
   output$lm_info <- renderPrint({
     Y <- input$COR1
     X <- input$COR2
     
     corPop <- simPop()
     corPop$GRP <- switch(input$STRAT,
                          SEX = c("Female","Male")[factor(corPop$SEX)],
                          AGE =c("<65","65+")[1+(corPop$AGE>=65)],
                          RTV =c("w/o RTV","w/ RTV")[1+corPop$RTV],
                          "All")
     
     invisible(corPop[,cat(paste0("(Pearson) Correlation Coefficient (",GRP,"): ",round(cor(get(Y),get(X)),2)),"\n"),by=GRP])
     
     if (!input$lm) return(invisible(NULL))

     invisible(corPop[,cat("\n**********",GRP,"**********",gsub("‘|’","\'",gsub("get\\(X\\)",paste0(X,paste0(rep(" ",6-nchar(X)),collapse="")),gsub("get\\(Y\\)",Y,capture.output(summary(lm(get(Y)~get(X),data=.SD)))))),sep="\n"),by=GRP])
     

     })
   
  output$pkPlot <- renderPlot({
    times <- c(do.call(c,lapply(1:(input$ADDL),function(n)(exp(seq(0,log(input$II/4+1),len=input$tps))-1)*4+(n-1)*input$II)),(exp(seq(log(input$II/4+1),log((input$II*2)/4+1),log(input$II)/input$tps))-1)*4+(input$ADDL-1)*input$II,input$II*(input$ADDL+1))
    simDat <- simPop()[,simInd(as.list(.SD),reactiveValuesToList(input),times),by="ID"]
    print(head(simDat))
    
    log <- scale_y_continuous()
    enh <- NULL
    if (input$log) log <- scale_y_log10()
    if (input$enh) enh <- geom_point()

    if (input$sum) return(ggplot(simDat[,list(TIME2=median(TIME),DV=median(DV),U=quantile(DV,0.95),L=quantile(DV,0.05),U2=quantile(DV,0.75),L2=quantile(DV,0.25)),by=c("TIME", "POP_LAB")],aes(TIME2,DV))+geom_ribbon(aes(ymin=L,ymax=U, fill=POP_LAB),alpha=1/4)+geom_ribbon(aes(ymin=L2,ymax=U2, fill=POP_LAB),alpha=1/4)+geom_line(aes(TIME,DV, col=POP_LAB),size=1)+theme_bw()+coord_cartesian(xlim=c(0,input$II*(input$ADDL+1)))+log+enh+labs(x="Time (h)" , y="Concentration (ug/L)", title="Rivaroxaban", col="", fill="")+theme(legend.position="bottom"))
    return(ggplot(simDat,aes(TIME,DV))+geom_line(aes(group=ID), col="black",alpha=1/4)+geom_line(data=simDat[,list(DV=median(DV)),by=c("TIME", "POP_LAB")],aes(group=POP_LAB, col=POP_LAB),size=1)+theme_bw()+theme(legend.position="bottom")+coord_cartesian(xlim=c(0,input$II*(input$ADDL+1)))+log+enh+labs(x="Time (h)" , y="Concentration (ug/L)",title="Rivaroxaban", col=""))#+scale_colour_manual(values=c("black", "grey50", "grey30", "grey70", "#d9230f")))
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