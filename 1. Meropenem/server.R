library(shiny)
library(ggplot2)
library(data.table)
library(rootSolve)

labDat <- data.table(read.csv("labs.csv",head=T,stringsAsFactors=F))[DRUG=="MEROPENEM",c("VAR","LAB"),with=F]

demoDat <- data.table(read.csv("Main_Population.csv"))[,c("BSA","CRCL","BMI","LBM") := round(.SD,2), .SDcols=c("BSA","CRCL","BMI","LBM")][AGE<=65 & AGE>=18]

  
popFun <- function(n=20) {
  popDat <- demoDat[][sample(1:length(ID),n)]
  popDat <- cbind(popDat,data.table(OM1=rnorm(n),OM2=rnorm(n),OM3=rnorm(n),OM4=rnorm(n)))
  popDat[order(ID)]
}

parFun <- function(simPop,input){
  merge(simPop,simPop[,simPar(as.list(.SD),input),by="ID"])
}


simPar <- function(indDat,input){
  with(c(indDat,input,par), {
    
    # Clearance
    CLi = (CL*((1/100)*pCL))*exp(OM1*sqrt(ETA1))*(CRCL/83)**(0.62)*(AGE/35)**(-0.34)
    # Central volume of distribution
    V2i = (V2*((1/100)*pV2))*exp(OM2*sqrt(ETA2))*(WT/70)**(0.99)
    # Peripheral volume of distribution
    V3i = (V2*((1/100)*pV3))*exp(OM3*sqrt(ETA3))
    # Inter-comparmental clearance
    Qi = (Q*((1/100)*pQ))*exp(OM4*sqrt(ETA4))
    # Dose
    A1 = DOSE

    k20<-CLi/V2i
    k23<-Qi/V2i
    k32<-Qi/V3i
    
    a<-k20+k23+k32
    
    l1<-(a+sqrt((a^2)-(4*k20*k32)))/2
    l2<-(a-sqrt((a^2)-(4*k20*k32)))/2
    C1<-((k32-l1)/(l2-l1))/V2i
    C2<-((k32-l2)/(l1-l2))/V2i
    R <- A1/pD1

    nest <- function(time,dno=1){
      DV <- R*(C1/l1*(1-exp(-l1*ifelse(time<=pD1,time,pD1)))*exp(-l1*(ifelse(time<=pD1,0,time-pD1)))+C2/l2*(1-exp(-l2*ifelse(time<=pD1,time,pD1)))*exp(-l2*(ifelse(time<=pD1,0,time-pD1))))
      if(dno>=5 | sum(time>II)==0) return(DV)
      DV[time>II] <- DV[time>II]+nest(time[time>II]-II,dno+1)
      DV
    }

    TMIC <- function() {
      root <- uniroot.all(function(x)(nest(x))-16,c(0,II))
      if (length(root)==1) { return(II - root[1])
      } else if (length(root)==0) { return(0)
      } else {
        return(diff(root))
      }
    }
    # (TMIC())
    list(CL=CLi,V2=V2i,V3=V3i,D1=pD1,Q=Qi,A1=A1,AUC=A1/CLi,CMAX=R*(C1/l1*(1-exp(-l1*pD1))+C2/l2*(1-exp(-l2*pD1))),TMAX=pD1,CSS=(A1/CLi)/II,TMIC=(TMIC()/II)*100,TOTAL=0)
  })
}

par <- list(CL=14.6,V2=10.8,V3=12.6,D1=0.5,Q=18.6,ETA1=0.118,ETA2=0.143,ETA3=0.290,ETA4=0.102,EPS=0.035, pV3=100, pQ=100, pCL=100, pV2=100)

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
     } 
    
    labX <- c(labDat$LAB,X)[match(X,c(labDat$VAR,X))]
    labY <- c(labDat$LAB,Y)[match(Y,c(labDat$VAR,Y))]
    pCur <- switch(input$STRAT,
                 SEX = ggplot(simPop(),aes_string(X,Y,col='c("Female","Male")[factor(SEX)]'))+layer1+labs(y=labY,x=labX,title="Meropenem"),
                 AGE = ggplot(simPop(),aes_string(X,Y,col='c("<45","45+")[1+(AGE>45)]'))+layer1+labs(y=labY,x=labX,title="Meropenem"),
                 pCur <- ggplot(simPop(),aes_string(X,Y))+layer1
    )
    if (Y =="TMIC") pCur <- pCur + coord_cartesian(ylim=c(0,100)) + scale_y_continuous(breaks=seq(0,100,5))

     if(input$lm) pCur <- pCur + geom_smooth(method="lm",se=FALSE)
     pCur+theme_bw()+labs(x=labX,y=labY,col="",title="Meropenem")
   })
   
   output$lm_info <- renderPrint({
     Y <- input$COR1
     X <- input$COR2
     

     
     corPop <- simPop()
     corPop$GRP <- switch(input$STRAT,
                            SEX = c("Female","Male")[factor(corPop$SEX)],
                            AGE =c("<45","45+")[1+(corPop$AGE>45)],
                            "All")
     
     invisible(corPop[,cat(paste0("(Pearson) Correlation Coefficient (",GRP,"): ",round(cor(get(Y),get(X)),2)),"\n"),by=GRP])
     
     if (!input$lm) return(invisible(NULL))

     invisible(corPop[,cat("\n**********",GRP,"**********",gsub("???|???","\'",gsub("get\\(X\\)",paste0(X,paste0(rep(" ",6-nchar(X)),collapse="")),gsub("get\\(Y\\)",Y,capture.output(summary(lm(get(Y)~get(X),data=.SD)))))),sep="\n"),by=GRP])
     

     })
   
  output$pkPlot <- renderPlot({
    times <- c(do.call(c,lapply(1:(input$ADDL),function(n)(exp(seq(0,log(input$II/4+1),len=input$tps))-1)*4+(n-1)*input$II)),(exp(seq(log(input$II/4+1),log((input$II*2)/4+1),log(input$II)/input$tps))-1)*4+(input$ADDL-1)*input$II,input$II*(input$ADDL+1))
    simDat <- simPop()[,simInd(as.list(.SD),reactiveValuesToList(input),times),by="ID"]
    
    log <- scale_y_continuous()
    enh <- NULL
    if (input$log) log <- scale_y_log10()
    if (input$enh) enh <- geom_point()
    if (input$sum) return(ggplot(simDat[,list(DV=median(DV),U=quantile(DV,0.95),L=quantile(DV,0.05),U2=quantile(DV,0.75),L2=quantile(DV,0.25)),by="TIME"],aes(TIME,DV))+geom_ribbon(aes(ymin=L,ymax=U),alpha=1/4,fill="#d9230f")+geom_ribbon(aes(ymin=L2,ymax=U2),alpha=1/4,fill="#d9230f")+geom_line(aes(TIME,DV),col="red",size=1)+theme_bw()+coord_cartesian(xlim=c(0,input$II*(input$ADDL+1)))+log+enh+labs(x="Time (h)" , y="Concentration (mg/l)", title="Meropenem")+geom_hline(yintercept=16,linetype=2,col="purple",alpha=1/2)+geom_text(data=data.frame(TIME=input$ADDL*input$II,DV=16),aes(label="MIC",group="hans"),col="purple",alpha=1/2))
    return(ggplot(simDat,aes(TIME,DV))+geom_line(aes(group=ID,col=ID),alpha=1/4)+geom_line(data=simDat[,list(DV=median(DV)),by="TIME"],col="red",size=1)+theme_bw()+theme(legend.position="none")+coord_cartesian(xlim=c(0,input$II*(input$ADDL+1)))+log+enh+labs(x="Time (h)" , y="Concentration (mg/l)", title="Meropenem")+geom_hline(yintercept=16,linetype=2,col="purple",alpha=1/2)+geom_text(data=data.frame(TIME=input$ADDL*input$II,DV=16),aes(label="MIC",group="hans"),col="purple",alpha=1/2))#+scale_colour_manual(values=c("black", "grey50", "grey30", "grey70", "#d9230f")))
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