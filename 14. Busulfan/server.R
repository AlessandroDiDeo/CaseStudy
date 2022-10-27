library(shiny)
library(ggplot2)
library(data.table)
library(rootSolve)

labDat <- data.table(read.csv("labs.csv",head=T,stringsAsFactors=F))[DRUG=="BUSULFAN",c("VAR","LAB"),with=F]

demoDat <- data.table(read.csv("Main_Population.csv", na.strings = "na"))[,c("BSA","CREAT","BMI") := round(.SD,2), .SDcols=c("BSA","CREAT","BMI")]


popFun <- function(n=20) {
  popDat <- demoDat[][sample(1:length(ID),n)]
  popDat <- cbind(popDat,data.table(OM1=rnorm(n), OM2=rnorm(n)), CHILD = cut(popDat$WT, breaks=c(-Inf, 9, Inf), labels = FALSE))
  popDat[order(ID)]
}

parFun <- function(simPop,input){
  merge(simPop,simPop[,simPar(as.list(.SD),input),by="ID"])
}


simPar <- function(indDat,input){
  with(c(indDat,input,par), {
    
    # Clearance
    CLi = (CL*((1/100)*pCL))*exp(OM1*sqrt(ETA1))*(WT/9)**(0.749)
    # Central volume of distribution
    V1i = (V1*((1/100)*pV1))*exp(OM2*sqrt(ETA2))*(WT/9)**(0.958)

    # Dose
    A1 = (c(DOSE, DOSE2)[CHILD])*WT
    
    C = 1/V1i   # Coefficient
    L = CLi/V1i # Exponent
    R = A1/pD1
    
    nest <- function(time,dno=1){
      DV <- R*(C/L *(1-exp(-L*ifelse(time<=pD1,time,pD1)))*exp(-L*(ifelse(time<=pD1,0,time-pD1))))
      if(dno>=5 | sum(time>II)==0) return(DV)
      DV[time>II] <- DV[time>II]+nest(time[time>II]-II,dno+1)
      DV
    }
    

    list(CL=CLi,V1=V1i,D1=pD1, A1=A1,AUC=A1/CLi,CMAX=R*(C/L*(1-exp(-L*pD1))),TMAX=pD1,CSS=(A1/CLi)/II,TOTAL=0)
  })
}

par <- list(CL=2.81,V1=1, D1=2, ETA1=0.0523, ETA2=0.0359 ,EPS=0.035, pV1=100, pCL=100)

simInd <- function(indDat,input,time= exp(seq(0,log(49),len=24))-1){
  with(c(indDat,input), {
    # Calculate rates
    C = 1/V1
    L = CL/V1

    R <- A1/D1
    
    nest <- function(time,dno){
      DV <- R*(C/L *(1-exp(-L*ifelse(time<=pD1,time,pD1)))*exp(-L*(ifelse(time<=pD1,0,time-pD1))))
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
                   SEX = ggplot(simPop(),aes_string(X,Y,col='c("Female","Male")[factor(SEX)]'))+layer1+labs(y=labY,x=labX,title="Busulfan"),
                   AGE = ggplot(simPop(),aes_string(X,Y,col='c("<45","45+")[1+(AGE>45)]'))+layer1+labs(y=labY,x=labX,title="Busulfan"),
                   WT = ggplot(simPop(),aes_string(X,Y,col='c("<9","9+")[1+(WT>9)]'))+layer1+labs(y=labY,x=labX,title="Busulfan"),
                   pCur <- ggplot(simPop(),aes_string(X,Y))+layer1
    )
    if (Y =="TMIC") pCur <- pCur + coord_cartesian(ylim=c(0,100)) + scale_y_continuous(breaks=seq(0,100,5))
    
    if(input$lm) pCur <- pCur + geom_smooth(method="lm",se=FALSE)
    pCur+theme_bw()+labs(x=labX,y=labY,col="",title="Busulfan")
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
    
    invisible(corPop[,cat("\n**********",GRP,"**********",gsub("‘|’","\'",gsub("get\\(X\\)",paste0(X,paste0(rep(" ",6-nchar(X)),collapse="")),gsub("get\\(Y\\)",Y,capture.output(summary(lm(get(Y)~get(X),data=.SD)))))),sep="\n"),by=GRP])
    
    
  })
  
  output$pkPlot <- renderPlot({
    times <- c(do.call(c,lapply(1:(input$ADDL),function(n)(exp(seq(0,log(input$II/4+1),len=input$tps))-1)*4+(n-1)*input$II)),(exp(seq(log(input$II/4+1),log((input$II*2)/4+1),log(input$II)/input$tps))-1)*4+(input$ADDL-1)*input$II,input$II*(input$ADDL+1))
    simDat <- simPop()[,simInd(as.list(.SD),reactiveValuesToList(input),times),by="ID"]
    
    log <- scale_y_continuous()
    enh <- NULL
    if (input$log) log <- scale_y_log10()
    if (input$enh) enh <- geom_point()
    if (input$sum) return(ggplot(simDat[,list(DV=median(DV),U=quantile(DV,0.95),L=quantile(DV,0.05),U2=quantile(DV,0.75),L2=quantile(DV,0.25)),by=c("TIME")],aes(TIME,DV))+geom_ribbon(aes(ymin=L,ymax=U),alpha=1/4,fill="#d9230f") +
                            geom_ribbon(aes(ymin=L2,ymax=U2),alpha=1/4,fill="#d9230f") + geom_line(aes(TIME,DV),col="red",size=1) +
                            theme_bw() + coord_cartesian(xlim=c(0,input$II*(input$ADDL+1))) + log + enh +
                            labs(x="Time (h)" , y="Concentration (mg/l)", title="Busulfan"))
    return(ggplot(simDat,aes(TIME,DV)) + geom_line(aes(group=ID,col=ID),alpha=1/4) +
             geom_line(data=simDat[,list(DV=median(DV)),by="TIME"],col="red",size=1) + theme_bw() +
             theme(legend.position="none")+coord_cartesian(xlim=c(0,input$II*(input$ADDL+1)))+log+enh+labs(x="Time (h)" , y="Concentration (mg/l)", title="Busulfan"))
  })
  
  output$popSum <- renderDataTable({
    popSum <- popDat()[,2:10,with=F]
    popSum <- data.table(Median=apply(popSum,2,median),
                         Mean=apply(popSum,2,mean),
                         SD=apply(popSum,2,sd),
                         "5th percentile"=apply(popSum,2,quantile,0.05),
                         "95th percentile"=apply(popSum,2,quantile,0.95))
    cbind(Var=names(popDat())[2:10],popSum[,round(.SD,2), .SDcols=1:5])
  },options=list(paging=F,ordering=F,searching=F,info=F))#,include.rownames=F,comment=T)
  output$dlDat <- downloadHandler(
    filename = 'popDat.csv',
    content =        function(file) {
      write.csv(simPop(), file)
    }
  )
})