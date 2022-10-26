library(shiny)
library(ggplot2)
library(data.table)

labDat <- data.table(read.csv("labs.csv",head=T,stringsAsFactors=F))[DRUG=="GABAPENTIN",c("VAR","LAB"),with=F]
demoDat <- data.table(read.csv("Main_Population.csv"))[,c("BSA","CRCL","BMI","LBM") := round(.SD,2), .SDcols=c("BSA","CRCL","BMI","LBM")][AGE>=18]

  
popFun <- function(n=20) {
 popDat <- demoDat[][sample(1:length(ID),n)]
 popDat <- cbind(popDat,data.table(OM1=rnorm(n),OM2=rnorm(n),OM3=rnorm(n),OM4=rnorm(n),OLD=as.numeric(popDat$AGE>=65)))
 popDat[order(ID)]
}

parFun <- function(simPop,input){
  merge(simPop,simPop[,simPar(as.list(.SD),input),by="ID"])
}


simPar <- function(indDat,input){
  with(c(indDat,input,par), {
    CLi = ((CL1^SEX*CL0^(1-SEX))*(pCL/100))*exp(OM1*sqrt(ETA1))*(CRCL/104)**(0.675)*(WT/79)**(0.228)
    V2i = ((V21^SEX*V20^(1-SEX))*(pV2/100))*exp(OM2*sqrt(ETA2))*(WT/79)**(0.544)*(AGE/30)**(-0.239)
    D1i = pD1*exp(OM3*sqrt(ETA3))#*(pD1/100)
    ALAG1i = ALAG1*exp(OM4*sqrt(ETA4))
    A1 <- DOSE*(OLD)+DOSE2*(1-OLD)
    
    l = CLi/V2i
    C = 1/V2i
    list(CL=CLi,V2=V2i,D1=D1i,ALAG1=ALAG1i,A1=A1,AUC=A1/CLi,CMAX=(A1/D1i)*(C/l)*(1-exp(-l*D1i)),TMAX=ALAG1i+D1i,TOTAL=0)
  })
}


par <- list(CL1=6.75,CL0=5.74,V21=86.3,V20=65.6,ALAG1=0.390,D1=6.86,R=1200/6.86,ETA1=0.04,ETA2=0.1,ETA3=0.04,ETA4=0.1,EPS=0.05, pCL=100, pV2=100)


oneComp <- function(indDat,input,time= exp(seq(0,log(49),len=24))-1){
  with(c(indDat,input), {
    nest <- function(time,dno){
      DV <- (A1/D1)*(C/l)*(1-exp(-l*ifelse(time<D1,time,D1)))*exp(-l*ifelse(time>D1,time-D1,0))
      if(dno>=ADDL | sum(time>II)==0) return(DV)
      DV[time>II] <- DV[time>II]+nest(time[time>II]-II,dno+1)
      DV
    }
    l = CL/V2
    C = 1/V2

        list(TIME=time+ALAG1,TIME1=time,DV=nest(time,1))#*exp(rnorm(1,0,sqrt(EPS))))
    
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
      print(X)
      layer1 <- geom_boxplot(aes(x=factor(""))) 
    } 
    
    labX <- c(labDat$LAB,X)[match(X,c(labDat$VAR,X))]
    labY <- c(labDat$LAB,Y)[match(Y,c(labDat$VAR,Y))]
    
    pCur <- switch(input$STRAT,
                 SEX = ggplot(simPop(),aes_string(X,Y,col='c("Female","Male")[factor(SEX)]'))+layer1+labs(y=labY,x=labX,title="Gabapentin"),
                 AGE = ggplot(simPop(),aes_string(X,Y,col='c("<65","65+")[1+(AGE>=65)]'))+layer1+labs(y=labY,x=labX,title="Gabapentin"),
                 pCur <- ggplot(simPop(),aes_string(X,Y))+layer1
    )

     if(input$lm) pCur <- pCur + geom_smooth(method="lm",se=FALSE)
     pCur+theme_bw()+labs(y=labY,x=labX,title="Gabapentin",col="")
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
    simDat <- simPop()[,oneComp(as.list(.SD),reactiveValuesToList(input),times),by="ID"]
    
    
#     ggplot(data=data.frame(x=4*10**seq(0,log10(40/4),len=20),y=1,x2=4*10**seq(log10(40/4),log10(40*2/4),log10(40)/20)),aes(x,y))+geom_point()+geom_point(aes(x2,y))
    
    log <- scale_y_continuous()
    enh <- NULL
    if (input$log) log <- scale_y_log10()
    if (input$enh) enh <- geom_point()
    if (input$sum) return(ggplot(simDat[,list(TIME=median(TIME1),DV=median(DV),U=quantile(DV,0.95),L=quantile(DV,0.05),U2=quantile(DV,0.75),L2=quantile(DV,0.25)),by="TIME1"],aes(TIME,DV))+geom_ribbon(aes(ymin=L,ymax=U),alpha=1/4,fill="#d9230f")+geom_ribbon(aes(ymin=L2,ymax=U2),alpha=1/4,fill="#d9230f")+geom_line(aes(TIME,DV),col="red",size=1)+theme_bw()+coord_cartesian(xlim=c(0,input$II*(input$ADDL+1)))+log+enh+labs(x="Time (h)" , y="Concentration (mg/l)", title="Gabapentin"))
    return(ggplot(simDat,aes(TIME,DV))+geom_line(aes(group=ID,col=ID),alpha=1/4)+geom_line(data=simDat[,list(TIME=median(TIME),DV=median(DV)),by="TIME1"],col="red",size=1)+theme_bw()+theme(legend.position="none")+coord_cartesian(xlim=c(0,input$II*(input$ADDL+1)))+log+enh+labs(x="Time (h)" , y="Concentration (mg/l)", title="Gabapentin"))#+scale_colour_manual(values=c("black", "grey50", "grey30", "grey70", "#d9230f")))
  })
  
    output$popSum <- renderDataTable({
      popSum <- popDat()[,2:11,with=F]
      popSum <- data.table(Median=apply(popSum,2,median),
                           Mean=apply(popSum,2,mean),
                           SD=apply(popSum,2,sd),
                           "5th percentile"=apply(popSum,2,quantile,0.05),
                           "95th percentile"=apply(popSum,2,quantile,0.95))
      temp <- cbind(Var=names(popDat())[2:11],popSum[,round(.SD,2), .SDcols=1:5])
      rbind(temp,data.frame(Var="Number of Patients above 65",Median=nrow(popDat()[AGE>=65]),Mean="",SD="", "5th percentile"="","95th percentile"="",check.names = FALSE))
      
    },options=list(paging=F,ordering=F,searching=F,info=F))#,include.rownames=F,comment=T)
  output$dlDat <- downloadHandler(
    filename = 'popDat.csv',
    content =        function(file) {
      write.csv(simPop(), file)
    }
  )
})