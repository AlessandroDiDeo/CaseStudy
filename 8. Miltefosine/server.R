library(shiny)
library(ggplot2)
library(data.table)

labDat <- data.table(read.csv("labs.csv",head=T,stringsAsFactors=F))[DRUG=="MILTEFOSINE",c("VAR","LAB"),with=F]
demoDat <- data.table(read.csv("Main_Population.csv"))[,c("BSA","CRCL","BMI","LBM") := round(.SD,2), .SDcols=c("BSA","CRCL","BMI","LBM")][AGE<65]

  
popFun <- function(n=20) {
  popDat <- demoDat[][sample(1:length(ID),n)]
  popDat <- cbind(transform(popDat, FFM = (WT/(((30.93*SEX)+(35.98*(0**SEX))) * (HT^2) + WT)) * (HT^2) * ((42.92*SEX)+(37.99*(0**SEX)))),data.table(OM1=rnorm(n),OM2=rnorm(n),OM3=rnorm(n),OM4=rnorm(n),OM5=rnorm(n)),CHILD = cut(popDat$AGE, breaks=c(-Inf, 12, 18, Inf), labels = FALSE))
  popDat$CHILD_LAB <- c("Children", "Adolescents", "Adults")[popDat$CHILD]
  popDat[order(ID)]
}

#FFM = (WT/(((30.93*SEX)+(35.98*(0**SEX))) * (HT^2) + WT)) * (HT^2) * ((42.92*SEX)+(37.99*(0**SEX)))

parFun <- function(simPop,input){
  merge(simPop,simPop[,simPar(as.list(.SD),input),by="ID"])
}


simPar <- function(indDat,input){
  with(c(indDat,input,par), {
    
    # Clearance
    CLi = (CL*((1/100)*pCL))*exp(OM1*sqrt(ETA2))*(FFM/53)**(0.75)
    # Central volume of distribution
    V2i = (V2*((1/100)*pV2))*exp(OM2*sqrt(ETA3))*(FFM/53)**(1)
    # Peripheral volume of distribution
    V3i = (V2*(pV3/100))#*(FFM/53)**(1)
    # Inter-comparmental clearance
    Qi = (Q*((1/100)*pQ))#*(FFM/53)**(0.75)
    # Absorption rate
    KAi = (KA*((1/100)*pKA))*exp(OM3*sqrt(ETA1))
    # Lag time
    ALAG1i = 0#(ALAG1*((1/100)*pALAG1))
    
    # Dose
    # A1 = DOSE*1*WT*(1-CHILD)+DOSE2*1*WT*CHILD
    A1 = (c(DOSE1*WT, DOSE2*WT, DOSE3)[CHILD])
    
    II = (c(II1, II2, II3)[CHILD])

    k20<-CLi/V2i
    k23<-Qi/V2i
    k32<-Qi/V3i
    
    a<-k20+k23+k32
    
    l1<-(a+sqrt((a^2)-(4*k20*k32)))/2
    l2<-(a-sqrt((a^2)-(4*k20*k32)))/2
    C1<-((k32-l1)/(l2-l1))/V2i
    C2<-((k32-l2)/(l1-l2))/V2i
    
    TMAX=(optimize(function(time) A1*(C1*KAi/(KAi-l1)*(exp(-l1*time)-exp(-KAi*time))+(C2*KAi/(KAi-l2)*(exp(-l2*time)-exp(-KAi*time)))), interval=c(0, 300), maximum=TRUE)$maximum)
    
    fmax <- function(test) A1*(C1*KA/(KA-l1)*(exp(-l1*test)-exp(-KA*test))+(C2*KA/(KA-l2)*(exp(-l2*test)-exp(-KA*test))))
    CMAX <- sum(sapply(TMAX+(II*(0:50)), fmax))
    
    list(CL=CLi,V2=V2i,V3=V3i,Q=Qi,KA=KAi,ALAG1=ALAG1i,A1=A1,AUC=A1/CLi,AUC24=(A1/CLi)*(24/II),CMAX=CMAX,TMAX=TMAX,TOTAL=0, II=II)
  })
}

# par <- list(CL=12.1,V2=31.3,V3=290,Q=30.7,KA=0.342,ALAG1=0.433,ETA1=0.3,ETA2=0.1,ETA3=0.1,EPS=0.04, pV3=100, pQ=100)
par <- list(KA=0.416, CL=3.99/24, V2=40.1, Q=0.0347/24, V3=1.75, ETA1=(18.2/100)^2, ETA2=(32.1/100)^2, ETA3=(34.1/100)^2, EPS=0.1156, pV3=100, pQ=100, pCL=100, pV2=100, pKA=100)

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
    
    #II per age
  
    nest <- function(time,dno){
      DV <- A1*(C1*KA/(KA-l1)*(exp(-l1*time)-exp(-KA*time))+(C2*KA/(KA-l2)*(exp(-l2*time)-exp(-KA*time))))
      if(dno>=((ADD*24)/II) | sum(time>II)==0) return(DV)
      DV[time>II] <- DV[time>II]+nest(time[time>II]-II,dno+1)
      DV
    }
    
    list(TIME=time+ALAG1,DV=nest(time,1), CHILD_LAB=CHILD_LAB)
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
                 SEX = ggplot(simPop(),aes_string(X,Y,col='c("Female","Male")[factor(SEX)]'))+layer1+labs(y=labY,x=labX,title="Miltefosine"),
                 # AGE = ggplot(simPop(),aes_string(X,Y,col='c("Children","Adults + Adolescents")[1+(AGE>=16)]'))+layer1+labs(y=labY,x=labX,title="Miltefosine"),
                 AGE = ggplot(simPop(),aes_string(X,Y,col='factor(c("Children", "Adolescents", "Adults")[cut(AGE, breaks=c(-Inf, 12, 18, Inf))], levels=c("Children", "Adolescents", "Adults"))'))+layer1+labs(y=labY,x=labX,title="Miltefosine"),
                 
                 pCur <- ggplot(simPop(),aes_string(X,Y))+layer1
    )
    
    
     if(input$lm) pCur <- pCur + geom_smooth(method="lm",se=FALSE)
     pCur+theme_bw()+labs(y=labY,x=labX,title="Miltefosine",col="")
   })
   
   output$lm_info <- renderPrint({
     Y <- input$COR1
     X <- input$COR2
     
     corPop <- simPop()
     corPop$GRP <- switch(input$STRAT,
                            SEX = c("Female","Male")[factor(corPop$SEX)],
                            AGE = c("Children", "Adolescents","Adults")[corPop$CHILD],
                            "All")
     
     invisible(corPop[,cat(paste0("(Pearson) Correlation Coefficient (",GRP,"): ",round(cor(get(Y),get(X)),2)),"\n"),by=GRP])
     
     if (!input$lm) return(invisible(NULL))

     invisible(corPop[,cat("\n**********",GRP,"**********",gsub("‘|’","\'",gsub("get\\(X\\)",paste0(X,paste0(rep(" ",6-nchar(X)),collapse="")),gsub("get\\(Y\\)",Y,capture.output(summary(lm(get(Y)~get(X),data=.SD)))))),sep="\n"),by=GRP])
     

     })
   
  output$pkPlot <- renderPlot({
    # times <- c(do.call(c,lapply(1:(input$ADDL),function(n)(exp(seq(0,log(max(input$II1, input$II2, input$II3)/4+1),len=input$tps))-1)*4+(n-1)*max(input$II1, input$II2, input$II3))),(exp(seq(log(max(input$II1, input$II2, input$II3)/4+1),log((max(input$II1, input$II2, input$II3)*2)/4+1),log(max(input$II1, input$II2, input$II3))/input$tps))-1)*4+(input$ADDL-1)*max(input$II1, input$II2, input$II3),max(input$II1, input$II2, input$II3)*(input$ADDL+1))
    
    # input <- list(II1=12, II2=24, II3=24, ADD=20, tps=12)
    
    times <- seq(0,input$ADD*24+24, len=input$tps*input$ADD)
    
    # times <- c(do.call(c,lapply(1:(input$ADDL),function(n)(exp(seq(0,log(max(input$II1, input$II2, input$II3)/4+1),len=input$tps))-1)*4+(n-1)*max(input$II1, input$II2, input$II3))),(exp(seq(log(max(input$II1, input$II2, input$II3)/4+1),log((max(input$II1, input$II2, input$II3)*2)/4+1),log(max(input$II1, input$II2, input$II3))/input$tps))-1)*4+(input$ADDL-1)*max(input$II1, input$II2, input$II3),max(input$II1, input$II2, input$II3)*(input$ADDL+1))
    
    
    simDat <- simPop()[,simInd(as.list(.SD),reactiveValuesToList(input),times),by="ID"]
    
    log <- scale_y_continuous()
    enh <- NULL
    if (input$log) log <- scale_y_log10()
    if (input$enh) enh <- geom_point()
    if (input$sum) return(ggplot(simDat[,list(DV=median(DV),U=quantile(DV,0.95),L=quantile(DV,0.05),U2=quantile(DV,0.75),L2=quantile(DV,0.25)),by=c("TIME", "CHILD_LAB")],aes(TIME/24,DV))+geom_ribbon(aes(ymin=L,ymax=U, fill=CHILD_LAB),alpha=1/4)+geom_line(aes(TIME/24,DV, col=CHILD_LAB),size=1)+theme_bw()+coord_cartesian(xlim=c(0,input$ADD+1))+log+enh+labs(x="Time (d)" , y="Concentration (ug/ml)", title="Miltefosine", fill="", col="")+theme(legend.position="bottom"))#)+geom_ribbon(aes(ymin=L2,ymax=U2, fill=CHILD_LAB),alpha=1/4)
    return(ggplot(simDat,aes(TIME/24,DV))+geom_line(aes(group=ID),alpha=1/4)+geom_line(data=simDat[,list(DV=median(DV)),by=c("TIME","CHILD_LAB")],aes(col=CHILD_LAB),size=1)+theme_bw()+theme(legend.position="bottom")+coord_cartesian(xlim=c(0,input$ADD+1))+log+enh+labs(x="Time (d)" , y="Concentration (ug/ml)", title="Miltefosine", col=""))
  })
  
    output$popSum <- renderDataTable({
      popSum <- popDat()[,2:11,with=F]
      popSum <- data.table(Median=apply(popSum,2,median),
                           Mean=apply(popSum,2,mean),
                           SD=apply(popSum,2,sd),
                           "5th percentile"=apply(popSum,2,quantile,0.05),
                           "95th percentile"=apply(popSum,2,quantile,0.95))
      temp <- cbind(Var=names(popDat())[2:11],popSum[,round(.SD,2), .SDcols=1:5])
      # rbind(temp,data.frame(Var="Number of Adults + Adolescents",Median=nrow(popDat()[AGE>=16]),Mean="",SD="", "5th percentile"="","95th percentile"="",check.names = FALSE))
      rbind(temp, transform(popDat()[,list(Var=paste0("Number of ", c("Children", "Adolescents", "Adults")[CHILD]), Median=.N),by="CHILD"][,-1,with=F],Mean="",SD="", "5th percentile"="","95th percentile"=""))
      
    },options=list(paging=F,ordering=F,searching=F,info=F))#,include.rownames=F,comment=T)
  output$dlDat <- downloadHandler(
    filename = 'popDat.csv',
    content =        function(file) {
      write.csv(simPop(), file)
    }
  )
})