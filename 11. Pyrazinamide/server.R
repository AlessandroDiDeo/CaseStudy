# function(input, output, session) {
library(data.table)
library(ggplot2)
# Combine the selected variables into a new data frame
par <- list(KA=3.56, CL=3.42, V2=29.2, ETA1=0.0351, ETA2=0.0251, ETA3=0.957, EPS1=1.89, EPS2=0.0907)


times <- seq(0,2*24,by=0.2)
labTab <- list(AUC="Area Under the Curve [mg/L/h]", CMAX="Maximum Concentration [mg/L]", TMAX = "Time of Maximum Concentration [h]", CL = "Clearance [L/h]", WT= "Weight [kg]", V2 = "Volume of Distribution [L]")

simInd <- function(indDat,input,time= exp(seq(0,log(49),len=24))-1){
  with(c(indDat,input), {
    if (is.logical(time)) time <- sort(c(0, runif(nrand, 0, 48)))
    # Calculate rates
    k<-CL/V2
    l1<-CL/V2
    C1<-1/V2
        # nest <- function(time,dno){
    DV <- A1*C1*(KA/(KA-l1))*(exp(-l1*time)-exp(-KA*time))
      # if(dno>=ADDL | sum(time>II)==0) return(DV)
      # DV[time>II] <- DV[time>II]+nest(time[time>II]-II,dno+1)
      # DV
    # }
    
    list(TIME=time ,DV=DV, label2=label2)
  })
}

labeler <- function(x){
  paste(paste(c(0,x[-length(x)]), x, sep="-"), "kg")
}

wt_bands <- function(x){
  wt_bands <- as.numeric(strsplit(x,"-")[[1]])
  wt_bands <- round(wt_bands[wt_bands>0 & wt_bands<=100 & !is.na(wt_bands)]/5)*5
  wt_bands <- unique(sort(c(wt_bands[wt_bands>0],100)))
  wt_bands <- list(wt_bands, length(wt_bands), labeler(wt_bands))
  return(wt_bands)
}

AUC <- function(TIME, DV) {
  sum(diff(TIME)*(DV[-1]+DV[-length(DV)]))/2
}

simPar <- function(indDat,input){
  with(c(indDat,input,par), {
    
    CLi = (CL+0.0545*(WT-48))*exp(OM1*sqrt(ETA1))
    V2i = (V2+0.433*(WT-48))*exp(OM2*sqrt(ETA2))
    KAi = KA
    # ALAG1i = ALAG1*exp(OM4*sqrt(ETA3))
    k<-CLi/V2i
    l1<-CLi/V2i
    C1<-1/V2i
    
    DOSE <- c(input[grep("DOSE[0-9]+",names(input))][order(names(input)[grep("DOSE[0-9]+",names(input))])],800,1200,1600,2000,5,5)[[label2]]
    # F1 <- 823/(1120 + DOSE)
    # A1 = DOSE*(WT**as.numeric(1))*F1
    A1 = DOSE
    
    
    # print(input[grep("DOSE[0-9]+",names(input))])
    
    T_TMAX <- (optimize(function(time) A1*C1*(KA/(KA-l1))*(exp(-l1*time)-exp(-KA*time)), interval=c(0, 300), maximum=TRUE)$maximum)
    T_CMAX <- A1*C1*(KA/(KA-l1))*(exp(-l1*T_TMAX)-exp(-KA*T_TMAX))
    
    list(CL=CLi,V2=V2i,KA=KAi,A1=A1,AUC=A1/CLi,CMAX=T_CMAX,TMAX=T_TMAX)
  })
}

# input <- list(DOSE=100, N=20)

popFun <- function(n=20, WT=70, WTBAND) {
  # print(wt_bands(WTBAND))
  popDat <- CJ(ID=1:n)
#   popDat <- cbind(popDat,data.table(OM1=rnorm(n),OM2=rnorm(n),OM3=rnorm(n),OM4=rnorm(n),OM5=rnorm(n), WT=WT*exp(rnorm(n,0,0.1))))
  popDat <- cbind(popDat,data.table(OM1=rnorm(n),OM2=rnorm(n),OM3=rnorm(n),OM4=rnorm(n),OM5=rnorm(n), WT=runif(n, 20, 100)))
  # popDat$label2 <- cut(popDat$WT, breaks=c(-Inf,5,20,50,Inf), labels=c("<5 kg", "5-20 kg", "20-50 kg", ">50 kg"))
  popDat$label2 <- cut(popDat$WT, breaks=c(-Inf,wt_bands(WTBAND)[[1]]), labels=wt_bands(WTBAND)[[3]])
  # print(head(popDat))
  popDat[order(ID)]
}

adultFun <- function(n=20, WT=70) {
  popDat <- CJ(ID=1:n)
    popDat <- cbind(popDat,data.table(OM1=rnorm(n),OM2=rnorm(n),OM3=rnorm(n),OM4=rnorm(n),OM5=rnorm(n), WT=runif(n, WT,WT+5)))#WT*exp(rnorm(n,0,0.1))))
  # popDat <- cbind(popDat,data.table(OM1=rnorm(n),OM2=rnorm(n),OM3=rnorm(n),OM4=rnorm(n),OM5=rnorm(n), WT=rnorm(n, 0, 100)))
  # popDat$label2 <- cut(popDat$WT, breaks=c(-Inf,5,20,50,Inf), labels=c("<5 kg", "5-20 kg", "20-50 kg", ">50 kg"))
  popDat$label2 <- factor("Adult")
  
  popDat[order(ID)]
}

parFun <- function(simPop,input){
  merge(simPop,simPop[,simPar(as.list(.SD),input),by="ID"])
}

r_factor <- function(list, n = 1) {
  factor(list, levels=list)[n]
}


# simInd <- function(indDat,input,time= exp(seq(0,log(49),len=24))-1) {
#   list(TIME=time,DV=exp(-1*time))
# }

popDat <- adultFun(1000, 75)
adultPop <- parFun(popDat,list(DOSE1 = 2000, ADM=1))
adultDat <- adultPop[,simInd(as.list(.SD),list(DOSE = 0.5),times),by="ID"]
adultSum <- adultDat[,list(M=median(DV), U=quantile(DV, 0.95), L=quantile(DV, 0.05)), by="TIME"]
adultTot <- melt(adultPop[,c("ID", "AUC", "TMAX", "CMAX"),with=F],id.vars="ID")



function(input, output, session) {  
  popDat <- reactive(popFun(input$N, input$WT, input$WTBAND))
  simPop <- reactive(parFun(popDat(),reactiveValuesToList(input)))

  
  updatewt <- observe({
    updateTextInput(session, "WTBAND", value = paste(c(0,wt_bands(input$WTBAND)[[1]]), collapse="-"))
    })
  
  
#   observeEvent(input$ADM,{
#     # output$doseControls()
#   })
  
  output$doseControls <- renderUI({
    out <- list()
    x <- wt_bands(input$WTBAND)
    for (i in (1:x[[2]])) {
      # out <- c(out,list(sliderInput(paste0("DOSE",i), paste0("Dose for patients ",x[[3]][i]," [",c("mg","mg/kg")[as.numeric(input$ADM)+1], "]"),  min = 1, max = 20*70**(1-as.numeric(input$ADM)), step = 1, value = 5*70**(1-as.numeric(input$ADM)))))
      out <- c(out,list(sliderInput(paste0("DOSE",i), paste0("Dose for patients ",x[[3]][i]," [mg]"),  min = 200, max = 4000, step = 50, value = ifelse(input$WTBAND=="0-40-55-70-100",c(800,1200,1600,2000)[i], 1500))))

    }
    # print(wt_bands(input$WTBAND))
    return(out)
  })
    

  
  simDat <- reactive(transform(simPop()[,simInd(as.list(.SD),reactiveValuesToList(input),times),by="ID"]))
  simTot <- reactive({
    transform(simPop(), label=floor(WT/5)*5, label2=cut(WT, breaks=c(-Inf,wt_bands(input$WTBAND)[[1]]), labels=wt_bands(input$WTBAND)[[3]]))
  })
  
  output$pkPlot <- renderPlot({
    if (!input$sum) {
      plot <- ggplot(simDat())+geom_ribbon(data=adultSum, aes(TIME, ymax=U, ymin=L, fill="Adult 70 kg reference"), col=NA, alpha=1/5)+geom_line(aes(TIME,DV,group=ID), alpha=1/3)+coord_cartesian(xlim=c(0,48))+labs(x="Time [h]", y="Concentration [mg/L]", fill="",col="")+theme(legend.position="top")+facet_wrap(~label2)+geom_line(data=simDat()[,list(M=median(DV)), by=c("TIME","label2")], aes(TIME, M), linetype=2)+geom_line(data=adultSum, aes(TIME, M, col="Adult 70 kg reference"), linetype=2)
    } else {
      plot <- ggplot(simDat()[,list(M=median(DV), U=quantile(DV, 0.95), L=quantile(DV, 0.05)), by=c("TIME","label2")])+geom_ribbon(data=adultSum, aes(TIME, ymax=U, ymin=L, fill="Adult 70 kg reference"), col=NA, alpha=1/5)+geom_ribbon(aes(TIME, ymin=L, ymax=U), alpha=1/3)+coord_cartesian(xlim=c(0,48))+labs(x="Time [h]", y="Concentration [mg/L]", fill="", col="")+theme(legend.position="top")+facet_wrap(~label2)+geom_line(aes(TIME, M), linetype=2)+geom_line(data=adultSum, aes(TIME, M, col="Adult 70 kg reference"), linetype=2)
      
    }
    if (input$log) plot <- plot+scale_y_log10()
    print(plot)
  })
  
  output$pkComp <- renderPlot({
    # print(simDat())
    if (input$sum) {
      test <- melt(simTot()[,c("label2", "label", "AUC", "CMAX", "TMAX"), with=F], id.vars=c("label", "label2"))
      #minn <- min(test[variable==input$label]$value)
      #test <- test[variable==input$label,list(N=.N, Y=min(value)), by="label"]
      #test$Y <- minn- minn
      # print(test)
    ggplot(simTot())+geom_ribbon(data=adultTot[,list(U=quantile(value, 0.875), L=quantile(value, 0.125),WT=c(15,105))],aes(x=WT, ymin=L, ymax=U, fill="75% CI of Adults (2000 mg)"), alpha=1/8)+geom_boxplot(aes_string("label+2.5",input$label,group="label", fill="label2"))+labs(x="", y=labTab[[match(input$label,names(labTab))]],fill="")+theme(legend.position="bottom")+geom_text(data=test, aes(label=N, x=label+2.5, y=Y*0.95))
    } else {
      ggplot(simTot())+geom_point(aes_string("WT",input$label,group="label", col="label2"))+labs(x="", y=labTab[[match(input$label,names(labTab))]],fill="", col="")+theme(legend.position="bottom")+geom_ribbon(data=adultTot[,list(U=quantile(value, 0.875), L=quantile(value, 0.125),WT=c(15,105))],aes(x=WT, ymin=L, ymax=U, fill="75% CI of Adults (2000 mg)"), alpha=1/8)
    }
  })
}
  
#variable==input$label
#variable==input$label