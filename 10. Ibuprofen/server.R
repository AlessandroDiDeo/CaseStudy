#Define server for Ibuprofen in Pre-Term Neonates Application (Example 3)
#Load package libraries
	library(shiny)
	library(ggplot2)
	library(deSolve)
	library(plyr)
	library(reshape2)
	library(compiler)
	library(doParallel)

	ICREF <- function(IC50, target) -((IC50*target)/(-100 + target))
	
#Code for functions and variables which are not reactive (not dependent on "input$X")
#Setting up cores to run parallel processes, thus increasing speed
#Set up a cluster of cores to run the application over
	cl <- makePSOCKcluster(3)
#detectCores() searches for the number of cores that the local machine has
#Contents with makePSOCKcluster brackets can be changed to a whole number if you
#want to assign an exact number	
	
#List packages that are required to be sent to each core for the parallel process
#The foreach package always needs to be included
#This example uses the .parallel argument in ddply which calls a function that uses
#lsoda from the deSolve package
	clusterEvalQ(cl, list(library(deSolve), library(foreach)))
	
#Registers the parallel backend with the foreach package (automatically loaded when
#doParallel is loaded)	
	registerDoParallel(cl)	

#ggplot2 theme for plotting
	theme_custom <- theme_set(theme_grey(18))
	
#TIME range - times where a concentration will be calculated
	TIME <- seq(from = 0, to = 72, by = 0.25)
	
#----------------------------------------------------------------------------------------	
#Define user-input dependent functions for output
	shinyServer(function(input, output) {
		
#Reactive expression to generate a reactive data frame
#This is called whenever the input changes
	all.data <- reactive({
	
		#Assign probabilities for 95% confidence intervals
			CIlow <- 0.25
			CIhi <- 0.75
		
		#Function for calculating median, upper and lower confidence intervals for x
		#Where x will be concentrations for R- and S-ibuprofen
		sumfuncx <- function(x) {
			stat1 <-  median(x)
			stat2 <-  quantile(x, probs=CIlow, names=F) 
			stat3 <-  quantile(x, probs=CIhi, names=F)
			stat4 <-  length(x)
			result <- c("median"=stat1, "low"=stat2, "hi"=stat3, "n"=stat4)
			result
		}
	
#----------------------------------------------------------------------------------------	
	#Create a dataframe with ID and parameter values for each individual	
	#Number of individuals
		n <- input$POP
		par.data <- seq(from = 1, to = n, by = 1)
		par.data <- data.frame(par.data)
		names(par.data)[1] <- "ID"
		
	#Make a parameter vector for input into DES function				
		par.data$K21 <- 0.014	#hour^-1
		
	#Define population values
		POPV1 <- 173*(1+input$V/100)		#mL/kg
		POPV2 <- 306*(1+input$V/100)		#mL/kg
		# POPKE1 <- 0.020		#hour^-1
		# POPKE2 <- 0.069		#hour^-1
		POPKE1 <- (173*0.020)/POPV1		#hour^-1
		POPKE2 <- (306*0.069)/POPV2		#hour^-1
		
	#Define population parameter variability
		ETA1 <- rnorm(n, mean = 0, sd = 26/100)
		ETA2 <- rnorm(n, mean = 0, sd = 95/100) #input$ETA2/100)
		ETA3 <- rnorm(n, mean = 0, sd = 58/100)
		ETA4 <- rnorm(n, mean = 0, sd = 26/100) #input$ETA4/100)
		
	#Simulate individual values	
		par.data$V1 <- POPV1*exp(ETA1)		#mL/kg
		par.data$V2 <- POPV2*exp(ETA2)		#mL/kg
		par.data$KE1 <- POPKE1*exp(ETA3)	#hour^-1
		par.data$KE2 <- POPKE2*exp(ETA4)	#hour^-1
	
#------------------------------------------------------------------------------------------	
	#Select between 3 bolus regimen or bolus/continuous reg	
		SELECT <- input$SELECT
				
	#Input dosing data for loading dose (LDOSE = loading dose)			
		LDOSE <- 1000*input$LDOSE/2		#microg/kg
		LTinf <- 0.25					#Loading dose infusion duration (hours)
		LRATE <- LDOSE/LTinf			#Loading infusion rate
		
		LTIMEinf <- c(0,LTinf,100)		#Vector marking infusion's time events, function
										#works long after infusion is finished
		LRATEinf <- c(LRATE,0,0)		#Vector marking infusion's rates
		
	#Define an interpolation function that returns rate when given time - "const"
		Lstep.doseinf <- approxfun(LTIMEinf, LRATEinf, method = "const")

		
	#Input dosing data for 2 bolus doses (BDOSE = bolus dose)				
		BDOSE <- 1000*input$BDOSE/2		#microg/kg
		BTinf <- 0.25					#Bolus dose infusion duration (hours)
		BRATE <- BDOSE/BTinf			#Bolus infusion rate
		
		II <- input$II
		# II <- 12
		nDo <- floor(100/II)
		
		BTIMEinf <- sort(c(0, seq(II, 100, by=II), seq(II, 100, by=II)+BTinf, 100+BTinf))
		BRATEinf <- c(0, rep(c(BRATE, 0),nDo), 0)
		
		print(BTIMEinf)
		print(BRATEinf)

		#BTIMEinf <- c(0,	24,BTinf+24,	48,BTinf+48,	100)	#Vector marking infusion's time events, function
		#															#works long after infusion is finished
		#BRATEinf <- c(0,	BRATE,0,		BRATE,0,		0)		#Vector marking infusion's rates
		
	#Define an interpolation function that returns rate when given time - "const"
		Bstep.doseinf <- approxfun(BTIMEinf, BRATEinf, method = "const")

		
	#Creating continuous infusion (for 72 hours)
		CDOSE <- 1000*input$CDOSE/2		#microg/kg		
		CTinf <- 72 					#Continuous infusion duration (hours)
		CRATE <- CDOSE/CTinf			#Continuous infusion rate
					
		CTIMEinf <- c(0,0.25,CTinf+0.25,100)	#Vector marking infusion's time events, function
												#works long after infusion is finished
		CRATEinf <- c(0,CRATE,0,0)				#Vector marking infusion's rates
				
	#Define an interpolation function that returns rate when given time - "const"
		Cstep.doseinf <- approxfun(CTIMEinf, CRATEinf, method = "const")

		
	#Input post-natal age
		AGE <- input$AGE	#hours
			
	#Differential equations describing R- and S-ibuprofen		
		DES <- function(T, A, THETA) {
			
			#If Loading-Bolus-Bolus regimen is selected
			if (SELECT == 1) {		
			RateL <- Lstep.doseinf(T)
			RateCB <- Bstep.doseinf(T)
			}
			
			#If Loading-Continuous Infusion regimen is selected
			if (SELECT == 2) {
			RateL <- Lstep.doseinf(T)
			RateCB <- Cstep.doseinf(T)
			}
								
			KE1 <- THETA[1]
			KE2 <- THETA[2]
			K21 <- THETA[3]
				
			dA <- vector(length = 2)					
			dA[1] = RateL + RateCB +K21*A[2] -KE1*A[1]		#S-ibuprofen compartment
			dA[2] = RateL + RateCB -K21*A[2] -(KE2+0.155*((T+AGE)/24))*A[2] #R-ibuprofen compartment
						
			list(dA)			
		}
	
	#Compile DES function
		DES.cmpf <- cmpfun(DES)
	
#----------------------------------------------------------------------------------------	
	#Function for simulating concentrations for the ith patient
		simulate.conc <- function(par.data) {
		
		#Set initial conditions in each compartment
			A_0 <- c(A1 = 0, A2 = 0)
		
		#List of parameter values	
			THETAlist <- c("KE1"=par.data$KE1,
							"KE2"=par.data$KE2,
							"K21"=par.data$K21)	
							
		#Run differential equation solver for simulated parameter data	
			sim.data <- lsoda(A_0, TIME, DES.cmpf, THETAlist)
			sim.data <- as.data.frame(sim.data)	
		}
		
	#Compile simulate.conc function	
		simulate.conc.cmpf <- cmpfun(simulate.conc)
		
	#Apply simulate.conc.cmpf function to each individual in par.data
	#ID is provided as a variable in which each value has simulate.conc.cmpf applied to
	#V1 and V2 are required to be here to preserve their values, so concentrations can be
	#calculated later
		sim.data <- ddply(par.data, .(ID,V1,V2), simulate.conc.cmpf, .parallel = TRUE)

	#Calculate concentration of R- and S-ibuprofen		
		sim.data$CONCS <- sim.data$A1/sim.data$V1
		sim.data$CONCR <- sim.data$A2/sim.data$V2
		
		statsCONCS <- ddply(sim.data, .(time), function(sim.data) sumfuncx(sim.data$CONCS))
		names(statsCONCS)[c(2,3,4)] <- c("Smedian","Slow","Shi")
		statsCONCR <- ddply(sim.data, .(time), function(sim.data) median(sim.data$CONCR))	
		names(statsCONCR)[2] <- "Rmedian"
		
	#Combine both datasets
		sim.data <- merge(statsCONCS,statsCONCR,by=c("time"),all=T)
	
	})	#Brackets closing "reactive" expression

#----------------------------------------------------------------------------------------
#Generate a plot of the data
#Also uses the inputs to build the plot
	output$plotCONC <- renderPlot({	
	
		plotobj <- ggplot(all.data())
		plotobj <- plotobj + geom_abline(aes(slope = 0, intercept = 2.5*4), linetype = "dashed", size = 1)
		plotobj <- plotobj + geom_abline(aes(slope = 0, intercept = 16.5*4), linetype = "dashed", size = 1)
		if (input$REF!=80) plotobj <- plotobj + geom_abline(aes(slope = 0, intercept = ICREF(2.5, input$REF)), linetype = "dashed", size = 1, col="darkred")
		if (input$REF!=80) plotobj <- plotobj + geom_abline(aes(slope = 0, intercept = ICREF(16.5, input$REF)), linetype = "dashed", size = 1, col="darkred")
		plotobj <- plotobj + geom_line(aes(x = time, y = Smedian, colour = 'S-ibuprofen'), size = 1)
		plotobj <- plotobj + geom_line(aes(x = time, y = Rmedian, colour = "R-ibuprofen"), size = 1)
		plotobj <- plotobj + geom_ribbon(aes(x = time, ymin = Slow, ymax = Shi), fill = "red", alpha = 0.3)
		plotobj <- plotobj + scale_y_continuous("Concentration (microg/mL) \n")
		plotobj <- plotobj + scale_x_continuous("\nTime (hours)", breaks=c(0,8,16,24,32,40,48,56,64,72))
		plotobj <- plotobj + scale_color_manual('', values = c('blue', 'red'), breaks = c('S-ibuprofen', 'R-ibuprofen'))	
		plotobj <- plotobj + annotate("text", x = 65, y = 3.5*4, label = "IC80 COX-1", colour = "black", size = 4)
		plotobj <- plotobj + annotate("text", x = 65, y = 17.5*4, label = "IC80 COX-2", colour = "black", size = 4)
		if (input$REF!=80) plotobj <- plotobj + annotate("text", x = 5, y = ICREF(2.5, input$REF)+4, label = paste0("IC",input$REF," COX-1"), size = 4, col="darkred")
		if (input$REF!=80) plotobj <- plotobj + annotate("text", x = 5, y = ICREF(16.5, input$REF)+4, label = paste0("IC",input$REF," COX-2"), size = 4, col="darkred")
		plotobj <- plotobj + theme_bw()
		print(plotobj)
	
	})	#Brackets closing "renderPlot" function
	
	output$RATE <- renderText({
		paste("Infusion rate =", signif(input$CDOSE/72, digits = 3) ,"mg/kg/hr")
	
	})	#Brackets closing "renderText" function
	
})	#Brackets closing "shinyServer" function
