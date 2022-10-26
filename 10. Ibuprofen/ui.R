#Define UI for Ibuprofen in Pre-Term Neonates Application (Example 3)
fixedPage(

	titlePanel("Ibuprofen"),
	#Sidebar panel with widgets
	sidebarLayout(
		sidebarPanel(
			
			#Selection box for dosing regimen
			radioButtons("SELECT",
						"Dosing Regimen:",
						choices = list("Loading-Bolus-Bolus" = 1,
										"Loading-Continuous Infusion" = 2),
						selected = 1),
			sliderInput("POP",
			            "Number of subjects:",
			            min = 10,
			            max = 600,
			            value = 50,
			            step = 10),
						
				br(),	#Add a blank break between widgets
										
				conditionalPanel(condition = "input.SELECT == 1",
						
			#Slider input for dose
			#Bolus dose for Loading-Bolus-Bolus regimen
			sliderInput("BDOSE",
						"Bolus Dose at 24 and 48 hours (mg/kg):",
						min = 0,
						max = 40,
						value = 5,
						step = 1)
						
				),	#Brackets closing "conditionalPanel"
				
				conditionalPanel(condition = "input.SELECT == 2",
										
			#Slider input for dose
			#Dose to be infused over 72 hours for Loading-Continuous Infusion regimen
			sliderInput("CDOSE",
						"Dose continuously infused over 72 hours (mg/kg):",
						min = 0,
						max = 80,
						value = 12,
						step = 1),

			#Numeric output for infusion rate
			textOutput("RATE")
				
				),	#Brackets closing "conditionalPanel"
			
			#Slider input for loading dose
			sliderInput("LDOSE",
						"Loading Dose (mg/kg):",
						min = 0,
						max = 40,
						value = 10,
						step = 1),
			
			sliderInput("II",
			            "Dosing Interval (h)",
			            min = 0,
			            max = 48,
			            value = 24,
			            step = 3),			
			#Slider input for age
			sliderInput("AGE",
						"Postnatal Age (Hours):",
						min = 0,
						max = 72,
						value = 24,
						step = 4),
			sliderInput("REF",
			            "Reference inhibitory concentration:",
			            min = 5,
			            max = 95,
			            value = 80,
			            step = 5),
			sliderInput("V",
			            "Change in Volume (%):",
			            min = -99,
			            max = 100,
			            value = 0,
			            step = 1),
			#Slider input for ETA1
			# sliderInput("ETA1",
			#             "S-ibuprofen volume of distribution BSV (CV%):",
			#             min = 0,
			#             max = 200,
			#             value = 26,
			#             step = 10),
			# #Slider input for ICREF
			
			#Slider input for ETA2
			# sliderInput("ETA2",
			#             "R-ibuprofen volume of distribution BSV (CV%):",
			#             min = 0,
			#             max = 200,
			#             value = 95,
			#             step = 10),
			
			#Slider input for ETA3
			# sliderInput("ETA3",
			#             "S-ibuprofen clearance BSV (CV%):",
			#             min = 0,
			#             max = 200,
			#             value = 58,
			#             step = 10),
			# 
			#Slider input for ETA4
			# sliderInput("ETA4",
			#             "R-ibuprofen clearance BSV (CV%):",
			#             min = 0,
			#             max = 200,
			#             value = 26,
			#             step = 10),
			
			br(),
			
			#Button to initiate simulation
			submitButton("Simulate"),
		
		align = "left"),	#Brackets closing "sidebarPanel"
							
		mainPanel(
		
			#Plot output for concentration-time profile
			plotOutput("plotCONC", height = 650, width = 900),
		
		align = "center")	#Brackets closing "mainPanel"
			
	)	#Brackets closing "sidebarLayout"

)	#Brackets closing "fixedPage"
