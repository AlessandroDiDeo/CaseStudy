library(shiny)

shinyUI(navbarPage(
  
  titlePanel("Pyrazinamide"),
  
  tabPanel("Plot",
           sidebarLayout(
             sidebarPanel(
               sliderInput("N", "Number of patients",  min = 1, max = 5000, step = 1, value = 10),
               # radioButtons("ADM", "Dosing in:", c("mg/kg" = 1, "mg" = 0), selected = 1, inline = TRUE),
               textInput("WTBAND", 'Weight band partitions (in kg)', "0-40-55-70-100"),
               uiOutput("doseControls"),
               # sliderInput("DOSE7", "Dose for patients <5kg [mg/kg]",  min = 1, max = 30, step = 1, value = 5),
               # sliderInput("DOSE6", "Dose for patients 5-20kg [mg/kg]",  min = 1, max = 30, step = 1, value = 5),
               # sliderInput("DOSE3", "Dose for patients 20-50kg [mg/kg]",  min = 1, max = 30, step = 1, value = 5),
               # sliderInput("DOSE4", "Dose for patients >50kg [mg/kg]",  min = 1, max = 30, step = 1, value = 5),
               # radioButtons("select_time","Sampling times:", c("Specified Interval" = 0 ,"User Specified" = 1, "Random" = 2), selected = 1),
               # sliderInput("II","Dose interval (hours)",min=4,max=48,value=12,step=1),
               # sliderInput("ADDL","Number of dose events",min=1,max=10,value=1,step=1)
               # sliderInput("WT","Population weight (± 5)[kg]",min=0.5,max=100,value=20,step=0.5),
               # textInput("simtimes", 'User specified sampling times', "0, 4, 20, 48"),
               selectInput("label","Compare:",choices=c("Area Under the Curve"="AUC","Maximum Concentration"="CMAX","Time of maximum Concentration"="TMAX"),selected="TMAX"),
               # numericInput("nrand", "Number of random samples", min=1, max=24, value=4),
               # submitButton("update", "Update Graph")
               checkboxInput("sum","Summarize", value=TRUE),
               checkboxInput("log","Log scale"),
               submitButton("Update Graph", icon=icon("refresh", class="fa-spin"))
             ),
             mainPanel(
               plotOutput("pkPlot"),
               plotOutput("pkComp")
             )
           )
  )
  
)
)

#0 - 10, 10 - 40, >40
