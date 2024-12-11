library(shiny)

shinyUI(navbarPage(
  
  titlePanel("Case Study 1"),
  
  tabPanel("Pharmacokinetics",
           sidebarLayout(
             sidebarPanel(
               sliderInput("N",
                           "Number of subjects",
                           min = 5,
                           max = 1200,
                           step = 5,
                           value = 100),
               sliderInput("DOSE1", "Dose (mg)",  min = 1, max = 200, step = 1, value = 10),
               sliderInput("II","Dose interval (hours)",min=4,max=48,value=12,step=1),
               sliderInput("ADDL","Number of dose events",min=1,max=10,value=1,step=1),
               checkboxInput("log","Log scale")
             ),
             mainPanel(
               plotOutput("pkPlot",height="300"),
               #plotOutput("AUC",height="300"),
               #plotOutput("Cmax",height="300")
             )
           )
  ),
  tabPanel("PKPD indices",
           sidebarLayout(
             sidebarPanel(
               selectInput("COR1","Parameter:",choices=c("Time above MIC"="TMIC", 
                                                         "Probability of target attainment" = "PTA"),selected="TMIC"),
               textInput("MIC", 'MIC values (in mg/L)', "0.125,0.25,0.5,1,2,4,8,16,32"),
               sliderInput("thres","Target T>MIC (%)",min=0,max=100,value=60,step=5)
             ),
             mainPanel(
               plotOutput("PKPD")
             )
           )
  )
)
)
