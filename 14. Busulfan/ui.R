library(shiny)

shinyUI(navbarPage(
  
  titlePanel("Busulfan"),
  
  tabPanel("Plot",
           sidebarLayout(
             sidebarPanel(
               sliderInput("DOSE", "Dose in children (< 9 kg)",  min = 0.1, max = 2, step = 0.1, value = 1000),
               sliderInput("DOSE", "Dose (> 9 kg)",  min = 0.1, max = 2, step = 0.1, value = 1000),
               # sliderInput("pCL",  "Percentage change in Clearance", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pV2", "Percentage change in Central Volume of Distribution", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pV3", "Percentage change in Peripheral Volume of Distribution", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pQ", "Percentage change in Intercompartmental Clearance", min = 10,  max = 500, step=5, value = 100),
               sliderInput("pD1", "Infusion duration (hours)", min = 0.1,  max = 24, step=0.1, value = 0.5),
               sliderInput("II","Dose interval",min=4,max=48,value=12,step=1),
               sliderInput("ADDL","Number of dose events",min=1,max=10,value=1,step=1),
               sliderInput("tps","Number of measurements within dosing interval",min=2,max=20,value=12,step=1),
               checkboxInput("sum","Summarize"),
               checkboxInput("log","Log scale"),
               checkboxInput("enh","Show Sampling points"),
               submitButton("Update Graph", icon=icon("refresh", class="fa-spin"))
             ),
             mainPanel(
               plotOutput("pkPlot",width="100%",height="600")
             )
           )
  ),
  tabPanel("Correlations",
           sidebarLayout(
             sidebarPanel(
               selectInput("COR1","Y-axis:",choices=c("Clearance"="CL","Central Volume of Distribution"="V2","Peripheral Volume of Distribution"="V3","Intercompartmental Clearance"="Q","Area Under the Curve"="AUC","Maximum Concentration"="CMAX","Time of maximum Concentration"="TMAX","Time Above MIC at SS" = "TMIC","Weight"="WT","Height"="HT","Age"="AGE","Serum Creatinine"="SCR","Albumin"="ALB","Lean Body Mass"="LBM","Body Mass Index"="BMI","Creatinine Clearance"="CRCL","BSA"="BSA"),selected="CL"),
               selectInput("COR2","X-axis:",choices=c("Clearance"="CL","Central Volume of Distribution"="V2","Peripheral Volume of Distribution"="V3","Intercompartmental Clearance"="Q","Area Under the Curve"="AUC","Maximum Concentration"="CMAX","Time of maximum Concentration"="TMAX","Time Above MIC at SS" = "TMIC","Weight"="WT","Height"="HT","Sex"="SEX","Age"="AGE","Serum Creatinine"="SCR","Albumin"="ALB","Lean Body Mass"="LBM","Body Mass Index"="BMI","Creatinine Clearance"="CRCL","BSA"="BSA","Total dataset"="TOTAL"),selected="WT"),
               checkboxInput("lm","Linear regression"),
               selectInput("STRAT","Stratify by:",choices=c("None"=0,Sex="SEX",Age="AGE")),
               submitButton("Update Graph", icon=icon("refresh", class="fa-spin"))
             ),
             mainPanel(
               plotOutput("corPlot"),
               verbatimTextOutput("lm_info")
             )
           )
  ),
  tabPanel("Dataset",
           sidebarLayout(
             sidebarPanel(
               sliderInput("N",
                           "Number of subjects",
                           min = 5,
                           max = 1200,
                           step = 5,
                           value = 20),
               downloadButton('dlDat', 'Download .csv')
             ),
             mainPanel(
               dataTableOutput('popSum')
             )
           )
  )
)
)
