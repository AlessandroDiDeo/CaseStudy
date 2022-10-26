library(shiny)

shinyUI(navbarPage(
  
  titlePanel("Amikacin"),
  
  tabPanel("Plot",
           sidebarLayout(
             sidebarPanel(
               sliderInput("DOSE", "Dose (mg/kg)",  min = 1, max = 200, step = 1, value = 20),
               # sliderInput("pCL",  "Percentage change in Clearance", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pV2", "Percentage change in Central Volume of Distribution", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pV3", "Percentage change in Peripheral Volume of Distribution", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pQ", "Percentage change in Intercompartmental Clearance", min = 10,  max = 500, step=5, value = 100),
               sliderInput("pD1", "Infusion duration (hours)", min = 0.1,  max = 12, step=0.05, value = 0.5),
               sliderInput("II","Dose interval (hours)",min=4,max=48,value=12,step=1),
               sliderInput("ADDL","Number of dose events",min=1,max=10,value=1,step=1),
               sliderInput("tps","Number of measurements within dosing interval",min=2,max=20,value=12,step=1),
               checkboxInput("sum","Summarize"),
               checkboxInput("log","Log scale"),
               checkboxInput("enh","Show Sampling points"),
               submitButton("Update Graph", icon=icon("refresh", class="fa-spin"))
             ),
             mainPanel(
               plotOutput("pkPlot",height="600")
             )
           )
  ),
  tabPanel("Correlations",
           sidebarLayout(
             sidebarPanel(
               selectInput("COR1","Y-axis:",choices=c("Clearance"="CL","Volume of Distribution"="V2","Area Under the Curve"="AUC","Maximum Concentration"="CMAX","Time of maximum Concentration"="TMAX","C1h/MIC Ratio"="CMIC","Weight"="WT","Height"="HT","Age"="AGE","Serum Creatinine"="SCR","Albumin"="ALB","Lean Body Mass"="LBM","Body Mass Index"="BMI","Creatinine Clearance"="CRCL","BSA"="BSA"),selected="CL"),
               selectInput("COR2","X-axis:",choices=c("Clearance"="CL","Volume of Distribution"="V2","Area Under the Curve"="AUC","Maximum Concentration"="CMAX","Time of maximum Concentration"="TMAX","C1h/MIC Ratio"="CMIC","Weight"="WT","Height"="HT","Sex"="SEX","Age"="AGE","Serum Creatinine"="SCR","Albumin"="ALB","Lean Body Mass"="LBM","Body Mass Index"="BMI","Creatinine Clearance"="CRCL","BSA"="BSA","Total dataset"="TOTAL"),selected="WT"),
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
