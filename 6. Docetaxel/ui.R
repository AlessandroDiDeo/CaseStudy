library(shiny)

shinyUI(navbarPage(
  
  titlePanel("Docetaxel"),
  
  tabPanel("Plot",
           sidebarLayout(
             sidebarPanel(
               sliderInput("DOSE1", "Dose C1236T Heteroz. (mg/m2)",  min = 5, max = 750, step = 5, value = 75),
               sliderInput("DOSE2", "Dose C1236T Homoz. (mg/m2)",  min = 5, max = 750, step = 5, value = 75),
               sliderInput("DOSE", "Dose Wild Type (mg/m2)",  min = 5, max = 750, step = 5, value = 75),
               # sliderInput("pCL",  "Percentage change in Clearance", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pV1", "Percentage change in Central Volume of Distribution", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pV2", "Percentage change in Peripheral Volume of Distribution 2", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pV3", "Percentage change in Peripheral Volume of Distribution 3", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pQ2", "Percentage change in Intercompartmental Clearance 2", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pQ3", "Percentage change in Intercompartmental Clearance 3", min = 10,  max = 500, step=5, value = 100),
               sliderInput("pD1", "Infusion duration (hours)", min = 0.1,  max = 12, step=0.1, value = 3),
               sliderInput("II","Dose interval (days)",min=0.5,max=7*4,value=0.5,step=0.5),
               sliderInput("ADDL","Number of dose events",min=1,max=10,value=2,step=1),
               sliderInput("tps","Number of measurements within dosing interval",min=2,max=30,value=12,step=1),
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
               selectInput("COR1","Y-axis:",choices=c("Clearance"="CL","Central Volume of Distribution"="V1","Peripheral Volume of Distribution 1"="V2","Peripheral Volume of Distribution 2"="V3","Intercompartmental Clearance 1"="Q2","Intercompartmental Clearance 2"="Q3","Steady State Concentration"="CSS","Area Under the Curve"="AUC","Maximum Concentration"="CMAX","Time of maximum Concentration"="TMAX","Weight"="WT","Height"="HT","Age"="AGE","Serum Creatinine"="SCR","Albumin"="ALB","Lean Body Mass"="LBM","Body Mass Index"="BMI","Creatinine Clearance"="CRCL","BSA"="BSA","Genotype"="GEN"),selected="CL"),
               selectInput("COR2","X-axis:",choices=c("Clearance"="CL","Central Volume of Distribution"="V1","Peripheral Volume of Distribution 1"="V2","Peripheral Volume of Distribution 2"="V3","Intercompartmental Clearance 1"="Q2","Intercompartmental Clearance 2"="Q3","Steady State Concentration"="CSS","Area Under the Curve"="AUC","Maximum Concentration"="CMAX","Time of maximum Concentration"="TMAX","Weight"="WT","Height"="HT","Sex"="SEX","Age"="AGE","Serum Creatinine"="SCR","Albumin"="ALB","Lean Body Mass"="LBM","Body Mass Index"="BMI","Creatinine Clearance"="CRCL","BSA"="BSA","Genotype"="GEN","Total dataset"="TOTAL"),selected="WT"),
               checkboxInput("lm","Linear regression"),
               selectInput("STRAT","Stratify by:",choices=c("None"=0,Sex="SEX",Age="AGE",Genotype="GEN")),
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
