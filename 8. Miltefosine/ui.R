library(shiny)

shinyUI(navbarPage(
  
  titlePanel("Miltefosine"),
  
  tabPanel("Plot",
           sidebarLayout(
             sidebarPanel(
               sliderInput("DOSE3", "Dose in Adults (mg)",  min = 10, max = 250, step = 10, value = 150),
               sliderInput("II3","Dose interval in Adults (hours)",min=4,max=48,value=24,step=1),
               
               sliderInput("DOSE2", "Dose in Adolescents (mg/kg)",  min = 0.5, max = 50, step = 0.1, value = 2.5),
               sliderInput("II2","Dose interval in Adolescents (hours)",min=4,max=48,value=24,step=1),
               
               sliderInput("DOSE1", "Dose in Children (mg/kg)",  min = 0.5, max = 50, step = 0.1, value = 2.5),
               sliderInput("II1","Dose interval in Children (hours)",min=4,max=48,value=24,step=1),
               
               # sliderInput("pCL",  "Percentage change in Clearance", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pV2", "Percentage change in Central Volume of Distribution", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pV3", "Percentage change in Peripheral Volume of Distribution", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pQ", "Percentage change in Intercompartmental Clearance", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pKA", "Percentage change in Absorption rate", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("ADDL","Number of dose events",min=1,max=80,value=20,step=1),
               sliderInput("ADD","Days of treatment",min=2,max=50,value=20,step=1),
               
               sliderInput("tps","Number of measurements per day",min=2,max=20,value=12,step=1),
               
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
               selectInput("COR1","Y-axis:",choices=c("Clearance"="CL","Central Volume of Distribution"="V2","Peripheral Volume of Distribution"="V2","Intercompartmental Clearance"="Q","Area Under the Curve"="AUC","Area Under the Curve (24h SS)"="AUC24","Maximum Concentration at SS"="CMAX","Time of maximum Concentration"="TMAX","Weight"="WT","Height"="HT","Age"="AGE","Serum Creatinine"="SCR","Albumin"="ALB","Lean Body Mass"="LBM","Body Mass Index"="BMI","Creatinine Clearance"="CRCL","BSA"="BSA"),selected="CL"),
               selectInput("COR2","X-axis:",choices=c("Clearance"="CL","Central Volume of Distribution"="V2","Peripheral Volume of Distribution"="V2","Intercompartmental Clearance"="Q","Area Under the Curve"="AUC","Area Under the Curve (24h SS)"="AUC24","Maximum Concentration at SS"="CMAX","Time of maximum Concentration"="TMAX","Weight"="WT","Height"="HT","Sex"="SEX","Age"="AGE","Serum Creatinine"="SCR","Albumin"="ALB","Lean Body Mass"="LBM","Body Mass Index"="BMI","Creatinine Clearance"="CRCL","BSA"="BSA","Total Dataset"="TOTAL"),selected="WT"),
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
                             max = 2000,
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
