library(shiny)

shinyUI(navbarPage(
  
  titlePanel("Mycophenolic acid"),
  
  tabPanel("Plot",
           sidebarLayout(
             sidebarPanel(
               sliderInput("DOSE", "Dose (mg) (w/o Cyclosporin)",  min = 100, max = 10000, step = 50, value = 1000),
               sliderInput("DOSE2", "Dose (mg) (w/ Cyclosporin)",  min = 100, max = 10000, step = 50, value = 1000),
               sliderInput("dCSA", "Cyclosporin dose (mg/kg)",  min = 0, max = 100, step = 5, value = 10),
               # sliderInput("pCL",  "Percentage change in Clearance", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pV2", "Percentage change in Central Volume of Distribution", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pV3", "Percentage change in Peripheral Volume of Distribution", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pQ", "Percentage change in Intercompartmental Clearance", min = 10,  max = 500, step=5, value = 100),
               # sliderInput("pKA", "Percentage change in Absorption rate", min = 10,  max = 500, step=5, value = 100),
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
               selectInput("COR1","Y-axis:",choices=c("Clearance"="CL","Central Volume of Distribution"="V2","Peripheral Volume of Distribution"="V2","Intercompartmental Clearance"="Q","Area Under the Curve"="AUC","Maximum Concentration"="CMAX","Trough Concentration"="CTROUGH","Time of maximum Concentration"="TMAX","Weight"="WT","Height"="HT","Age"="AGE","Serum Creatinine"="SCR","Albumin"="ALB","Lean Body Mass"="LBM","Body Mass Index"="BMI","Creatinine Clearance"="CRCL","BSA"="BSA","Cyclosporine"="CSA"),selected="CL"),
               selectInput("COR2","X-axis:",choices=c("Clearance"="CL","Central Volume of Distribution"="V2","Peripheral Volume of Distribution"="V2","Intercompartmental Clearance"="Q","Area Under the Curve"="AUC","Maximum Concentration"="CMAX","Time of maximum Concentration"="TMAX","Weight"="WT","Height"="HT","Sex"="SEX","Age"="AGE","Serum Creatinine"="SCR","Albumin"="ALB","Lean Body Mass"="LBM","Body Mass Index"="BMI","Creatinine Clearance"="CRCL","BSA"="BSA","Cyclosporine"="CSA","Total Dataset"="TOTAL"),selected="WT"),
               checkboxInput("lm","Linear regression"),
               selectInput("STRAT","Stratify by:",choices=c("None"=0,Sex="SEX",Age="AGE",CSA="CSA")),
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
