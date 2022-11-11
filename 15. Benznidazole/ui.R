library(shiny)
library(shinydashboard)
library(deSolve)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)

source('ddw_kineto_models.R')

dashboardPage(title = "T. Cruzi model",
              dashboardHeader(title = h3("T. Cruzi model", style = "font-family:'helvetica'; color:'aqua';"), titleWidth = NULL
              ),
              dashboardSidebar(width = NULL,
                               sidebarMenu(
                                 menuItem("Graphs",  icon = icon("line-chart"), selected=TRUE, tabName="pkpd_sim1"),
                                 menuItem("Model scheme",  icon = icon("map-o"), selected=FALSE, tabName = 'scheme')
                               ),
                               fluidPage(h5('', style = "padding:260px;"),
                                         hr(style = "border-color:gray;"),
                                         h5('Model by: S. D\'Agate, P. Healy and O. Della Pasqua', style = "color:gray;"),
                                         h5('Created by: S. D\'Agate and T. Leja', style = "color:gray;"))
              ),
              dashboardBody(
                tabItems(
                  # First tab content
                  tabItem(tabName = "pkpd_sim1",
                          fluidRow(
                            column(width=9,
                                   tabBox(title = '', selected = 'Classic output', width = NULL,
                                          tabPanel(title = 'Classic output', width=NULL,
                                                   solidHeader = F, plotOutput("plot_sim1", height=627)),
                                          tabPanel(title = 'Detailed output', width=NULL,
                                                   solidHeader = F, plotOutput("plot_sim2", height=627))
                                   )),
                            column(width=3,
                                   box(width = NULL, 
                                       radioButtons('mode_switch', label = h4(strong('Choose experimental design')),
                                                    choices = 
                                                      list('Standard experimental design (constant concentration)' = 1,
                                                           'Simulated experimental design' = 2),
                                                    selected = 1)),
                                   tabBox(width=NULL,
                                          tabPanel('Dose',
                                                   conditionalPanel(condition = 'input.mode_switch == 1',
                                                                    fluidRow(column(6, numericInput("conc_1", HTML("Concentration [&mu;M]"), 
                                                                                                    value = 10)))),
                                                   sliderInput("interval_1", "Dose interval [days]", 
                                                               min = 0, max = 7, value = c(0,3),
                                                               step = 1, ticks = F),
                                                   conditionalPanel(condition = 'input.mode_switch == 2',
                                                                    fluidRow(
                                                                      column(6, numericInput("dose_1", HTML("Dose [mg]"), 
                                                                                             value = 10)),
                                                                      column(6, numericInput("molmass", HTML("Molar mass [g/mol]"), 
                                                                                             value = 260.249)),
                                                                      column(6, numericInput("ndose_1", "Doses/day", 
                                                                                             value = 1))
                                                                    ),
                                                                    hr(),
                                                                    sliderInput("CL", HTML('CL [L/hr]'), 
                                                                                min = 0.1, max = 200, 
                                                                                value = 1.73, step = 0.01, ticks = F),
                                                                    sliderInput("V", HTML('V [L]'), 
                                                                                min = 0.1, max = 100, 
                                                                                value = 89.6, step = 0.01, ticks = F),
                                                                    sliderInput("KA", HTML('K<sub>A</sub> [hr<sup>-1</sup>]'), 
                                                                                min = 0.1, max = 100, 
                                                                                value = 1.15, step = 0.01, ticks = F)
                                                   )),
                                          tabPanel('Response',
                                                   fluidRow(column(12, h4(strong(HTML('K<sub>1</sub>')), align = 'center')),
                                                            column(6, sliderInput("ic50k1_1", HTML('IC<sub>50</sub> [&mu;M]'), 
                                                                                  min = 0, max = 200, value = 4.26, step=0.01,
                                                                                  ticks=F)),
                                                            column(6, sliderInput("tek1", HTML('T<sub>delay</sub> [hr]'), 
                                                                                  min = 0.0001, max = 168, value = 0.0001, 
                                                                                  step=1, ticks=F))
                                                   ),
                                                   hr(),
                                                   fluidRow(column(12, h4(strong(HTML('K<sub>d</sub>')), align = 'center')),
                                                            column(6, sliderInput("emax_d", HTML('E<sub>max</sub>'), 
                                                                                  min = 0, max = 200, value = 77.10, step=0.1,
                                                                                  ticks=F)),
                                                            column(6, sliderInput("ic50kd_1", HTML('EC<sub>50</sub> [&mu;M]'), 
                                                                                  min = 0, max = 200, value = 14.73, step=0.01,
                                                                                  ticks=F)),
                                                            column(6, sliderInput("tekd", HTML('T<sub>delay</sub> [hr]'), 
                                                                                  min = 0.0001, max = 168, value = 0.0001, 
                                                                                  step=1, ticks=F))
                                                   ),
                                                   hr(),
                                                   fluidRow(column(12, h4(strong(HTML('K<sub>in</sub>')), align = 'center')),
                                                            column(6, sliderInput("ic50kin_1", HTML('IC<sub>50</sub> [&mu;M]'),
                                                                                  min = 0, max = 500, value = 500, step=0.01,
                                                                                  ticks=F)),
                                                            column(6, sliderInput("tekin", HTML('T<sub>delay</sub> [hr]'), 
                                                                                  min = 0.0001, max = 168, value = 0.0001, 
                                                                                  step=1, ticks=F))
                                                   ),
                                                   hr(),
                                                   fluidRow(column(12, h4(strong(HTML('K<sub>out</sub>')), align = 'center')),
                                                            column(6, sliderInput("ic50kout_1", HTML('IC<sub>50</sub>'), 
                                                                                  min = 0, max = 200, value = 7.24, step=0.01,
                                                                                  ticks=F)),
                                                            column(6, sliderInput("tekout", HTML('T<sub>delay</sub> [hr]'), 
                                                                                  min = 0.0001, max = 168, value = 0.0001, 
                                                                                  step=1, ticks=F))
                                                   )
                                          )
                                   )
                            )
                          )
                  ),
                  tabItem(tabName = 'scheme',
                          fluidRow(
                            column(width=9, align = 'center',
                                   box(title = '', width = NULL, solidHeader = F, collapsed = F, 
                                       plotOutput("scheme", height = NULL)))))
                )
              )
)
