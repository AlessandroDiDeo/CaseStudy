function(input, output, session) {
  
  library(shiny)
  library(shinydashboard)
  library(deSolve)
  library(reshape2)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  
  source('ddw_kineto_models.R')
  
  mod1 <- reactive({
    
    if (input$mode_switch == 2){
      cond_CL	<- input$CL
      cond_V	<- input$V
      cond_KA	<- input$KA
    } else {
      cond_CL	<- 0
      cond_V	<- 1
      cond_KA	<- 1e+10
    }
    
    parameters = c(
      KPROT0      = 1.24 * 10^-3,
      KOUT0       = 4.63 * 10^-3,
      N           = 3.30,
      K1_1        = 6.64 * 10^-2,
      # K1_2        = 6.82E-02, 
      # K1_3        = 8.04E-02,
      GAMMA       = 1,
      GR1         = 7.98 * 10^-2,
      # GR2         = 5.29E-02,
      # GR3         = 5.11E-02, 
      NMAX        = 6.70 * 10^+2, 
      BETA        = 7.79 * 10^-1,
      KDEATH0     = 0.00016,
      
      IC50KIN     = input$ic50kin_1,
      IC50K1      = input$ic50k1_1,
      IC50KOUT    = input$ic50kout_1,
      IC50KD    	= input$ic50kd_1,
      EMAXKD      = input$emax_d,
      
      TEKD        = input$tekd,
      TEKIN       = input$tekin,
      TEK1        = input$tek1,
      TEKOUT      = input$tekout,
      
      CL          = cond_CL,
      V           = cond_V,
      KA          = cond_KA
      
    )
    
    S = c(
      S1 = (2.5*10^2) - ((2.5*10^2) * (47/100)),
      S2 = (2.5*10^2 * (47/100)),
      S3 = 0,
      S4 = 0,
      S5 = 0,
      S6 = 0,
      S7 = 0,
      S8 = (2.5*10^2) * (47/100) * 3.30,
      S9 = 0,
      CONC = 0
    )
    
    
    if (input$mode_switch == 2) {
      
      times = seq(0, input$interval_1[2]*24, by = 1)
      
      out1 = lsoda(y 		= S, 
                   times 	= times, 
                   func 	= ddw_kineto_mod1, 
                   parms 	= parameters,
                   rtol	= 1e-10,
                   events	= list(data=data.frame(var="S9",
                                                 time=seq(input$interval_1[1]*24,input$interval_1[2]*24, 
                                                          24/input$ndose_1), 
                                                 value=input$dose_1/input$molmass*1000/input$V, 
                                                 method=c("add"))))
    } else {
      
      times = seq(0, input$interval_1[2]*24, by = 1)
      
      out1 = lsoda(y 		= S, 
                   times 	= times, 
                   func 	= ddw_kineto_mod1, 
                   parms 	= parameters,
                   rtol	= 1e-10,
                   events	= list(data=data.frame(var="S9",
                                                 time=0, 
                                                 value=input$conc_1, 
                                                 method=c("rep"))))
    }
    
    
    return(out1)
  })
  
  output$plot_sim1 <- renderPlot({
    if (input$mode_switch == 1) {
      dat_tot = data.frame(mod1()[,c('S1', 'S2', 'S4', 'S5', 'S6', 'S7')])
      dat_tot_inf = data.frame(mod1()[,c('S2', 'S4', 'S5', 'S6', 'S7')])
      tot_cell = rowSums(dat_tot)
      tot_cell_inf = rowSums(dat_tot_inf)
      dat = data.frame(mod1()[,c(1, 9)], tot_cell, tot_cell_inf/tot_cell*100, mod1()[,c('CONC')],
                       stringsAsFactors=F)
      df = melt(dat, id.vars="time", variable.name="state")
      levels(df$state) <- c("Amastigotes/Cell", 'Total number of cells', 'Infected cells [%]', 
                            "Drug Concentration [uM]")
      df$value[df$state == 'Amastigotes/Cell'] <- df$value[df$state == 'Amastigotes/Cell']/tot_cell
      
      p = ggplot(data=df, aes(x=time, y=value, col=state)) + 
        geom_line() + 
        facet_wrap(~state, nrow=3, scales='free') +
        labs(x='Time [hr]', y='') +
        scale_color_discrete(name="") +
        geom_hline(data = data.frame(state = "Amastigotes/Cell", thresh = 1), aes(yintercept = thresh), linetype = 'dotted') +
        geom_hline(data = data.frame(state = "Drug Concentration [uM]", thresh = 38.4), aes(yintercept = thresh), linetype = 'dotted') +
        scale_x_continuous(breaks = seq(4, input$interval_1[2]*24, by = round(input$interval_1[2]*24/18))) +
        theme_bw() + 
        theme(
          axis.title.y = element_text(size=10, vjust=1),
          axis.title.x = element_text(size=10, vjust=1),
          axis.text.y = element_text(size=8),
          axis.text.x = element_text(size=8),
          strip.background = element_rect(fill='grey85',colour='grey85'),
          legend.justification=c(0,0), legend.position = 'bottom',
          legend.key = element_blank(),
          legend.key.size = unit(10, 'mm'))
    }
    else {
      dat_tot = data.frame(mod1()[,c('S1', 'S2', 'S4', 'S5', 'S6', 'S7')])
      dat_tot_inf = data.frame(mod1()[,c('S2', 'S4', 'S5', 'S6', 'S7')])
      tot_cell = rowSums(dat_tot)
      tot_cell_inf = rowSums(dat_tot_inf)
      dat = data.frame(mod1()[,c(1, 9)], tot_cell, tot_cell_inf/tot_cell*100, mod1()[,c('CONC')], 
                       stringsAsFactors=F)
      df = melt(dat, id.vars="time", variable.name="state")
      levels(df$state) <- c("Amastigotes/Cell", 'Total number of cells', 
                            'Infected cells [%]', "Drug Concentration [uM]")
      df$value[df$state == 'Amastigotes/Cell'] <- df$value[df$state == 'Amastigotes/Cell']/tot_cell
      
      p = ggplot(data=df, aes(x=time, y=value, fill = state), alpha = 0.3) + 
        geom_line(aes(col=state)) + 
        facet_wrap(~state, nrow=3, scales='free') +
        geom_hline(data = data.frame(state = factor("Amastigotes/Cell", levels = c("Amastigotes/Cell", 'Total number of cells', 
                                                                                   'Infected cells [%]', "Drug Concentration [uM]")), thresh = 1), 
                   aes(yintercept = thresh), linetype = 'dotted') +
        geom_hline(data = data.frame(state = "Drug Concentration [uM]", thresh = 38.4), aes(yintercept = thresh), linetype = 'dotted') +
        labs(x='Time [hr]', y='') +
        scale_color_discrete(name="") +
        scale_fill_discrete(name="") +
        scale_x_continuous(breaks = seq(4, (input$interval_1[2] + 3)*24, by = round((input$interval_1[2] + 3)*24/18))) +
        geom_ribbon(data = df[df$state != 'Infected cells [%]', ], 
                    aes(x = time, ymin = value/1.19, ymax = value*1.19), alpha = 0.3) +
        geom_ribbon(data = df[df$state == 'Infected cells [%]', ], 
                    aes(x = time, ymin = value/1.19, ymax = pmin(100, value*1.19)), alpha = 0.3) +
        theme_bw() + 
        theme(
          axis.title.y = element_text(size=10, vjust=1),
          axis.title.x = element_text(size=10, vjust=1),
          axis.text.y = element_text(size=8),
          axis.text.x = element_text(size=8),
          strip.background = element_rect(fill='grey85',colour='grey85'),
          legend.justification=c(0,0), legend.position = 'bottom',
          legend.key = element_blank(),
          legend.key.size = unit(10, 'mm'))
    }
    print(p)
  })
  
  output$plot_sim2 <- renderPlot({
    if (input$mode_switch == 1) {
      dat = data.frame(mod1()[,c(1:3, 5:8)], mod1()[,c(9)], mod1()[,c(4)], mod1()[,c('CONC')],
                       stringsAsFactors=F)
      df = melt(dat, id.vars="time", variable.name="state")
      levels(df$state) <- c("1. H9c2", '2. iH9c2', '3. S1', '4. S2', '5. S3', '6. S4', 
                            "7. Amastigotes", "8. Trypomastigotes", 'Drug Concentration [uM]')
      
      p = ggplot(data=df, aes(x=time, y=value, col=state)) + 
        geom_line() + 
        facet_wrap(~state, nrow=3, scales='free') +
        labs(x='Time [hr]', y='') +
        scale_color_discrete(name="") +
        scale_x_continuous(breaks = seq(4, input$interval_1[2]*24, by = round(input$interval_1[2]*24/18))) +
        geom_hline(data = data.frame(state = "Drug Concentration [uM]", thresh = 38.4), aes(yintercept = thresh), linetype = 'dotted') +
        theme_bw() + 
        theme(
          axis.title.y = element_text(size=10, vjust=1),
          axis.title.x = element_text(size=10, vjust=1),
          axis.text.y = element_text(size=8),
          axis.text.x = element_text(size=8),
          strip.background = element_rect(fill='grey85',colour='grey85'),
          legend.justification=c(0,0), legend.position = 'bottom',
          legend.key = element_blank(),
          legend.key.size = unit(10, 'mm'))
    }
    else {
      dat = data.frame(mod1()[,c(1:3, 5:8)], mod1()[,c(9)], mod1()[,c(4)], mod1()[,c('CONC')],
                       stringsAsFactors=F)
      df = melt(dat, id.vars="time", variable.name="state")
      levels(df$state) <- c("1. H9c2", '2. iH9c2', '3. S1', '4. S2', '5. S3', '6. S4', 
                            "7. Amastigotes", "8. Trypomastigotes", 'Drug Concentration [uM]')
      
      p = ggplot(data=df, aes(x=time, y=value)) + 
        geom_line(aes(col=state)) + 
        facet_wrap(~state, nrow=3, scales='free') +
        labs(x='Time [hr]', y='') +
        scale_color_discrete(name="") +
        scale_fill_discrete(name="") +
        geom_hline(data = data.frame(state = "Drug Concentration [uM]", thresh = 38.4), aes(yintercept = thresh), linetype = 'dotted') +
        scale_x_continuous(breaks = seq(4, (input$interval_1[2]+3)*24, by = round((input$interval_1[2]+3)*24/18))) +
        geom_ribbon(aes(x = time, ymin = value/1.19, ymax = value*1.19, fill = state), alpha = 0.3) +
        theme_bw() + 
        theme(
          axis.title.y = element_text(size=10, vjust=1),
          axis.title.x = element_text(size=10, vjust=1),
          axis.text.y = element_text(size=8),
          axis.text.x = element_text(size=8),
          strip.background = element_rect(fill='grey85',colour='grey85'),
          legend.justification=c(0,0), legend.position = 'bottom',
          legend.key = element_blank(),
          legend.key.size = unit(10, 'mm'))
    }
    print(p)
  })
  
  output$scheme <- renderImage({
    pixelRatio_curr <- session$clientData$pixelratio
    width_curr  <- session$clientData$output_scheme_width
    height_curr <- session$clientData$output_scheme_height
    return(list(
      src = "images/scheme.png",
      contentType = "image/png",
      alt = "Model scheme", heigth = height_curr, width = width_curr))
  }, deleteFile = F)
  
}