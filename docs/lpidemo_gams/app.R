
library(shiny)

### functions ####

# function to generate population time series
N <- function(x){
    popsize <- vector("numeric", 6)
    popsize[1] = 10
    for(i in 2:6){
        popsize[i] = popsize[i-1]*x
    }
    return(popsize)
}

# step 1: smooth the time series and predict full time series
library(mgcv)
smooth <- function(x) {
    time = 1:6
    # model each time series with a smoother (GAM)
    m <- gam(x ~ s(time, k = 3), method = "REML")
    # predict the gam over the time steps of the time series
    Npred <- predict.gam(m)
    return(Npred)
}

# function to calculate annual growth rate as ratio
dt <- function(x){
    dt <- vector("numeric", length(x))
    for(t in 2:length(x)){
        dt[t] <- x[t]/x[t-1]
    }
    dt = dt[-1]
    return(dt)
}

# function to calculate LPI from mean dt
calclpi <- function(dt){
    # calculate index value
    lpi = c(1) # initial value is 1 
    for(i in 2:length(dt)){
        lpi[i] <- lpi[i-1]*10^dt[i] }
    return(lpi)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("LPI demo"),
    
    fluidPage(
        fluidRow(
            column(4, 
                   p("Select the rate of change in the declining or growing populations below."),
                   sliderInput("change",
                               "Rate of change:",
                               min = 0,
                               max = 1,
                               value = 0.8, 
                               step = 0.1)
            ),
            column(4, 
                   h4("Population time series"),
                   p("Population time series are generated as:"),
                   withMathJax(helpText("$$N_t = rate*N_{t-1}$$"))
            ),
            column(4, plotOutput("abundances", height = '300px', width = '300px'))
        ),
        fluidRow(
            column(3, 
                   h4("Step 1: Smooth and predict time series"),
                   p("Population time series are smoothed with a GAM and predicted over the series' time interval. Each point is an observation of the abundance of the population at a given time step. Smooth trend is shown as a line overlayed on these data points."),
                   withMathJax(helpText("$$N_{pred} = s(time)$$"))
            ),
            column(3, 
                   h4("Step 2: Calculate annual growth rates of each population"),
                   p("Population growth rates are calculated as the ratio of abundances at each time step compared to the previous time step:"),
                   withMathJax(helpText("$$dt = N_{pred}/{N_{pred(t-1)}}$$"))
            ),
            column(3, 
                   h4("Step 3: Get annual (geometric) mean growth rate"),
                   p("For each time step, we then take the geometric mean of the annual growth rates of the three populations:"),
                   withMathJax(helpText("$$\\overline{dt} = (dt_A * dt_B * dt_C)^{1/3}$$"))
            ),
            column(3, 
                   h4("Step 4: Calculate index"),
                   p("We then calculate the index, which converts these growth rates relative to baseline value of 1 at the first time step."),
                   withMathJax(helpText("$$I_t = I_{t-1}*10^{\\overline{dt}}$$"))
            )
        ),
        fluidRow(
            column(3, plotOutput("gams")),
            column(3, plotOutput("growthrates")),
            column(3, plotOutput("meangrowthrate")),
            column(3, plotOutput("lpi"))
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    time <- 1:6
    
    # make populations
    A <- reactive(N(1/input$change))
    B <- N(1)
    C <- reactive(N(input$change))
    
    output$abundances <- renderPlot({
        
        # plot them
        plot(time, A(), type = "l", ylim = c(0, max(A())), ylab = "Abundance", xlab = "time")
        lines(B)
        lines(C())
        
    })
    
    output$gams <- renderPlot({
        
        # plot the smoothed gams
        plot(time, smooth(A()), ylim = c(0, max(smooth(A()))), 
             ylab = "Smoothed abundance", xlab = "time", type = "l")
        lines(B)
        lines(smooth(C()))
        
    })
    
    output$growthrates <- renderPlot({
        
        # calculate growth rates
        dtA <- dt(A()) 
        dtB <- dt(B)
        dtC <- dt(C())
        
        # plot them
        plot(time[-1], dtA, type = "l", ylim = c(0,max(dtA)), 
             ylab = "Growth rate", xlab = "time")
        lines(x = time[-1], y = dtB)
        lines(x = time[-1], dtC)
        
    })
    
    output$meangrowthrate <- renderPlot({
        
        # calculate growth rates
        dtA <- dt(A()) 
        dtB <- dt(B)
        dtC <- dt(C())
        
        # calculate mean growth rate
        df <- cbind(dtA, dtB, dtC)
        means <- apply(df, 1, function(x) prod(x, na.rm = TRUE)^(1/3))
        df <- as.data.frame(cbind(df, means))
        
        # plot them
        plot(time[-1], df$means, type = "l", ylab = "Growth rate", xlab = "time")
        
    })
    
    output$lpi <- renderPlot({
        
        # calculate growth rates
        dtA <- dt(A()) 
        dtB <- dt(B)
        dtC <- dt(C())
        
        # calculate mean growth rate
        df <- cbind(dtA, dtB, dtC)
        means <- apply(df, 1, function(x) prod(x, na.rm = TRUE)^(1/3))
        df <- as.data.frame(cbind(df, means))
        
        lpi <- calclpi(log10(df$means))
        plot(lpi ~ time[-1], type = "l")
        
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
