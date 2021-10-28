
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
            column(12, p("Select the rate of change in the declining or growing populations below.")
                   ),
        ),
        fluidRow(
            column(6,
                   sliderInput("change",
                               "Rate of change:",
                               min = 0,
                               max = 1,
                               value = 0.5, 
                               step = 0.1)
            )
            ),
        fluidRow(
            column(3, 
                   h4("Population abundance time series"), 
                   br(),
                   p("Say we have 3 populations (A, B, C): one that has grown, one that has remained stable, and one that has declined."),
                   plotOutput("abundances")
            ),
            column(3, 
                   h4("Step 1: Calculate annual growth rates of each population"),
                   p("Population growth rates are calculated as the ratio of abundances at each time step compared to the previous time step:"),
                   withMathJax(helpText("$$dt = N/{N_{t-1}}$$")),
                   plotOutput("growthrates")
            ),
            column(3, 
                   h4("Step 2: Get annual (geometric) mean growth rate"),
                   p("For each time step, we then take the geometric mean of the annual growth rates of the three populations:"),
                   withMathJax(helpText("$$\\overline{dt} = (dt_A * dt_B * dt_C)^{1/3}$$")),
                   plotOutput("meangrowthrate")
            ),
            column(3, 
                   h4("Step 3: Calculate index"),
                   p("We then calculate the index, which converts these growth rates relative to baseline value of 1 at the first time step."),
                   withMathJax(helpText("$$I_t = I_{t-1}*10^{\\overline{dt}}$$")),
                   plotOutput("lpi")
                   )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    time = 1:6
    
    # make populations
    A <- reactive(N(1/input$change))
    B <- N(1)
    C <- reactive(N(input$change*1))
    
    output$abundances <- renderPlot({
        
        # plot them
        plot(time, A(), type = "l", ylim = c(0, max(A())), ylab = "Abundance", xlab = "time")
        lines(B)
        lines(C())
        
    })
    
    output$growthrates <- renderPlot({
        
        # calculate growth rates
        dtA <- dt(A()) 
        dtB <- dt(B)
        dtC <- dt(C())
        
        # plot them
        plot(time[-1], dtA, type = "l", ylim = c(0,max(dtA)), 
             ylab = "Growth rate", xlab = "time")
        lines(x= time[-1], y = dtB)
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
        plot(time[-1], df$means, type = "l", ylim = c(0,2), 
             ylab = "Growth rate", xlab = "time")
        
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
        plot(lpi ~ time[-1], type = "l", ylim = c(0,2)
             )
        
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
