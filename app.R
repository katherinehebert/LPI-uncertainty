#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(dplyr)
library(tidyr)
library(shiny)
library(ggplot2)
library(patchwork)
library(thematic)
source('~/Documents/GitHub/LPI-Sensitivity/scripts/sim_mech.R')

plot_function <- function(df){
    ggplot(df) + 
        geom_line(aes(x = time, y = N, col = popID, group = popID)) + 
        coord_cartesian(ylim = c(10, 190)) +
        facet_wrap(~set, dir = "h") +
        theme(legend.position = "none") 
}

# styling
library(bslib)
cute_theme <- bs_theme(
    bg = "#FFFFFF", fg = "#003f5c", primary = "#003f5c", 
    base_font = font_google("Roboto"),
    heading_font = font_google("Roboto")
)

## CARRYING CAPACITY SCENARIOS ##  ------
K_increase = 100 + 5*c(0:9)
K_stable = rep(100, 10)
K_decline = 100 - 5*c(0:9)
K_list = list(K_decline, K_stable, K_increase)

# Define UI for application that draws a histogram
ui <- fluidPage(theme = cute_theme,

    # Application title
    titlePanel("Simulating interacting populations"),

     # Sidebar with a slider input for number of bins 
     sidebarLayout(
        sidebarPanel(
            sliderInput("lambda_i",
                        "Growth rate (populations i):",
                        min = 0.5,
                        max = 2,
                        value = 1.5),
            sliderInput("lambda_j",
                        "Growth rate (populations j):",
                        min = 0.5,
                        max = 2,
                        value = 1.5),
            sliderInput("alpha_ji",
                        "Interaction effect of J on I:",
                        min = -0.5,
                        max = 0.5,
                        value = -0.01),
            sliderInput("alpha_ij",
                        "Interaction effect of I on J:",
                        min = -0.5,
                        max = 0.5,
                        value = -0.01),
            sliderInput("process",
                        "Process error:",
                        min = 0,
                        max = 1,
                        value = 0.1),
            sliderInput("obs",
                        "Observation error:",
                        min = 0,
                        max = 50,
                        value = 10),
            sliderInput("lag",
                        "Lag in interaction:",
                        min = 0,
                        max = 5,
                        value = 0)
         ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("populationPlot", height = "700px")
        )
    )
)

thematic_shiny()

# Define server logic required to draw a histogram
server <- function(input, output) {
    

    
    output$populationPlot <- renderPlot({
        
        sim_list <- list()
        for(i in 1:3){
            sim_list[[i]] <- sim_mech(
            n_pairs = 10, timesteps = 10,
            N0i = 100, N0j = 100,
            lambda_i = input$lambda_i, 
            lambda_j = input$lambda_j,
            alpha_ij = input$alpha_ij, 
            alpha_ji = input$alpha_ji,
            process = input$process, 
            observation = input$obs,
            K = K_list[[i]],
            lag_value = input$lag,
            save_figs = FALSE
        )
}
        A <- plot_function(sim_list[[1]]) + labs(title = "Decline")
        B <- plot_function(sim_list[[2]]) + labs(title = "Stable")
        C <- plot_function(sim_list[[3]]) + labs(title = "Increase")
        A / B / C + plot_annotation(tag_levels = "a")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
