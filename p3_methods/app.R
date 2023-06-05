library(shiny)
library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)

ggfacet <- function(df, cols, nrow = 1) {
  df2 <- 
    df %>% 
    pivot_longer(all_of(cols), names_to = "name", values_to = "value") %>%
    mutate(name = factor(name, levels = cols))
  gg <- 
    ggplot(data = df2, aes(x = x, y = y, col = value)) + 
    geom_point(alpha = 1, cex = 2) + 
    facet_wrap(~name, nrow = nrow) +
    theme_bw() + 
    coord_equal() +
    scale_color_viridis_c(option = "viridis") + 
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
          panel.grid.major  = element_blank(), panel.grid.minor = element_blank())
  return(gg)
}


# Define your data frame here
df <- read.csv(here("p3_methods", "example_df.csv"))

# Define the UI for the Shiny app
ui <- fluidPage(
  titlePanel("Plot Generator"),
  sidebarLayout(
    sidebarPanel(
      selectInput("K", "K", choices = c(1, 2), selected = 2),
      selectInput("phi", "phi", choices = c(0.5, 1.0), selected = 1.0),
      selectInput("m", "m", choices = c(0.25, 1), selected = 1.0),
      selectInput("r", "r", choices = c(0.3, 0.6), selected = 0.6),
      selectInput("H", "H", choices = c(0.05, 0.5), selected = 0.5),
      selectInput("seed", "seed", choices = c(1, 2, 3), selected = 1)
    ),
    mainPanel(
      plotOutput("plot", width = "1200px", height = "600px")
    )
  )
)

# Define the server logic
server <- function(input, output) {
  output$plot <- renderPlot({
    filtered_df <- df %>%
      filter(K == input$K, phi == input$phi, m == input$m, seed == input$seed, H == input$H, r == input$r)
    
    g <- ggfacet(filtered_df, c("env1", "z1", "X0_0", "X0_1", "X0_2", "X0_3",
                                "env2", "z2", "X0_4", "X0_5", "X0_6", "X0_7",
                                "X0_8", "X0_9", "X0_10", "X0_11", "X0_12", "X0_13"), nrow = 3)
    
    print(g)
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

