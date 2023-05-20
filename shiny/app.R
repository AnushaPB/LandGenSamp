library(shiny)
library(gridExtra)
library(here)
library(tidyverse)
library(viridis)

# Define UI for app that draws histograms ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Landscape genomic simulation results"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Select methods
      selectInput("methods", 
                  label = "Methods",
                  choices = c("lfmm", "rda", "mmrr", "gdm"),
                  selected = "lfmm",
                  multiple = TRUE), 
      
      # Select input for sampling type
      selectInput("sampling",
                  label = "Sampling Type",
                  choices = c("individual", "site"),
                  selected = "individual"),
      
      # Select input for the statistic
      uiOutput("statSelector"),
      
      # Additional selection boxes based on methods
      uiOutput("additionalSelections"),
      
      # Additional selection boxes for K, m, phi, r, H
      uiOutput("additionalVariables"),
      
      width = 3
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Plots side by side ----
      plotOutput(outputId = "distPlots", width = "100%")
      
    )
  )
)

# Define server logic required to draw histograms ----
server <- function(input, output) {
  
  # Render additional selection boxes based on the selected methods
  output$additionalSelections <- renderUI({
    methods <- input$methods
    
    if ("lfmm" %in% methods) {
      tagList(
        selectInput("lfmmMethod",
                    "LFMM Method",
                    choices = c("ridge", "lasso"),
                    selected = "ridge"),
        selectInput("K_selection",
                    "K selection method",
                    choices = c("tess", "tracy.widom", "find.clusters", "quick.elbow"),
                    selected = "tess"),
        selectInput("padj",
                    "p adjust method",
                    choices = c("none", "fdr", "holm", "bonferroni"),
                    selected = "fdr"),
        selectInput("sig",
                    "alpha level",
                    choices = c("0.05", "0.1"),
                    selected = "0.05")
      )
    } else if ("rda" %in% methods) {
      tagList(
        selectInput("padj",
                    "P Adjust Method",
                    choices = c("none", "fdr", "holm", "bonferroni"),
                    selected = "fdr"),
        selectInput("correctPC",
                    "PCA correction",
                    choices = c("yes", "no"),
                    selected = "no")
      )
    } else {
      NULL
    }
  })
  
  # Render the statistic selector based on the selected methods
  output$statSelector <- renderUI({
    methods <- input$methods
    selectInput(
      "stat",
      "Statistic",
      choices = stat_options[[methods[1]]],  # Use the first selected method for choices
      selected = stat_options[[methods[1]]][1]  # Use the first selected method for default
    )
  })
  
  # Render additional variables selection based on the selected methods
  output$additionalVariables <- renderUI({
    tagList(
      checkboxGroupInput("K", "K", choices = c("low", "high"), selected = c("low", "high")),
      checkboxGroupInput("m", "m", choices = c("low", "high"), selected = c("low", "high")),
      checkboxGroupInput("phi", "phi", choices = c("low", "high"), selected = c("low", "high")),
      checkboxGroupInput("r", "r", choices = c("low", "high"), selected = c("low", "high")),
      checkboxGroupInput("H", "H", choices = c("low", "high"), selected = c("low", "high"))
    )
  })
  
  # Render the plots based on the selected methods and statistic
  # Render the plots based on the selected methods and statistic
  output$distPlots <- renderPlot({
    methods <- input$methods
    stat_name <- input$stat
    
    # Create an empty list to store the plots
    plots <- list()
    
    # Loop through each selected method
    for (method in methods) {
      x <- get_data(method, input$sampling)
      x <- filter_additional_variables(x, input, method)
      
      # Filter data based on selected method and additional variables
      if (method == "lfmm") {
        if (!is.null(input$lfmmMethod) && !is.null(input$K_selection)) {
          x <- 
            x %>%
            filter(method == input$lfmmMethod) %>%
            filter(K_selection == input$K_selection)
        }
      }
      
      if (method == "rda") {
        if (!is.null(input$correctPC)) {
          input_logic <- input$correctPC == "yes"
          x <- 
            x %>% filter(correctPC == input_logic) 
        }
      }
      
      if (method == "rda" || method == "lfmm") {
        if (!is.null(input$padj) && !is.null(input$sig)) {
          x <- 
            x %>%
            filter(padj == input$padj) %>%
            filter(sig == input$sig)
        }
      }
      
      # Generate the plot for the current method
      plot <- MEGAPLOT(x, stat_name = stat_name, dig = 2)
      
      # Add the plot to the list
      plots[[method]] <- plot
    }
    
    # Arrange the plots side by side
    grid.arrange(grobs = plots, ncol = length(plots))
  }, height = 945 * 0.95, width = 1525.5 * 0.95)
  
  
}

# Function to filter data based on additional variables
filter_additional_variables <- function(filtered_data, input, method) {
  # Helper function to filter based on low/high selection
  filter_low_high <- function(df, var_name, low_val, high_val) {
    if ("low" %in% input[[var_name]] && "high" %in% input[[var_name]]) {
      # Keep both low and high values
      return(df)
    } else if ("low" %in% input[[var_name]]) {
      return(df[df[[var_name]] == low_val, ])
    } else if ("high" %in% input[[var_name]]) {
      return(df[df[[var_name]] == high_val, ])
    } else {
      return(df)
    }
  }
  
  # Define the filter ranges for each variable
  filter_ranges <- list(
    K = c(1, 2),
    m = c(0.25, 1.00),
    phi= c(0.1, 0.5),
    r = c(0.3, 0.6),
    H = c(0.05, 0.5)
  )
  
  # Apply filters using purrr::map
  filtered_data <- reduce(names(filter_ranges), function(data, var) {
    filter_low_high(data, var, filter_ranges[[var]][1], filter_ranges[[var]][2])
  }, .init = filtered_data)
  
  return(filtered_data)
}

shinyApp(ui = ui, server = server)

