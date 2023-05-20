# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Landscape genomic simulation results"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Select method
      selectInput("method", 
                  label = "Method",
                  choices = c("lfmm", "rda", "mmrr", "gdm"),
                  selected = "lfmm"), 
      
      # Select input for sampling type
      selectInput("sampling",
                  label = "Sampling Type",
                  choices = c("individual", "site"),
                  selected = "individual"),
      
      # Select input for the statistic
      uiOutput("statSelector"),
      
      # Additional selection boxes based on method
      uiOutput("additionalSelections"),
      
      # Additional selection boxes for K, m, phi, r, H
      uiOutput("additionalVariables"),
      
      width = 3
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "distPlot", width = 9)
      
    )
  )
)
# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Render additional selection boxes based on the selected method
  output$additionalSelections <- renderUI({
    method <- input$method
    
    if (method == "lfmm") {
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
    } else if (method == "rda") {
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
  
  # Render the statistic selector based on the selected method
  output$statSelector <- renderUI({
    method <- input$method
    selectInput(
      "stat",
      "Statistic",
      choices = stat_options[[method]],
      selected = stat_options[[method]][1]
    )
  })
  
  # Render additional variables selection based on the selected method
  output$additionalVariables <- renderUI({
    tagList(
      checkboxGroupInput("K", "K", choices = c("low", "high"), selected = c("low", "high")),
      checkboxGroupInput("m", "m", choices = c("low", "high"), selected = c("low", "high")),
      checkboxGroupInput("phi", "phi", choices = c("low", "high"), selected = c("low", "high")),
      checkboxGroupInput("r", "r", choices = c("low", "high"), selected = c("low", "high")),
      checkboxGroupInput("H", "H", choices = c("low", "high"), selected = c("low", "high"))
    )
  })
  
  # Render the plot based on the selected method and statistic
  output$distPlot <- renderPlot({
    x <- get_data(input$method, input$sampling)
    x <- filter_additional_variables(x, input)
    stat_name <- input$stat
    
    # Filter data based on selected method and additional variables
    if (input$method == "lfmm"){
      x <- 
        x %>%
        filter(method == input$lfmmMethod) %>%
        filter(K_selection == input$K_selection)
    }
    
    if (input$method == "rda"){
      input_logic <- input$correctPC == "yes"
      x <- 
        x %>% filter(correctPC == input_logic) 
    }
    
    if (input$method == "rda" | input$method == "lfmm"){
      x <- 
        x %>%
        filter(padj == input$padj) %>%
        filter(sig == input$sig)
    }
    
    MEGAPLOT(x, stat_name = stat_name, dig = 2)
  }, height = 945*.95, width = 1525.5*.95)
  
}

# Function to filter data based on additional variables
filter_additional_variables <- function(filtered_data, input) {
  # Helper function to filter based on low/high selection
  filter_low_high <- function(df, var_name, low_val, high_val) {
    if ("low" %in% input[[var_name]] && "high" %in% input[[var_name]]) {
      # Keep both low and high values
      return(df)
    } else if ("low" %in% input[[var_name]]) {
      return(df[df[[var_name]] == low_val, ])
    } else if ("high" %in% input[[var_name]]) {
      return(df[df[[var_name]] == high_val, ])
    }
  }
  
  # Define the filter ranges for each variable
  filter_ranges <- list(
    K = c(1, 2),
    m = c(0.25, 1.00),
    phi = c(0.1, 0.5),
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

    