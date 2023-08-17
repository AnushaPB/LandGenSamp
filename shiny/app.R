library("tidyverse")
library("here")
library("shiny")
library("viridis")
source(here("p4_analysis", "analysis_functions.R"))
source(here("shiny", "app_functions.R"))

# Load data
combos <- expand_grid(method = c("lfmm", "rda", "gdm", "mmrr", "mmrr2", "gdm2"), sampling = c("individual", "site"))
clean_data <- pmap(combos, get_app_data)
names(clean_data) <- pmap_chr(combos, ~paste0(.x,"_", .y))

# Make key for stats
stat_options <- 
  list(lfmm = colnames(clean_data$lfmm_ind)[-c(1:9)],
       rda = colnames(clean_data$rda_ind)[-c(1:9)],
       mmrr = colnames(clean_data$mmrr_ind)[-c(1:9)],
       gdm = colnames(clean_data$gdm_ind)[-c(1:9)],
       mmrr2 = colnames(clean_data$mmrr2_ind)[-c(1:9)],
       gdm2 = colnames(clean_data$gdm2_ind)[-c(1:9)])


stat_convert <- function(stat){
  if (stat == "nloci detected") return("TOTALN")
  if (stat == "number of latent factors") return("K.1")
  if (stat == "Ratio Absolute Error") return("RAE")
  if (stat == "Ratio Bias") return("ratio_err")
  if (stat == "IBD Absolute Error") return("geo_ae")
  if (stat == "IBD Bias") return("geo_err")
  if (stat == "IBE Absolute Error") return("env_ae")
  if (stat == "IBE Bias") return("comboenv_err")
  if (stat == "Ratio of IBE to IBD") return("ratio")
  if (stat == "IBD Coefficient") return("geo_coeff")
  if (stat == "IBE Coefficient") return("comboenv_coeff")
  if (stat == "Proportion of Null Models") return("ratio")
  return(stat)
}



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
                  label = "method",
                  choices = c("lfmm", "rda", "mmrr", "mmrr2", "gdm", "gdm2"),
                  selected = "lfmm"), 
      
      # Select input for sampling type
      selectInput("sampling",
                  label = "sampling type",
                  choices = c("individual", "site"),
                  selected = "individual"),
      
      # Select whether to summarize
      selectInput("plot_type",
                   label = "plot",
                   choices = c("megaplot", "summary"),
                   selected = "megaplot"),
      
      # Select input for the statistic
      uiOutput("statSelector"),
      
      # Additional selection boxes based on method
      uiOutput("additionalSelections"),
      
      # Additional selection boxes for K, m, phi, r, H
      uiOutput("additionalVariables"),
      
      width = 3
      
    ),
    
    
    # Main panel for displaying outputs ----
    # Main panel for displaying outputs ----
    mainPanel(
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      ), 
      
      plotOutput(outputId = "distPlot", width = "100%")
    
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
        selectInput("lfmm_method",
                    "lfmm method",
                    choices = c("ridge", "lasso"),
                    selected = "ridge"),
        selectInput("K_method",
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
                    "P adjust Method",
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
      "statistic",
      choices = stat_options[[method]],
      selected = stat_options[[method]][1]
    )
  })
  
  # Render additional variables selection based on the selected method
  output$additionalVariables <- renderUI({
    tagList(
      checkboxGroupInput("K", "population size (K)", choices = c("low", "high"), selected = c("low", "high")),
      checkboxGroupInput("m", "migration (m)", choices = c("low", "high"), selected = c("low", "high")),
      checkboxGroupInput("phi", "selection strength (phi)", choices = c("low", "high"), selected = c("low", "high")),
      checkboxGroupInput("r", "environmental correlation (r)", choices = c("low", "high"), selected = c("low", "high")),
      checkboxGroupInput("H", "environmental autocorrelation (H)", choices = c("low", "high"), selected = c("low", "high"))
    )
  })
  
  # Render the plot based on the selected method and statistic
  output$distPlot <- renderPlot({
    x <- clean_data[[paste0(input$method, "_", input$sampling)]]
    x <- filter_additional_variables(x, input)
    stat_name <- stat_convert(input$stat)
    
    # Filter data based on selected method and additional variables
    if (input$method == "lfmm"){
      x <- 
        x %>%
        filter(lfmm_method == input$lfmm_method) %>%
        filter(K_method == input$K_method)
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
    
    colpal <- "cividis"
    error_stat <- grepl("*_ae*", stat_name) | grepl("*err*", stat_name)
    if (grepl("TPR", stat_name)) colpal <- "plasma"
    if (grepl("FDR", stat_name) | error_stat) {colpal <- "viridis"; direction <- -1} else direction <- 1
    if (stat_name %in% c("TOTALN", "geo_coeff", "comboenv_coeff", "ratio", "K_factor")) colpal <- "cividis"
    if (stat_name %in% c("ratio_err", "geo_coeff_err", "comboenv_coeff_err")) divergent <- TRUE else divergent <- FALSE
    
    if (input$plot_type == "megaplot") print(MEGAPLOT(x, stat_name = stat_name, dig = 2, colpal = colpal, divergent = divergent, direction = direction))
    if (input$plot_type == "summary") print(summary_hplot(x, stat_name = stat_name, dig = 2, colpal = colpal, divergent = divergent, direction = direction))
    
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
    phi = c(0.5, 1),
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

    