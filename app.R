library(shiny)
library(tidyverse)
library(here)
source(here("p4_analysis", "analysis_functions.R"))


# Functions for plotting ---------------------------------------------------------------------------------------------
get_data <- function(method, sampling) {
  if (sampling == "individual") sampling <- "indsampling"
  if (sampling == "site") sampling <- "sitesampling"
  
  file_name <- paste0(method, "_", sampling, "_results.csv")
  file_path <- here("p3_methods", "outputs", file_name)
  
  if (method == "lfmm") df <- format_lfmm(file_path)
  if (method == "rda") df <- format_rda(file_path)
  if (method == "mmrr") df <- format_mmrr(file_path)
  if (method == "gdm") df <- format_gdm(file_path)
  
  return(df)
}

map2(expand_grid(list("individual", "site")))

MEGAPLOT <- function(df, stat_name, minv = NULL, maxv = NULL, aggfunc = mean, colpal = "plasma", direction = 1, divergent = FALSE, na.rm=TRUE, dig = 3){
  
  agg <- 
    df %>%
    # convert to data.frame
    data.frame() %>%
    # remove full data
    filter(sampstrat != "full") %>%
    # rename stat_name column as stat for simplified use
    rename("stat" = all_of(stat_name)) %>%
    # group by all params
    group_by(K, phi, m, H, r, nsamp, sampstrat) %>%
    # summarize
    custom_agg(aggfunc, na.rm) %>%
    # create new group for plotting
    mutate(group = paste0("K=", K,
                          " phi=", phi,
                          " m=", m,
                          "\nH=", H,
                          " r=", r))
  
  p <- ggplot(agg, aes(nsamp, sampstrat)) +
    geom_tile(aes(fill = stat)) + 
    geom_text(aes(label = round(stat, dig), hjust = 0.5)) +
    theme_bw() +
    coord_fixed() + 
    facet_wrap( ~ group, nrow = 4) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), legend.position = "none",
          axis.title.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.ticks.y=element_blank(),
          axis.text.x = element_text(color = "grey50", size = 14),
          axis.text.y = element_text(color = "gray50", size = 14), 
          plot.margin=unit(rep(0.4,4),"cm"),
          strip.text.x = element_text(size = 18),
          strip.text.y = element_text(size = 18),
          strip.background = element_blank(),
          strip.text = element_text(color = "black")) 
  
  if(divergent){
    p <- p + scale_fill_gradient2(low = "#2066AC", mid="#F7F7F7", high = "#d79232", midpoint = 0, limits=c(minv, maxv))
  } else {
    p <- p + scale_fill_viridis(limits=c(minv, maxv), option = colpal, direction = direction)
  }
  
  return(p)
  
}


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Hello Shiny!"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "bins",
                  label = "Number of bins:",
                  min = 1,
                  max = 50,
                  value = 30)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      plotOutput(outputId = "distPlot")
      
    )
  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlot({
    
    x    <- faithful$waiting
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    hist(x, breaks = bins, col = "#007bc2", border = "white",
         xlab = "Waiting time to next eruption (in mins)",
         main = "Histogram of waiting times")
    
  })
  
}


shinyApp(ui = ui, server = server)