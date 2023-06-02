library(here)
library(tidyverse)

source(here("general_functions.R"))
params <- params %>% filter(params$it == 0)

ggfacet <- function(df, cols, nrow = 1) {
  df2 <- 
    df %>% 
    pivot_longer(all_of(cols), names_to = "name", values_to = "value") %>%
    mutate(name = factor(name, levels = cols))
  gg <- 
    ggplot(data = df2, aes(x = x, y = y, col = value)) + 
    geom_point() + 
    facet_wrap(~name, nrow = nrow) 
  return(gg)
}

plot_params <- function(i, params){
  genf <- get_data(i, params, "dos")
  gsdf <- get_data(i, params, "gsd")
  s <- sample(nrow(gsdf), 1000)
  gsd <- gsdf[s,]
  gen <- genf[s,]
  
  tp <- data.frame(name = names(params[i,]), val = unlist(params[i,])) 
  tc <- paste(paste0(tp$name, ": ", tp$val), collapse = " | ")

  df <- bind_cols(gsd, gen/2)
  df <- df[,c("env1", "z1", "0_0", "0_1", 
              "env2", "z2", "0_2", "0_3",
              "0_5", "0_6", "0_7", "0_8", 
              "x", "y")]
  plt <- 
    ggfacet(df, c("env1", "z1", "0_0", "0_1", 
                "env2", "z2", "0_2", "0_3",
                "0_5", "0_6", "0_7", "0_8"), nrow = 3) + 
    theme_bw() + 
    coord_equal() +
    scale_color_viridis_c(option = "viridis") + 
    ggtitle(tc) + 
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  
  print(plt)
  return(data.frame(params[i,], df))
}

pdf("example_plots.pdf")
dfs <- purrr::map(1:nrow(params), ~plot_params(.x, params))
dev.off()

dfs <- dfs %>% bind_rows()
write.csv(dfs, "example_df.csv", row.names = FALSE)