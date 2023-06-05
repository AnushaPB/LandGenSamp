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
    geom_point(alpha = 1, cex = 0.5) + 
    facet_wrap(~name, nrow = nrow) +
    theme_bw() + 
    coord_equal() +
    scale_color_viridis_c(option = "viridis") + 
    ggtitle(tc) + 
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
          panel.grid.major  = element_blank(), panel.grid.minor = element_blank())
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
  df <- df[,c("env1", "z1", "X0_0", "X0_1", "X0_2", "X0_3",
              "env2", "z2", "X0_4", "X0_5", "X0_6", "X0_7",
              "X0_8", "X0_9", "X0_10", "X0_11", "X0_12", "X0_13", 
              "x", "y")]
  plt <- 
    ggfacet(df, c("env1", "z1", "X0_0", "X0_1", "X0_2", "X0_3",
                "env2", "z2", "X0_4", "X0_5", "X0_6", "X0_7",
                "X0_8", "X0_9", "X0_10", "X0_11", "X0_12", "X0_13"), nrow = 3)
  
  print(plt)
  return(data.frame(params[i,], df))
}

pdf("example_plots.pdf")
dfs <- purrr::map(1:nrow(params), ~plot_params(.x, params))
dev.off()

dfs <- dfs %>% bind_rows()
write.csv(dfs, "example_df.csv", row.names = FALSE)

df <- read.csv(here("p3_methods", "example_df.csv"))
head(df)

df %>%
  filter(K == 1, phi == 0.5, m == 0.25, seed == 1, H == 0.05, r == 0.3) %>%
  ggfacet(c("env1", "z1", "X0_0", "X0_1", "X0_2", "X0_3",
            "env2", "z2", "X0_4", "X0_5", "X0_6", "X0_7",
            "X0_8", "X0_9", "X0_10", "X0_11", "X0_12", "X0_13"), nrow = 3)
