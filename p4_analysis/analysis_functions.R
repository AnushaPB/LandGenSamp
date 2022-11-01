MEGAPLOT <- function(df, stat_name, minv = NULL, maxv = NULL, aggfunc = "mean", colpal = "plasma", direction = 1, divergent = FALSE, na.rm=TRUE){
  
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
    mutate(group = paste0("K=", agg$K,
                       " phi=", agg$phi,
                       " m=", agg$m,
                       "\nH=", agg$H,
                       " r=", agg$r))
  
  p <- ggplot(agg, aes(nsamp, sampstrat)) +
    geom_tile(aes(fill = stat)) + 
    geom_text(aes(label = signif(stat, digits = sigdig), hjust = 0.5)) +
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
          strip.background =element_blank(),
          strip.text = element_text(color = "black")) 
  
  if(divergent){
    p <- p + scale_fill_gradient2(low = "#2066AC", mid="#F7F7F7", high = "#d79232", midpoint = 0, limits=c(minv, maxv))
  } else {
    p <- p + scale_fill_viridis(limits=c(minv, maxv), option = colpal, direction = direction)
  }

  return(p)
  
}

summary_hplot <- function(df, stat_name = "stat", na.rm = TRUE, colpal = "plasma", sigdig=2, aggfunc = "mean", minv = NULL, maxv = NULL, direction = 1, divergent = FALSE){
 
  # Summarize dataframe
  agg <- df %>% 
    # convert to data.frame
    data.frame() %>%
    # remove full data
    filter(sampstrat != "full") %>%
    # rename stat_name column as stat for simplified use
    rename("stat" = all_of(stat_name)) %>%
    # select columns
    select(K, phi, m, H, r, nsamp, sampstrat, stat) %>%
    # convert from wide to long format
    pivot_longer(c(K, phi, m, H, r)) %>%
    # turn value into numeric from factor
    mutate(value = as.numeric(as.character(value))) %>%
    # group by the parameter (name)
    group_by(name) %>% 
    # create new low/high variable (the group_by name makes sure this is done for each parameter seperately)
    mutate(low_high = case_when(
        value == min(value) ~ "low",
        value == max(value) ~ "high",
        TRUE ~ "NA")) %>% 
    # group again by the new group (low_high), the param (name), and nsamp/sampstrat
    group_by(name, low_high, nsamp, sampstrat) %>%
    # summarize by group
    custom_agg(aggfunc, na.rm) %>%
    # order nsamp for plotting
    mutate(nsamp = factor(nsamp, levels = unique(nsamp)[order(unique(nsamp))])) %>%
    # order low_high for plotting
    mutate(low_high = factor(low_high, ordered = TRUE, levels = c("low", "high"))) %>%
    # create friendly names
    mutate(param = case_when(name == 'K' ~ 'population size',
                             name == 'm' ~ 'migration',
                             name == 'phi' ~ "selection strength",
                             name == "H" ~ "spatial autocorrelation",
                             name == "r" ~ "correlation", 
                             TRUE ~ "NA")) %>%
    # order param for plotting
    mutate(param = factor(param, ordered = TRUE, levels = c("population size", 
                                                   "migration", 
                                                   "selection strength", 
                                                   "spatial autocorrelation", 
                                                   "correlation")))
  
  # define max and min for plotting
  if(is.null(maxv)){ maxv <- max(agg$stat, na.rm = TRUE)}
  if(is.null(minv)){ minv <- min(agg$stat, na.rm = TRUE)}
  
  # plot results
  p <- ggplot(agg, aes(nsamp, sampstrat)) +
    geom_tile(aes(fill = stat)) + 
    geom_text(aes(label = signif(stat, digits = sigdig), hjust = 0.5)) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), legend.position = "none",
          axis.title.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.ticks.y=element_blank(),
          axis.text.x = element_text(color = "grey50", size = 14),
          axis.text.y = element_text(color = "gray50", size = 14), 
          plot.margin=unit(rep(0.4,4),"cm"),
          strip.text.x = element_text(size = 18),
          strip.text.y = element_text(size = 18),
          strip.background =element_blank(),
          strip.text = element_text(color = "black")) +
    coord_fixed() + 
    facet_grid(low_high ~ param)
  
  if(divergent){
    p <- p + scale_fill_gradient2(low = "#2066AC", mid="#F7F7F7", high = "#d79232", midpoint = 0, limits=c(minv, maxv))
  } else {
    p <- p + scale_fill_viridis(limits = c(minv, maxv), option = colpal, direction = direction)
  }
  
  return(p)
}



summary_vplot <- function(df, stat = "stat", plot.type = "rain_plot", varlist = c("K", "phi", "H", "r", "m", "sampstrat"), colpal = "plasma", nrow = 2){
  
  df$stat <- df[,stat]
  
  plot_list <- map(varlist, get(plot.type), df, stat, colpal)
  
  pl <- do.call("grid.arrange", c(plot_list, nrow = nrow))
  
  return(pl)
}


vplot <- function(var, df, stat, colpal = "plasma"){
  p <- ggplot(df, aes(fill=get(var), y=stat, x=factor(nsamp))) + 
    #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
    scale_fill_viridis(discrete=T, option=colpal, name = var) +
    xlab("nsamp") +
    theme_bw()+
    geom_boxplot(width = 0.4, position = position_dodge(width = 0.9))+
    #geom_point(aes(col=K),position = position_dodge(width = 0.9), alpha=0.1)+
    scale_colour_viridis(discrete=T, option=colpal, name = var)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=21))
  
  return(p)
}


rain_plot <- function(var, df, stat, colpal = "plasma"){
  p <- ggplot(df, aes(fill=get(var), y=stat, x=factor(nsamp))) + 
    ## add half-violin from {ggdist} package
    ggdist::stat_halfeye(
      ## custom bandwidth
      adjust = .5, 
      ## adjust height
      width = .6, 
      ## move geom to the right
      justification = -.2, 
      ## remove slab interval
      .width = 0, 
      point_colour = NA,
      position = position_dodge(width = 0.9)
    ) + 
    geom_boxplot(
      width = .12, 
      ## remove outliers
      outlier.color = NA, ## `outlier.shape = NA` works as well
      position = position_dodge(width = 0.9)
    ) +
    scale_fill_viridis(discrete=T, 
                       option=colpal, 
                       name = var) +
    
    scale_colour_viridis(discrete=T, 
                         option=colpal, 
                         name = var) +
    xlab("nsamp") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=21)) + 
    ## remove white space on the left
    coord_cartesian(xlim = c(1.2, NA))
  
  return(p)
}

custom_agg <- function(agg, aggfunc = "mean", na.rm = TRUE){
  if(aggfunc == "mean") agg <- agg %>% summarize(stat = mean(stat, na.rm = na.rm), .groups = "keep")
  
  if(aggfunc == "count_na") agg <- agg %>% summarize(stat = sum(is.na(stat)), .groups = "keep")
  
  if(aggfunc == "count_null") agg <- agg %>% summarize(stat = sum(is.na(stat)), .groups = "keep")
  
  if(aggfunc == "prop_na") agg <- agg %>% summarize(stat = mean(is.na(stat)), .groups = "keep")
  
  if(aggfunc == "var") agg <- agg %>% summarize(stat = var(stat, na.rm = na.rm), .groups = "keep")
  
  if(aggfunc == "rmse") {
    agg <- moddf %>%
      mutate(stat = stat^2) %>%
      group_by(K, phi, m, H, r, nsamp, sampstrat) %>%
      summarize(stat = sqrt(mean(stat, na.rm=na.rm)), .groups = "keep")
  }
  
  return(agg)
}

var_to_fact <- function(df){
  vars <- c("K", "phi", "m", "seed", "H", "r", "nsamp", "sampstrat", "it")
  df[, vars] <- data.frame(lapply(df[, vars], as.factor))
  return(df)
}

format_mmrr <- function(path, full = FALSE){
  df <- read.csv(path)
  
  # remove rows for full
  if(!full) df <- df[df$sampstrat != "full", ]
  
  df <- var_to_fact(df)
  
  df$comboenv_err <- (df$env1_err + df$env2_err)/2
  df$comboenv_coeff <- (df$env1_coeff + df$env2_coeff)/2
  
  df$env_ae <- abs(df$comboenv_err)
  df$geo_ae <- abs(df$geo_err)
  df$ratio_ae <- abs(df$ratio_err)
  
  df$tpIBE <- ((df$env1_p < 0.05) + (df$env2_p < 0.05))/2
  df$tpIBD <- as.numeric(df$geo_p < 0.05)
  
  return(df)
}

format_gdm <- function(path, full = FALSE){
  df <- read.csv(path)
  
  # remove rows for full
  if(!full) df <- df[df$sampstrat != "full", ]
  
  df <- var_to_fact(df)
  
  df$null <- df$geo_coeff == "NULL" & df$env1_coeff == "NULL" & df$env2_coeff == "NULL"
  df[df$null, c("env1_err", 
                "env2_err", 
                "geo_err", 
                "env1_coeff", 
                "env2_coeff", 
                "geo_coeff", 
                "ratio_err", 
                "ratio")] <- NA
  
  df <- var_to_fact(df)
  
  df[, c("env1_err", 
         "env2_err", 
         "geo_err", 
         "env1_coeff", 
         "env2_coeff", 
         "geo_coeff", 
         "ratio_err", 
         "ratio")] <- apply(df[, c("env1_err", 
                             "env2_err", 
                             "geo_err", 
                             "env1_coeff", 
                             "env2_coeff", 
                             "geo_coeff", 
                             "ratio_err", 
                             "ratio")], 2, as.numeric)
  
  df$comboenv_err <- (df$env1_err + df$env2_err)/2
  df$comboenv_coeff <- (df$env1_coeff + df$env2_coeff)/2
  
  df$env_ae <- abs(df$comboenv_err)
  df$geo_ae <- abs(df$geo_err)
  df$ratio_ae <- abs(df$ratio_err)
  
  return(df)
}

format_lfmm <- function(path, full = FALSE){
  
  df <- read.csv(path)
  
  # remove rows for full
  if(!full) df <- df[df$sampstrat != "full", ]
  
  # convert to factors
  df <- var_to_fact(df)
  
  #IMPORTANT NOTE ABOUT LFMM DATA:
  #CURRENTLY IF THERE WERE NOT LOCI IDENTIFIED AS SIGNIFICANT TPR=0 AND FDR = NA (BECAUSE DIVIDE BY 0 = INF = NA). 
  #NOT SURE HOW TO DEAL WITH THIS (I.E. LEAVE AS NA AND EXCLUDE FROM CALCULATIONS OR CONVERT NAS TO 0)
  #I THINK CONVERTING FROM NA TO 0 IS THE MOST ACCURATE (SINCE THE FDR IS 0)
  df$TPRCOMBO[is.na(df$TPRCOMBO)] <- 0
  df$FDRCOMBO[df$TPRCOMBO == 0 & is.na(df$FDRCOMBO)] <- 0
  df$FPRCOMBO[df$FPRCOMBO == 0 & is.na(df$FPRCOMBO)] <- 0
  
  df$roc <- (df$roc1 + df$roc2)/2
  df$pr <- (df$pr1 + df$pr2)/2
  # CHECK THIS: remove rows with all NA
  all_na <- apply(df, 1, function(x) all(is.na(x)))
  df <- df[!all_na,]
  return(df)
}

lfmm_plotter <- function(x, ...) {
  capture.output(grob <- summary_hplot(x, ...))
  method <- unique(x$method)
  K_selection <- unique(x$K_selection)
  title <- paste0("K selection: ", K_selection, " | ", "LFMM method: ", method)
  grob2 <- grid.arrange(grob, top = grid::textGrob(title, x = 0, hjust = 0, gp=gpar(fontsize=20)))
  return(grob2)
}

rda_plotter <- function(x, ...) {
  capture.output(grob <- summary_hplot(x, ...))
  correctPC <- as.logical(unique(x$correctPC))
  if(correctPC) title <- "pRDA" else title <- "RDA"
  grob2 <- grid.arrange(grob, top = grid::textGrob(title, x = 0, hjust = 0, gp=gpar(fontsize=20)))
  return(grob2)
}

format_rda <- function(path, padj = "fdr", alpha = 0.05, full = FALSE){
  
  df <- read.csv(path)
  
  # remove rows for full
  if(!full) df <- df[df$sampstrat != "full", ]
  
  # convert to factors
  df <- var_to_fact(df)
  
  #IMPORTANT NOTE ABOUT RDA DATA:
  #CURRENTLY IF THERE WERE NOT LOCI IDENTIFIED AS SIGNIFICANT TPR=0 AND FDR = NA (BECAUSE DIVIDE BY 0 = INF = NA). 
  #NOT SURE HOW TO DEAL WITH THIS (I.E. LEAVE AS NA AND EXCLUDE FROM CALCULATIONS OR CONVERT NAS TO 0)
  #I THINK CONVERTING FROM NA TO 0 IS THE MOST ACCURATE (SINCE THE FDR IS 0)
  df$TPR[is.na(df$TPR)] <- 0
  df$FDR[df$TPR == 0 & is.na(df$FDR)] <- 0
  df$TPRCOMBO[is.na(df$TPRCOMBO)] <- 0
  df$FDRCOMBO[df$TPRCOMBO == 0 & is.na(df$FDRCOMBO)] <- 0
  df$FPRCOMBO[df$FPRCOMBO == 0 & is.na(df$FPRCOMBO)] <- 0
  
  
  # subset by padj method and alpha
  df <- df[df$padj == padj & alpha == alpha,]
  
  return(df)
}

run_lmer <- function(df, stat, filepath = NULL, seed = 22){
  
  if(!is.null(seed)) set.seed(22)
  
  # mixed effect model
  moddf <- df[df$sampstrat != "full",]
  moddf$stat <- moddf[, stat]
  
  # make nsamp into continuous variable
  moddf$nsamp <- as.numeric(as.character(moddf$nsamp))
  
  # mixed model
  mod <- lmerTest::lmer(stat ~ nsamp + sampstrat + K + m + phi + H + r + (1 | seed), 
                          moddf, na.action = "na.omit", subset = NULL, weights = NULL, offset = NULL)
  
  
  
  
  # print anova result
  print(pretty_anova(mod, filepath))
  
  # print tukey test
  print(pretty_tukey(mod, filepath))
  
}

pretty_tukey <- function(mod, filepath = NULL){
  em <- emmeans::emmeans(mod, pairwise ~ sampstrat, adjust = "tukey")
  em_df <- data.frame(em$contrasts)
  
  d <- max(abs(c(min(em_df$estimate), max(em_df$estimate))))
  
  em_tb <- em_df %>%
    dplyr::select(-df) %>%
    gt::gt() %>%
    cols_label(
      contrast = "Contrast",
      estimate = "Estimate",
      SE = "SE",
      z.ratio = "Z ratio",
      p.value = "p"
    ) %>%
    fmt_number(
      columns = c(estimate, SE, z.ratio),
      decimals = 4,
      suffixing = TRUE
    ) %>%
    tab_style(
      style = list(
        cell_text(weight = "bold")
      ),
      locations = cells_body(
        rows = p.value < 0.05,
        columns = c(p.value, contrast)
      )
    ) %>%
    tab_footnote(
      footnote = "p < 0.05",
      placement = "right",
      locations = cells_body(
        columns = c(p.value, contrast),
        rows = p.value < 0.05 & p.value > 0.01
      )
    ) %>%
    tab_footnote(
      footnote = "p < 0.01",
      placement = "right",
      locations = cells_body(
        columns = c(p.value, contrast),
        rows = p.value < 0.01 & p.value > 0.001
      )
    ) %>%
    tab_footnote(
      footnote = "p < 0.001",
      placement = "right",
      locations = cells_body(
        columns = c(p.value, contrast),
        rows = p.value < 0.001 
      )
    ) %>% fmt_scientific(
      columns = p.value,
      rows = p.value <= 0.009,
      decimals = 1
    ) %>%
    opt_footnote_marks(
      marks = c("***", "**", "*")
      ) %>%
    gtExtras::gt_hulk_col_numeric(
      estimate, 
      trim = TRUE, 
      na.color = "white",
      domain = c(-d, d)
    ) 
  
  if(!is.null(filepath)) write.csv(em_df, gsub(".csv", "_tukey.csv", filepath), row.names = FALSE)
  return(em_tb)
}

pretty_anova <- function(mod, filepath = NULL){
  
  aov <- anova(mod)
  
  effects <- fixef(mod)
  effects_sampstrat <- sum(abs(effects[-which(names(effects) %in% c("Intercept", "nsamp", "K2", "m1", "phi0.5", "H0.5", "r0.6"))]))
  effects <- c(effects["nsamp"], sampstrat = effects_sampstrat, effects[c("K2", "m1", "phi0.5", "H0.5", "r0.6")])
  aov_df <- data.frame(Variable = rownames(aov), FixedEffects = effects, aov)
  
  aov_df$Pr..F. <- signif(aov_df$Pr..F., 2)
  
  d <- max(abs(c(min(aov_df$FixedEffects), max(aov_df$FixedEffects))))
  
  aov_tb <- aov_df %>%
    gt::gt() %>%
    cols_label(
      FixedEffects = "Fixed Effects",
      Variable = "Predictors",
      Sum.Sq = "Sum Sq",
      Mean.Sq = "Mean Sq",
      NumDF = "NumDF",
      DenDF = "DenDF",
      F.value = "F value",
      Pr..F. = "Pr(>F)"
    ) %>%
    fmt_number(
      columns = c(2:4, 6),
      decimals = 4,
      suffixing = TRUE
    ) %>%
    tab_style(
      style = list(
        cell_text(weight = "bold")
      ),
      locations = cells_body(
        columns = c(Pr..F., Variable),
        rows = Pr..F. < 0.05
      )
    ) %>%
    tab_footnote(
      footnote = "p < 0.05",
      placement = "right",
      locations = cells_body(
        columns = Pr..F.,
        rows = Pr..F. < 0.05 & Pr..F. > 0.01
      )
    ) %>%
    tab_footnote(
      footnote = "p < 0.01",
      placement = "right",
      locations = cells_body(
        columns = Pr..F.,
        rows = Pr..F. < 0.01 & Pr..F. > 0.001
      )
    ) %>%
    tab_footnote(
      footnote = "p < 0.001",
      placement = "right",
      locations = cells_body(
        columns = c(Pr..F.),
        rows = Pr..F. < 0.001 
      )
    ) %>% fmt_scientific(
      columns = Pr..F.,
      rows = Pr..F. <= 0.009,
      decimals = 1
    ) %>%
    opt_footnote_marks(
      marks = c("***", "**", "*")
    ) %>%
    gtExtras::gt_hulk_col_numeric(
      FixedEffects, 
      trim = TRUE, 
      na.color = "white",
      domain = c(-d, d)
    ) 
  
  if(!is.null(filepath)) write.csv(aov_df, gsub(".csv", "_lmer.csv", filepath), row.names = FALSE)
  return(aov_tb)
}
