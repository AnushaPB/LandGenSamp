MEGAPLOT <- function(df, stat_name, minv = NULL, maxv = NULL, aggfunc = mean, colpal = "plasma", direction = 1, divergent = FALSE, na.rm=TRUE, dig = 3){
  
  agg <- make_ggdf(df, stat_name = stat_name, aggfunc = aggfunc, na.rm = na.rm)
  
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

make_ggdf <- function(df, stat_name, aggfunc, na.rm = TRUE){
  df %>%
    # convert to data.frame
    data.frame() %>%
    # remove full data
    filter(sampstrat != "full") %>%
    # rename stat_name column as stat for simplified use
    dplyr::rename("stat" = all_of(stat_name)) %>%
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
}

summary_hplot <- function(df, stat_name = "stat", na.rm = TRUE, colpal = "plasma", dig=3, aggfunc = mean, minv = NULL, maxv = NULL, direction = 1, divergent = FALSE, title = NULL, full = FALSE){
 
  # Summarize dataframe
  agg <- df %>% 
    # convert to data.frame
    data.frame() %>%
    # deal with full data
    filter(sampstrat != "full") %>%
    # rename stat_name column as stat for simplified use
    dplyr::rename("stat" = all_of(stat_name)) %>%
    # select columns
    dplyr::select(K, phi, m, H, r, nsamp, sampstrat, stat) %>%
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
    mutate(param = case_when(name == 'K' ~ 'population\nsize',
                             name == 'm' ~ 'migration\nrate',
                             name == 'phi' ~ "selection\nstrength",
                             name == "H" ~ "spatial\nautocorrelation",
                             name == "r" ~ "environmental\ncorrelation", 
                             TRUE ~ "NA")) %>%
    # order param for plotting
    mutate(param = factor(param, ordered = TRUE, levels = c("population\nsize", 
                                                   "migration\nrate", 
                                                   "selection\nstrength", 
                                                   "spatial\nautocorrelation", 
                                                   "environmental\ncorrelation"))) %>%
    # round
    mutate(stat = round(stat, dig))
  
  p <- heat_plot(agg, stat_name = NULL, minv = minv, maxv = maxv, title = title, facet = TRUE, dig = dig, 
                 colpal = colpal, direction = direction, divergent = divergent)
  
  return(p)
}

heat_plot <- function(df, stat_name = NULL, minv = NULL, maxv = NULL, title = NULL, facet = FALSE, dig = 2, 
                      colpal = "plasma", direction = 1, divergent = FALSE){
  
  if (!is.null(stat_name)) df$stat <- df[[stat_name]]
  
  # define max and min for plotting
  if (is.null(maxv)) maxv <- max(df$stat, na.rm = TRUE)
  if (is.null(minv)) minv <- min(df$stat, na.rm = TRUE)
  
  # plot results
  p <- ggplot(df, aes(nsamp, sampstrat)) +
    geom_tile(aes(fill = stat)) + 
    geom_text(aes(label = round(stat, dig), hjust = 0.5)) +
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
          strip.text = element_text(color = "#282828"),
          plot.title = element_text(size=20),
          legend.key.size = unit(1, 'cm'), 
          legend.title = element_text(size=18), 
          legend.text = element_text(size=16)) +
    coord_fixed() 
  
  if (facet) p <- p + facet_grid(low_high ~ param)
  
  if (divergent){
    p <- p + scale_fill_gradient2(name = stat_name, low = "#2066AC", mid="#F7F7F7", high = "#d79232", midpoint = 0, limits=c(minv, maxv))
  } else {
    p <- p + scale_fill_viridis(name = stat_name, limits = c(minv, maxv), option = colpal, direction = direction)
  }
  
  if (!is.null(title)) p <- p + ggtitle(title)
  
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

custom_agg <- function(agg, aggfunc = mean, na.rm = TRUE){
  
  if(is.character(aggfunc)){
    if(aggfunc == "count_na") return(agg %>% summarize(stat = sum(is.na(stat)), .groups = "keep"))
    
    if(aggfunc == "count_null") return(agg %>% summarize(stat = sum(is.na(stat)), .groups = "keep"))
    
    if(aggfunc == "prop_na") return(agg %>% summarize(stat = mean(is.na(stat)), .groups = "keep"))
    
    if(aggfunc == "rmse") {
      return(agg %>%
               mutate(stat = stat^2) %>%
               group_by(K, phi, m, H, r, nsamp, sampstrat) %>%
               summarize(stat = sqrt(mean(stat, na.rm=na.rm)), .groups = "keep"))
    }
  }

  return(agg %>% summarize(stat = aggfunc(stat, na.rm = na.rm), .groups = "keep"))
}

var_to_fact <- function(df){
  mutate_at(df, c("K", "phi", "m", "seed", "H", "r", "nsamp", "sampstrat", "it"), as.factor)
}

format_mmrr <- function(path, full = FALSE){
  df <- read.csv(path)
  df <- extra_mmrr_stats(df)
  
  #give sampling strategies simpler names
  df <- df %>%
    
    # remove full rows
    filter(sampstrat != "full") %>%
    
    # rename sampstrats
    mutate(sampstrat = case_when(
      sampstrat == "envgeo" ~ "EG",
      sampstrat == "rand" ~ "R",
      sampstrat == "trans" ~ "T",
      sampstrat == "grid" ~ "G",
      sampstrat == "equi" ~ "EQ",
      TRUE ~ "NA")) %>%
    
    # convert to factors
    var_to_fact() %>%
    
    # make combo variables
    mutate(RAE = ratio_ae) 
  
  # check number of rows
  if ("trans" %in% df$sampstrat) stopifnot(nrow(df) %% (960 * 4 * 4) == 0)
  if ("equi" %in% df$sampstrat) stopifnot(nrow(df) %% (960 * 3 * 3) == 0)
  
  
  return(df)
}

format_gdm <- function(path, full = FALSE){
  df <- read.csv(path)
  
  #give sampling strategies simpler names
  df <- df %>%
    
    # remove full rows
    filter(sampstrat != "full") %>%
    
    # rename sampstrats
    mutate(sampstrat = case_when(
      sampstrat == "envgeo" ~ "EG",
      sampstrat == "rand" ~ "R",
      sampstrat == "trans" ~ "T",
      sampstrat == "grid" ~ "G",
      sampstrat == "equi" ~ "EQ",
      TRUE ~ "NA")) %>%
    
    # convert to factors
    var_to_fact() %>%
    
    # create NULL column (using geo_coeff, but env_coeff cols will also be NULL)
    mutate(null = case_when(geo_coeff == "NULL" ~ 1, TRUE ~ 0)) %>%

    # convert NULLs to NAs
    mutate_at(c("env1_err", "env2_err", "geo_err", 
                "env1_coeff", "env2_coeff", "geo_coeff", 
                "ratio_err", "ratio"), ~ na_if(., "NULL")) %>%
    
    # convert to numeric
    mutate_at(c("env1_err", "env2_err", "geo_err", 
                "env1_coeff", "env2_coeff", "geo_coeff", 
                "ratio_err", "ratio"), as.numeric) %>%

    # make combo variables
    mutate(comboenv_err = (env1_err + env2_err)/2,
           comboenv_coeff = (env1_coeff + env2_coeff)/2,
           env_ae = abs(comboenv_err),
           geo_ae = abs(geo_err),
           ratio_ae = abs(ratio_err),
           RAE = ratio_ae) 
  
  # check number of rows
  if ("trans" %in% df$sampstrat) stopifnot(nrow(df) %% (960 * 4 * 4) == 0)
  if ("equi" %in% df$sampstrat) stopifnot(nrow(df) %% (960 * 3 * 3) == 0)
  
  return(df)
}


mmrr_gdm_plotter <- function(x, stat, ...) {
  capture.output(grob <- summary_hplot(x, stat_name = stat, ...))
  
  if(stat == "ratio_ae") stat <- "Ratio Absolute Error"
  if(stat == "env_ae") stat <- "IBE Absolute Error"
  if(stat == "geo_ae") stat <- "IBD Absolute Error"
  
  if(stat == "ratio_err") stat <- "Ratio Bias"
  if(stat == "comboenv_err") stat <- "IBE Bias"
  if(stat == "geo_err") stat <- "IBD Bias"
  
  title <- paste0("Statistic: ", stat)
  grob2 <- grid.arrange(grob, top = grid::textGrob(title, x = 0, hjust = 0, gp=gpar(fontsize=20)))
  return(grob2)
}

format_lfmm <- function(path, p_filter = TRUE){
  
  df <- read.csv(path)
  
  #give sampling strategies simpler names
  df <- df %>%
    
    # remove full rows
    filter(sampstrat != "full") %>%

    # rename sampstrats
    mutate(sampstrat = case_when(
      sampstrat == "envgeo" ~ "EG",
      sampstrat == "rand" ~ "R",
      sampstrat == "trans" ~ "T",
      sampstrat == "grid" ~ "G",
      sampstrat == "equi" ~ "EQ",
      TRUE ~ "NA")) %>%
    
    # convert to factors
    var_to_fact() %>%
    
    #IMPORTANT NOTE ABOUT LFMM DATA:
    #CURRENTLY IF THERE WERE NOT LOCI IDENTIFIED AS SIGNIFICANT TPR=0 AND FDR = NA (BECAUSE DIVIDE BY 0 = INF = NA). 
    #NOT SURE HOW TO DEAL WITH THIS (I.E. LEAVE AS NA AND EXCLUDE FROM CALCULATIONS OR CONVERT NAS TO 0)
    #I THINK CONVERTING FROM NA TO 0 IS THE MOST ACCURATE (SINCE THE FDR IS 0)
    mutate(TPRCOMBO = case_when(is.na(TPRCOMBO) ~ 0, TRUE ~ TPRCOMBO)) %>%
    mutate(FDRCOMBO = case_when(TPRCOMBO == 0 & is.na(FDRCOMBO) ~ 0, TRUE ~ FDRCOMBO)) %>%
    mutate(FPRCOMBO = case_when(TPRCOMBO == 0 & is.na(FPRCOMBO) ~ 0, TRUE ~ FPRCOMBO)) %>%
    
    # make ROC and PR columns
    mutate(roc = (roc1 + roc2)/2, pr = (pr1 + pr2)/2) %>%
    
    # make new columns for TPR and FDR with nice names
    mutate(TPR = TPRCOMBO, FDR = FDRCOMBO) 
  
    # filter tested parameters
    
    #filter(padj == "fdr" | is.na(padj)) %>% 
    #filter(sig == 0.05 | is.na(sig)) 
    #filter(K_selection == "tess" | K_selection == "full")
  
  if (p_filter) {
    df <- 
    df %>% 
    # Keep NA to retain NULL models
    filter(padj == "fdr" | is.na(padj)) %>%  
    filter(sig == 0.05 | is.na(sig)) 
  }
  
  # check number of rows
  if ("T" %in% df$sampstrat) stopifnot((nrow(df) %% (960 * 4 * 4)) == 0)
  if ("EQ" %in% df$sampstrat) stopifnot((nrow(df) %% (960 * 3 * 3)) == 0)
  
  return(df)
}

lfmm_plotter <- function(x, stat, ...) {
  capture.output(grob <- summary_hplot(x, stat_name = stat, ...))
  method <- unique(x$method)
  K_selection <- unique(x$K_selection)
  
  if(stat == "K.1") stat <- "K"
  if(stat == "TOTALN") stat <- "Total number of loci"
  if(stat == "TPRCOMBO") stat <- "TPR"
  if(stat == "FDRCOMBO") stat <- "FDR"
  
  if(K_selection == "full") K_selection <- "tracy.widom (full)"
  subtitle <- paste0("K selection: ", K_selection, " | LFMM method: ", method)
  title <- paste0("Statistic: ", stat)
  grob2 <- grid.arrange(grob, top = grid::textGrob(subtitle, x = 0, hjust = 0, gp=gpar(fontsize=20)))
  grob3 <- grid.arrange(grob2, top = grid::textGrob(title, x = 0, hjust = 0, gp=gpar(fontsize=24)))
  return(grob3)
}

rda_plotter <- function(x, stat, ...) {
  capture.output(grob <- summary_hplot(x, stat_name = stat, ...))
  correctPC <- as.logical(unique(x$correctPC))
  if(correctPC) title <- "pRDA" else title <- "RDA"

  if(stat == "TOTALN") stat <- "Total number of loci"
  if(stat == "TPRCOMBO") stat <- "TPR"
  if(stat == "FDRCOMBO") stat <- "FDR"
  
  title <- paste0("Statistic: ", stat, " | ", title)
  grob2 <- grid.arrange(grob, top = grid::textGrob(title, x = 0, hjust = 0, gp=gpar(fontsize=20)))
  return(grob2)
}


format_rda <- function(path, p_filter = TRUE, full = FALSE){
  
  df <- read.csv(path)
  
  #give sampling strategies simpler names
  df <- df %>%
    
    # remove full rows
    filter(sampstrat != "full") %>%
    
    # rename sampstrats
    mutate(sampstrat = case_when(
      sampstrat == "envgeo" ~ "EG",
      sampstrat == "rand" ~ "R",
      sampstrat == "trans" ~ "T",
      sampstrat == "grid" ~ "G",
      sampstrat == "equi" ~ "EQ",
      TRUE ~ "NA")) %>%
    
    # convert to factors
    var_to_fact() %>%
    
    #IMPORTANT NOTE ABOUT RDA DATA:
    #CURRENTLY IF THERE WERE NOT LOCI IDENTIFIED AS SIGNIFICANT TPR=0 AND FDR = NA (BECAUSE DIVIDE BY 0 = INF = NA). 
    #NOT SURE HOW TO DEAL WITH THIS (I.E. LEAVE AS NA AND EXCLUDE FROM CALCULATIONS OR CONVERT NAS TO 0)
    #I THINK CONVERTING FROM NA TO 0 IS THE MOST ACCURATE (SINCE THE FDR IS 0)
    mutate(TPRCOMBO = case_when(is.na(TPRCOMBO) ~ 0, TRUE ~ TPRCOMBO)) %>%
    mutate(FDRCOMBO = case_when(TPRCOMBO == 0 & is.na(FDRCOMBO) ~ 0, TRUE ~ FDRCOMBO)) %>%
    mutate(FPRCOMBO = case_when(TPRCOMBO == 0 & is.na(FPRCOMBO) ~ 0, TRUE ~ FPRCOMBO)) %>%

    # make new columns for TPR and FDR with nice names
    mutate(TPR = TPRCOMBO, FDR = FDRCOMBO) %>%
    
    # remove all rows with all NA values
    filter(if_all(everything(), ~ !is.na(.))) 
  
  
  if (p_filter) {
    df <- 
      df %>% 
      filter(padj == "fdr" | is.na(padj)) %>%  
      filter(sig == 0.05 | is.na(sig)) 
  }
    
  
  # check number of rows
  if ("trans" %in% df$sampstrat) stopifnot(nrow(df) %% (960 * 4 * 4) == 0)
  if ("equi" %in% df$sampstrat) stopifnot(nrow(df) %% (960 * 3 * 3) == 0)
  
  return(df)
}

run_lmer <- function(df, stat, filepath = NULL, seed = 22, f = NULL){
  
  # check number of rows
  if ("trans" %in% df$sampstrat) stopifnot(nrow(df) == (960 * 4 * 4))
  if ("equi" %in% df$sampstrat) stopifnot(nrow(df) == (960 * 3 * 3))
  
  if(!is.null(seed)) set.seed(22)
  
  # mixed effect model
  moddf <- df[df$sampstrat != "full",]
  moddf$stat <- moddf[, stat]
  
  # make nsamp into continuous variable
  moddf$nsamp <- as.numeric(as.character(moddf$nsamp))
  
  # mixed model
  if (is.null(f)) f <- formula("stat ~ nsamp + sampstrat + K + m + phi + H + r + (1 | seed)")
  mod <- lmerTest::lmer(f, moddf, na.action = "na.omit", subset = NULL, weights = NULL, offset = NULL)
  
  # anova 
  gt1 <- pretty_anova(mod, filepath)
  
  # tukey test
  gt2 <- pretty_tukey(mod, filepath)
  
  print(gt1)
  print(gt2)
}

pretty_tukey <- function(mod, filepath = NULL){
  em <- emmeans::emmeans(mod, pairwise ~ sampstrat, adjust = "tukey")
  em_df <- data.frame(em$contrasts)
  
  d <- max(abs(c(min(em_df$estimate), max(em_df$estimate))))
  
  em_tb <- em_df %>%
    dplyr::select(-df) %>%
    gt::gt() %>%
    fmt_number(
      columns = c(2, 3, 4),
      decimals = 4,
      suffixing = TRUE
    ) %>%
    cols_label(
      contrast = "Contrast",
      estimate = "Estimate",
      SE = "SE",
      p.value = "p"
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
    fmt_scientific(
      columns = p.value,
      rows = p.value <= 0.009,
      decimals = 1
    ) %>%
    data_color(
      columns = estimate,
      fn = scales::col_numeric(
        palette = c("#8787fcff","#fcfcfcff", "#f6b26bff"),
        domain = c(-d,d),
        na.color = "white",
      )
    ) 
  #%>%
    #tab_header(
    #  title = md("Tukey test"),
    #  subtitle = md("*pairwise ~ sampstrat*")
    #)
  
  if (any(colnames(em_df) %in% "z.ratio")) em_tb <- em_tb %>% cols_label(z.ratio = "Z ratio")
  if (any(colnames(em_df) %in% "t.ratio")) em_tb <- em_tb %>% cols_label(t.ratio = "t ratio")
  
  if (!is.null(filepath)) write.csv(em_df, gsub(".csv", "_tukey.csv", filepath), row.names = FALSE)
  return(em_tb)
}

pretty_anova <- function(mod, filepath = NULL){
  
  aov <- anova(mod)
  
  effects <- fixef(mod)
  effects_sampstrat <- sum(abs(effects[-which(names(effects) %in% c("Intercept", "nsamp", "K2", "m1", "phi0.5", "H0.5", "r0.6"))]))
  effects <- c(effects["nsamp"], sampstrat = effects_sampstrat, effects[c("K2", "m1", "phi0.5", "H0.5", "r0.6")])
  effects <- na.omit(effects)
  aov_df <- data.frame(Variable = rownames(aov), FixedEffects = effects, aov)
  
  aov_df$Pr..F. <- signif(aov_df$Pr..F., 2)
  
  d <- max(abs(c(min(aov_df$FixedEffects), max(aov_df$FixedEffects))))
  
  aov_df <- 
    aov_df %>%
    mutate(Variable = case_when(Variable == "nsamp" ~ "Sample number",
                                Variable == "K" ~ "Population size", 
                                Variable == "phi" ~ "Selection strength",
                                Variable == "m" ~ "Migration",
                                Variable == "sampstrat" ~ "Sampling strategy",
                                Variable == "H" ~ "Spatial autocorrelation",
                                Variable == "r" ~ "Environmental correlation"))
  
  aov_tb <- aov_df %>%
    filter(Variable != "Sampling strategy") %>%
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
    ) %>% 
    fmt_scientific(
      columns = Pr..F.,
      rows = Pr..F. <= 0.009,
      decimals = 1
    ) %>%
    opt_footnote_marks(
      marks = c("***", "**", "*")
    ) %>%
    data_color(
      columns = FixedEffects,
      fn = scales::col_numeric(
        palette = c("#8787fcff","#fcfcfcff", "#f6b26bff"),
        domain = c(-d,d),
        na.color = "white",
      )
    ) #%>%
    #tab_header(
    #  title = md("Linear mixed effect model"),
    #  subtitle = md("*statistic ~ nsamp + sampstrat + K + m + phi + H + r + (1 | seed)*")
    #)
  
  if(!is.null(filepath)) write.csv(aov_df, gsub(".csv", "_lmer.csv", filepath), row.names = FALSE)
  return(aov_tb)
}


arrange_figures <- function(df1, df2, stat, title1 = NULL, title2 = NULL, ...){
  p1 <- summary_hplot(df1, stat, title = title1, ...)
  p2 <- summary_hplot(df2, stat, title = title2, ...)
  plot <- ggarrange(p1, p2, common.legend = TRUE, nrow = 2, legend = "right")
  return(plot)
}

format_data <- function(method, sampling, p_filter = FALSE) {
  if (sampling == "individual") sampling <- "indsampling"
  if (sampling == "site") sampling <- "sitesampling"
  
  file_name <- paste0(method, "_", sampling, "_results.csv")
  file_path <- here("p3_methods", "outputs", file_name)
  
  if (method == "lfmm") df <- format_lfmm(file_path, p_filter = p_filter)
  if (method == "rda") df <- format_rda(file_path, p_filter = p_filter)
  if (method == "mmrr") df <- format_mmrr(file_path)
  if (method == "gdm") df <- format_gdm(file_path)
  
  return(df)
}


extra_mmrr_stats <- function(df){
  #Create dataframe with all variable combos
  params <- expand.grid(K = c(1, 2), 
                        phi = c(0.1, 0.5),
                        m = c(0.25, 1.0),
                        seed = c(1, 2, 3),
                        H = c(0.05, 0.5),
                        r = c(0.3, 0.6),
                        it = 0:9)
  
  new_ls <- pmap(params, mmrr_stats, df)
  new_df <- bind_rows(new_ls)
  stopifnot(nrow(new_df) == nrow(df))
  
  return(new_df)
}

mmrr_stats <- function(K, phi, m, seed, H, r, it, df){
  sig = 0.05
  subdf <- df[df$K == K & df$phi == phi & df$m == m & df$seed == seed & df$H == H & df$r == r & df$it == it, ]
  
  full <- subdf[subdf$sampstrat == "full",]
  sub <- subdf[subdf$sampstrat != "full",]
  
  # calculate ratio
  full$ratio <- (abs(full$env1_coeff) + abs(full$env2_coeff))/abs(full$geo_coeff)
  sub$ratio <- (abs(sub$env1_coeff) + abs(sub$env2_coeff))/abs(sub$geo_coeff)
  sub$ratio_err <- sub$ratio - full$ratio
  sub$ratio_ae <- abs(sub$ratio_err)
  
  # calculate ratio only including sig values
  # multiply because if the condition is TRUE = 1, env1_coeff = env1_coeff and if FALSE = 0, env1_coeff = 0
  full$env1_coeff_sig <- full$env1_coeff * as.numeric(full$env1_p < sig)
  full$env2_coeff_sig <- full$env2_coeff * as.numeric(full$env2_p < sig)
  full$geo_coeff_sig <- full$geo_coeff * as.numeric(full$geo_p < sig)
  sub$env1_coeff_sig <- sub$env1_coeff * as.numeric(sub$env1_p < sig)
  sub$env2_coeff_sig <- sub$env2_coeff * as.numeric(sub$env2_p < sig)
  sub$geo_coeff_sig <- sub$geo_coeff * as.numeric(sub$geo_p < sig)
  full$ratio_sig <- (abs(full$env1_coeff_sig) + abs(full$env2_coeff_sig))/abs(full$geo_coeff_sig)
  sub$ratio_sig <- (abs(sub$env1_coeff_sig) + abs(sub$env2_coeff_sig))/abs(sub$geo_coeff_sig)
  sub$ratio_sig_err <- sub$ratio_sig - full$ratio_sig
  sub$ratio_sig_ae <- abs(sub$ratio_sig_err)
  
  # calculate IBE error
  sub$env1_err1 <- sub$env1_coeff - full$env1_coeff
  sub$env2_err1 <- sub$env2_coeff - full$env2_coeff
  sub$env1_err2 <- abs(sub$env1_coeff) - abs(full$env1_coeff)
  sub$env2_err2 <- abs(sub$env2_coeff) - abs(full$env2_coeff)
  
  sub$env1_sig_err1 <- sub$env1_coeff_sig - full$env1_coeff_sig
  sub$env2_sig_err1 <- sub$env2_coeff_sig - full$env2_coeff_sig
  sub$env1_sig_err2 <- abs(sub$env1_coeff_sig) - abs(full$env1_coeff_sig)
  sub$env2_sig_err2 <- abs(sub$env2_coeff_sig) - abs(full$env2_coeff_sig)
  
  sub$comboenv_err1 <- (sub$env1_sig_err1 + sub$env2_sig_err1)/2
  sub$comboenv_err2 <- (sub$env1_sig_err2 + sub$env2_sig_err2)/2
  sub$comboenv_sig_err1 <- (sub$env1_sig_err1 + sub$env2_sig_err1)/2
  sub$comboenv_sig_err2 <- (sub$env1_sig_err2 + sub$env2_sig_err2)/2
  
  sub$comboenv_ae1 <- (abs(sub$env1_sig_err1) + abs(sub$env2_sig_err1))/2
  sub$comboenv_ae2 <- (abs(sub$env1_sig_err2) + abs(sub$env2_sig_err2))/2
  sub$comboenv_sig_ae1 <- (abs(sub$env1_sig_err1) + abs(sub$env2_sig_err1))/2
  sub$comboenv_sig_ae2 <- (abs(sub$env1_sig_err2) + abs(sub$env2_sig_err2))/2
  
  # no negative or insignificant values for IBD, so the abs()/sig doesn't matter
  sub$geo_err <- sub$geo_coeff - full$geo_coeff
  sub$geo_ae <- abs(sub$geo_err)
  
  #prop signif
  sub$prop_comboenv_sig <- ((sub$env1_p < sig) + (sub$env2_p < sig))/2
  sub$prop_geo_sig <- (sub$geo_p < sig)
  full$prop_comboenv_sig <- ((full$env1_p < sig) + (full$env2_p < sig))/2
  full$prop_geo_sig <- (full$geo_p < sig)
  
  sub$prop_comboenv_sigfull <- sub$prop_comboenv_sig == full$prop_comboenv_sig 
  sub$prop_geo_sigfull <- sub$prop_geo_sig == full$prop_geo_sig
  
 # sub$prop_comboenv_sigfull1 <-(((sub$env1_p < sig) == (full$env1_p < sig)) + ((sub$env2_p < sig) == (full$env2_p < sig)))/2
  #sub$prop_comboenv_sigfull2 <-(((sub$env1_p < sig) & (full$env1_p < sig)) + ((sub$env2_p < sig) & (full$env2_p < sig)))/2
  #sub$prop_comboenv_sigfull3 <-(((sub$env1_p < sig) & (full$env1_p < sig)) + ((sub$env2_p < sig) & (full$env2_p < sig)))/(full$env1_p < sig + full$env2_p < sig)
  sub$prop_comboenv_sigfull_env1 <-((sub$env1_p < sig) & (full$env1_p < sig))/(full$env1_p < sig)
  sub$prop_comboenv_sigfull_env2 <-((sub$env2_p < sig) & (full$env2_p < sig))/(full$env2_p < sig)
  sub$prop_comboenv_sigfull_combo <- (sub$prop_comboenv_sigfull_env1 + sub$prop_comboenv_sigfull_env2)/2
  new_df <- bind_rows(sub, full)
  stopifnot(nrow(new_df) == nrow(subdf))
  
  # negative value
  
  return(new_df)
}
