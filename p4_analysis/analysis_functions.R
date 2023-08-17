
MEGAPLOT <- function(df, stat, minv = NULL, maxv = NULL, aggfunc = mean, colpal = NULL, direction = 1, na.rm=TRUE, dig = 3, pretty_names = TRUE){
  stat_name <- stat
  agg <- make_ggdf(df, stat_name = stat_name, aggfunc = aggfunc, na.rm = na.rm)
  
  if (pretty_names) stat_name <- make_pretty_names(stat_name)
  
  # define max and min for plotting
  if (is.null(maxv)) maxv <- max(agg$stat, na.rm = TRUE)
  if (is.null(minv)) minv <- min(agg$stat, na.rm = TRUE)
  
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
          strip.text = element_text(color = "black"),
          plot.title = element_text(size = 30, face = "bold"))
  
  # get default color palettes based on stats
  if (is.null(colpal)){
    pal <- get_colpal(stat_name)
    colpal <- pal$colpal
    direction <- pal$direction
  }
  
  if(colpal == "divergent"){
    p <- p + scale_fill_gradient2(low = "#2066AC", mid="#F7F7F7", high = "#d79232", midpoint = 0, limits=c(minv, maxv))
  } else {
    p <- p + scale_fill_viridis(limits = c(minv, maxv), option = colpal, direction = direction)
  }

  p <- p + ggtitle(stat_name)
  
  return(p)
  
}

make_pretty_names <- function(stat_name) {
  new_name <- ""
  if (stat_name == "K_factor") new_name <- "Number of latent factors"
  if (stat_name == "TOTALN") new_name <- "Total number of loci"
  if (stat_name == "TPRCOMBO") new_name <- "TPR"
  if (stat_name == "FDRCOMBO") new_name <- "FDR"
  if (stat_name == "FDRCOMBO") new_name <- "FDR"
  if (grepl("*_coeff_err", stat_name) | grepl("*coefferr", stat_name)) new_name <- "Bias"
  if (grepl("*_coeff_err_ae", stat_name) | grepl("*coefferrae", stat_name)) new_name <- "AE"
  if (grepl("*TPR", stat_name)) new_name <- "TPR"
  if (grepl("*FDR", stat_name)) new_name <- "FDR"
  if (grepl("*geo*", stat_name)) new_name <- paste("IBD", new_name)
  if (grepl("*env*", stat_name)) new_name <- paste("IBE", new_name)
  if (grepl("relaxed", stat_name)) new_name <- paste(new_name, "relaxed")
  if (grepl("strict", stat_name)) new_name <- paste(new_name, "strict")
  if (new_name == "") new_name <- stat_name
  return(new_name)
}

get_colpal <- function(stat_name){
  if (is.null(stat_name)) return(list(colpal = "cividis", direction = 1))
  if (stat_name == "K_factor" | stat_name == "TOTALN") return(list(colpal = "cividis", direction = 1))
  if (grepl("*_coeff_err_ae", stat_name) | grepl("*FDR*", stat_name)) return(list(colpal = "viridis", direction = -1))
  if (grepl("*_coeff_err", stat_name)) return(list(colpal = "divergent"))
  if (grepl("*TPR*", stat_name)) return(list(colpal = "plasma", direction = 1))
  return(list(colpal = "cividis", direction = 1))
}

make_lmer_path <- function(method, sampling, stat){
  stat <- gsub("_", "", stat)
  filepath <- paste0(method, "_", sampling, "_", stat, ".csv")
  return(here(p4path, filepath))
}

make_ggdf <- function(df, stat_name, aggfunc = mean, na.rm = TRUE){
  ggdf <- 
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
    # mutate to numeric
    mutate_at(c("K", "phi", "m", "H", "r"), ~as.numeric(as.character(.))) %>%
    # mutate params to get low and high (not sure why but mutate_at wasn't working)
    mutate(
      K = case_when(K == min(.[["K"]]) ~ "L", K == max(.[["K"]]) ~ "H"),
      m = case_when(m == min(.[["m"]]) ~ "L", m == max(.[["m"]]) ~ "H"),
      phi = case_when(phi == min(.[["phi"]]) ~ "L", phi == max(.[["phi"]]) ~ "H"),
      H = case_when(H == min(.[["H"]]) ~ "L", H == max(.[["H"]]) ~ "H"),
      r = case_when(r == min(.[["r"]]) ~ "L", r == max(.[["r"]]) ~ "H"),
    ) %>% 
  # create new group for plotting
    mutate(group = paste0("K=", K,
                          " phi=", phi,
                          " m=", m,
                          "\nH=", H,
                          " r=", r))
  
  return(ggdf)
}

summary_hplot <- function(df, stat_name = "stat", na.rm = TRUE, colpal = NULL, dig=3, aggfunc = mean, minv = NULL, maxv = NULL, direction = 1, title = NULL, full = FALSE, pretty_names = TRUE){
 
  agg <- make_agg(df, stat_name = stat_name, na.rm = na.rm, dig = dig, aggfunc = aggfunc, minv = minv, maxv = maxv, full = full)
  
  # get default color palettes based on stats
  if (is.null(colpal)){
    pal <- get_colpal(stat_name)
    colpal <- pal$colpal
    direction <- pal$direction
  }
  
  p <- heat_plot(agg, stat_name = NULL, minv = minv, maxv = maxv, title = title, facet = TRUE, dig = dig, 
                 colpal = colpal, direction = direction)
  
  return(p)
}

make_agg <- function(df, stat_name = "stat", na.rm = TRUE, dig = 3, aggfunc = mean, minv = NULL, maxv = NULL, full = FALSE){
  
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
  
  return(agg)
}

agg_mm <- function(x, stat){
  result <-
    map(x, \(df, y = stat) {
      agg <- make_agg(df, stat_name = stat)
      list(min = min(agg[["stat"]], na.rm = TRUE), max = max(agg[["stat"]], na.rm = TRUE))
    }) %>%
    transpose() %>%
    map(unlist)
  
  return(c(min = min(result$min, na.rm = TRUE), max = max(result$max, na.rm = TRUE)))
}


heat_plot <- function(df, stat_name = NULL, minv = NULL, maxv = NULL, title = NULL, facet = FALSE, dig = 2, 
                      colpal = NULL, direction = 1){
  
  if (!is.null(stat_name)) df$stat <- df[[stat_name]]
  
  # define max and min for plotting
  if (is.null(maxv)) maxv <- max(df$stat, na.rm = TRUE)
  if (is.null(minv)) minv <- min(df$stat, na.rm = TRUE)
  
  # plot results
  p <- ggplot(df, aes(nsamp, sampstrat)) +
    geom_tile(aes(fill = stat)) + 
    geom_text(aes(label = round(stat, dig), hjust = 0.5), size = 5) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), legend.position = "none",
          axis.title.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.ticks.y=element_blank(),
          axis.text.x = element_text(color = "grey50", size = 16),
          axis.text.y = element_text(color = "gray50", size = 16), 
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
  
  # get default color palettes based on stats
  if (is.null(colpal)){
    pal <- get_colpal(stat_name)
    colpal <- pal$colpal
    direction <- pal$direction
  }
  
  if (colpal == "divergent"){
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
  df <- extra_ibeibd_stats(df)
  
  if (!full) df <- df %>% filter(sampstrat != "full") 
  #give sampling strategies simpler names
  df <- df %>%
    
    # rename sampstrats
    mutate(sampstrat = case_when(
      sampstrat == "envgeo" ~ "EG",
      sampstrat == "rand" ~ "R",
      sampstrat == "trans" ~ "T",
      sampstrat == "grid" ~ "G",
      sampstrat == "equi" ~ "EQ",
      sampstrat == "full" ~ "full",
      TRUE ~ "NA")) %>%
    
    # convert to factors
    var_to_fact() 
  
  # make combos
  if (any(colnames(df) == "env1_coeff")) {
    df <- 
      df %>%  
      # make combo variables
      mutate(comboenv_coeff_err = (env1_coeff_err + env2_coeff_err)/2,
             comboenv_coeff = (env1_coeff + env2_coeff)/2,
             env_coeff_err = comboenv_coeff_err,
             env_coeff_err_ae = abs(comboenv_coeff_err),
             # needs to be row means to remove NAs
             env_p_TPR = rowMeans(.[,c("env1_p_TPR", "env2_p_TPR")], na.rm = TRUE),
             env1_TPR_propNA = is.na(env1_p_TPR),
             env2_TPR_propNA = is.na(env2_p_TPR),
             env_TPR_propNA = (env1_TPR_propNA + env2_TPR_propNA)/2)

  }
  
  # check number of rows
  if ("trans" %in% df$sampstrat) stopifnot(nrow(df) %% (960 * 4 * 4) == 0)
  if ("equi" %in% df$sampstrat) stopifnot(nrow(df) %% (960 * 3 * 3) == 0)
  
  
  return(df)
}

format_gdm <- function(path, full = FALSE){
  df <- read.csv(path)
  
  df <- extra_ibeibd_stats(df)
  
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
    mutate(null_geo_coeff = case_when(is.na(geo_coeff) ~ 1, TRUE ~ 0),
           null_geo_p = case_when(is.na(geo_p) ~ 1, TRUE ~ 0),
           null_env1_p = case_when(is.na(env1_p) ~ 1, TRUE ~ 0),
           null_env2_p = case_when(is.na(env2_p) ~ 1, TRUE ~ 0),
           null_env_p = null_env1_p + null_env2_p
           )

  if (any(colnames(df) == "env1_coeff")) {
    df <- df %>%  
      # make combo variables
      mutate(env_coeff_err = (env1_coeff_err + env2_coeff_err)/2,
             env_coeff = (env1_coeff + env2_coeff)/2,
             env_coeff_err_ae = abs(env_coeff_err),
             env_p_TPR = (env1_p_TPR + env2_p_TPR)/2) 
  }
  
  # create NA to 0
  df <- 
    df %>%
    mutate(env_p_TPR0 = case_when(is.na(env_p_TPR) ~ 0, .default = env_p_TPR),
           geo_p_TPR0 = case_when(is.na(geo_p_TPR) ~ 0, .default = geo_p_TPR))
  
  # check number of rows
  if ("trans" %in% df$sampstrat) stopifnot(nrow(df) %% (960 * 4 * 4) == 0)
  if ("equi" %in% df$sampstrat) stopifnot(nrow(df) %% (960 * 3 * 3) == 0)
  
  return(df)
}


mmrr_gdm_plotter <- function(x, stat, ...) {
  capture.output(grob <- summary_hplot(x, stat_name = stat, ...))
  title <- paste0("Statistic: ", make_pretty_names(stat))
  grob2 <- grid.arrange(grob, top = grid::textGrob(title, x = 0, hjust = 0, gp=gpar(fontsize=20)))
  return()
}

format_lfmm <- function(path, p_filter = TRUE){
  
  df <- read.csv(path)
  
  #give sampling strategies simpler names
  df <-
    df %>%
    
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
    #I THINK CONVERTING FROM NA TO 0 IS THE MOST ACCURATE (SINCE THE FDR IS 0)
    mutate(TPRCOMBO = case_when(is.na(TPRCOMBO) ~ 0, TRUE ~ TPRCOMBO)) %>%
    mutate(FDRCOMBO = case_when(TPRCOMBO == 0 & is.na(FDRCOMBO) ~ 0, TRUE ~ FDRCOMBO)) %>%
    mutate(FPRCOMBO = case_when(TPRCOMBO == 0 & is.na(FPRCOMBO) ~ 0, TRUE ~ FPRCOMBO)) %>%
  
    # filter by MAF
    filter(maf == 0.05) %>%
    
    # convert all to strict/relaxed
    mutate(all = case_when(all == TRUE ~ "relaxed", all == FALSE ~ "strict")) %>%
    
    # select stats for analysis
    select(K, phi, m, seed, H, r, it, sampstrat, nsamp, K_method, lfmm_method, maf, sig, padj, all, TPRCOMBO, FDRCOMBO, TOTALN, K_factor) %>%

    # widen stats
    pivot_wider(names_from = all, 
                values_from = c("TPRCOMBO", "FDRCOMBO", "TOTALN"), 
                names_glue = "{.value}_{all}") %>%
    
    # TEMPORARY: TOTALN for all = FALSE/strict was incorrectly calculated (didn't use unique loci)
    mutate(TOTALN = TOTALN_relaxed) %>%
    select(-TOTALN_relaxed, -TOTALN_strict)
  
  
  if (p_filter) {
    df <- 
      df %>% 
      # Keep NA to retain NULL models
      filter(padj == "fdr" | is.na(padj)) %>%  
      filter(sig == 0.05 | is.na(sig)) 
  }
  
  # filter by method for row count
  df2 <- df %>% filter(lfmm_method == "ridge")
  
  # check number of rows
  if ("T" %in% df$sampstrat) stopifnot((nrow(df2) %% (96 * 4 * 4)) == 0)
  if ("EQ" %in% df$sampstrat) stopifnot((nrow(df2) %% (96 * 3 * 3)) == 0)
  
  return(df)
}

lfmm_plotter <- function(x, stat, ...) {
  capture.output(grob <- summary_hplot(x, stat_name = stat, ...))
  lfmm_method <- unique(x$lfmm_method)
  stat <- make_pretty_names(stat)
  subtitle <- paste0("LFMM method: ", lfmm_method)
  title <- paste0(stat)
  grob2 <- grid.arrange(grob, top = grid::textGrob(subtitle, x = 0, hjust = 0, gp=gpar(fontsize=20)))
  grob3 <- grid.arrange(grob2, top = grid::textGrob(title, x = 0, hjust = 0, gp=gpar(fontsize=24)))
  return(grob3)
}

rda_plotter <- function(x, stat, ...) {
  capture.output(grob <- summary_hplot(x, stat_name = stat, ...))
  correctPC <- as.logical(unique(x$correctPC))
  if(correctPC) rda_method <- "partial (PC correction)" else rda_method <- "standard (no PC correction)"
  stat <- make_pretty_names(stat)
  subtitle <- paste0("RDA method: ", rda_method)
  title <- paste0(stat)
  grob2 <- grid.arrange(grob, top = grid::textGrob(subtitle, x = 0, hjust = 0, gp=gpar(fontsize=20)))
  grob3 <- grid.arrange(grob2, top = grid::textGrob(title, x = 0, hjust = 0, gp=gpar(fontsize=24)))
  return(grob3)
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
    mutate(FPRCOMBO = case_when(TPRCOMBO == 0 & is.na(FPRCOMBO) ~ 0, TRUE ~ FPRCOMBO),
           FDRCOMBO_rdap = case_when(TPRCOMBO_rdap == 0 & is.na(FDRCOMBO_rdap) ~ 0, TRUE ~ FDRCOMBO_rdap)) %>%

    # make new columns for TPR and FDR with nice names
    mutate(TPR = TPRCOMBO, FDR = FDRCOMBO)  %>%
    
    # fix all
    mutate(all = case_when(is.na(all) ~ FALSE, .default = all)) %>%
    
    # filter by MAF 
    filter(maf == 0.05) %>%
    
    # convert all to strict/relaxed
    mutate(all = case_when(all == TRUE ~ "relaxed", all == FALSE ~ "strict")) %>%
    
    # select stats for analysis
    select(K, phi, m, seed, H, r, it, sampstrat, nsamp, correctPC, maf, sig, padj, all, TPRCOMBO, FDRCOMBO, TOTALN) %>%
    
    # widen stats
    pivot_wider(names_from = all, 
                values_from = c("TPRCOMBO", "FDRCOMBO"), 
                names_glue = "{.value}_{all}") 
  
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
  
  if (!is.null(seed)) set.seed(seed)
  
  # mixed effect model
  df <- data.frame(df)
  moddf <- df[df$sampstrat != "full",]
  moddf$stat <- moddf[, stat]
  
  # check if all 1 value
  if (length(unique(moddf$stat)) == 1) {
    aov_df <- data.frame(Variable = NA, FixedEffects = NA, Sum.Sq = NA, Mean.Sq = NA, NumDF = NA, DenDF = NA, F.value = NA, Pr..F. = NA)
    if(!is.null(filepath)) write.csv(aov_df, gsub(".csv", "_lmer.csv", filepath), row.names = FALSE)
    em_df <- data.frame(contrast = NA, estimate = NA, SE = NA, df = NA, z.ratio = NA, p.value = NA)
    if(!is.null(filepath)) write.csv(em_df, gsub(".csv", "_tukey.csv", filepath), row.names = FALSE)
    warning(paste("\nAll values of", stat, "are", unique(moddf$stat), "; no model is returned"))
    return()
  }
  
  # make nsamp into continuous variable
  moddf$nsamp <- as.numeric(as.character(moddf$nsamp))
  
  # mixed model
  if (is.null(f)) f <- formula("stat ~ nsamp + sampstrat + K + m + phi + H + r + (1 | seed)")
  mod <- lmerTest::lmer(f, moddf, na.action = "na.omit", subset = NULL, weights = NULL, offset = NULL)
  
  # anova 
  gt1 <- pretty_anova(mod, filepath, stat = stat)
  
  # tukey test
  gt2 <- pretty_tukey(mod, filepath, stat = stat)
  
  print(gt1)
  print(gt2)
  return()
}

pretty_tukey <- function(mod, filepath = NULL, stat = "stat"){
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
    ) %>%
    tab_header(
      title = md(paste0("Tukey test for ",  make_pretty_names(stat))),
      subtitle = md("*pairwise ~ sampstrat*")
    )
  
  if (any(colnames(em_df) %in% "z.ratio")) em_tb <- em_tb %>% cols_label(z.ratio = "Z ratio")
  if (any(colnames(em_df) %in% "t.ratio")) em_tb <- em_tb %>% cols_label(t.ratio = "t ratio")
  
  if (!is.null(filepath)) write.csv(em_df, gsub(".csv", "_tukey.csv", filepath), row.names = FALSE)
  return(em_tb)
}

pretty_anova <- function(mod, filepath = NULL, stat = "stat"){
  
  aov <- anova(mod)
  
  effects <- fixef(mod)
  effects_sampstrat <- sum(abs(effects[-which(names(effects) %in% c("Intercept", "nsamp", "K2", "m1", "phi1", "H0.5", "r0.6"))]))
  effects <- c(effects["nsamp"], sampstrat = effects_sampstrat, effects[c("K2", "m1", "phi1", "H0.5", "r0.6")])
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
      marks = c("*", "**", "***")
    ) %>%
    data_color(
      columns = FixedEffects,
      fn = scales::col_numeric(
        palette = c("#8787fcff","#fcfcfcff", "#f6b26bff"),
        domain = c(-d,d),
        na.color = "white",
      )
    ) %>%
    tab_header(
      title = md("Linear mixed effect model"),
      subtitle = md(paste0(make_pretty_names(stat), " ~ nsamp + sampstrat + K + m + phi + H + r + (1 | seed)"))
    )
  
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
  path <- here("p3_methods", "outputs", file_name)
  
  if (method == "lfmm") df <- format_lfmm(path, p_filter = p_filter)
  if (method == "rda") df <- format_rda(path, p_filter = p_filter)
  if (method == "mmrr") df <- format_mmrr(path)
  if (method == "mmrr2") df <- format_mmrr(path)
  if (method == "gdm") df <- format_gdm(path)
  if (method == "gdm2") df <- format_gdm(path)
  
  return(df)
}

extra_ibeibd_stats <- function(df){
  #Create dataframe with all variable combos
  params <- expand.grid(K = unique(df$K), 
                        phi = unique(df$phi),
                        m = unique(df$m),
                        seed = unique(df$seed),
                        H = unique(df$H),
                        r = unique(df$r),
                        it = unique(df$it))
  
  new_ls <- pmap(params, ibeibd_stats, df)
  new_df <- bind_rows(new_ls)
  stopifnot(nrow(new_df) == nrow(df))
  
  return(new_df)
}

ibeibd_stats <- function(K, phi, m, seed, H, r, it, df){
  subdf <- df[df$K == K & df$phi == phi & df$m == m & df$seed == seed & df$H == H & df$r == r & df$it == it, ]
  
  full <- subdf[subdf$sampstrat == "full",]
  sub <- subdf[subdf$sampstrat != "full",]
  
  if (any(colnames(df) == "env1_coeff")){
    
    # calculate SMAPE (deal with zeros)
    sub$env1_smape <- smape(full$env1_coeff, sub$env1_coeff)
    sub$env2_smape <- smape(full$env2_coeff, sub$env2_coeff)
    sub$env_smape <- (sub$env1_smape + sub$env2_smape)/2
    
    # calculate FDR
    sub$env1_p_FDR <- (!(full$env1_p < 0.05) & sub$env1_p < 0.05)/(sub$env1_p < 0.05)
    sub$env1_p_FDR[is.na(sub$env1_p_FDR)] <- 0 
    sub$env2_p_FDR <- (!(full$env2_p < 0.05) & sub$env2_p < 0.05)/(sub$env2_p < 0.05)
    sub$env2_p_FDR[is.na(sub$env2_p_FDR)] <- 0 
    sub$env_p_FDR <- (sub$env1_p_FDR + sub$env2_p_FDR)/2
    sub$geo_p_FDR <- (!(full$geo_p < 0.05) & sub$geo_p < 0.05)/(sub$geo_p < 0.05)
    sub$geo_p_FDR[is.na(sub$geo_p_FDR)] <- 0 
    
    # calculate IBE
    full$IBE <- abs(full$env1_coeff) + abs(full$env2_coeff)
    sub$IBE <- abs(sub$env1_coeff) + abs(sub$env2_coeff)
    sub$IBE_err_ae <- abs(full$IBE - sub$IBE)
    sub$IBE_smape <- smape(full$IBE, sub$IBE)
    
  } 
    
  sub$geo_rpd <- 2 * abs(sub$geo_coeff - full$geo_coeff)/(abs(full$geo_coeff) + abs(sub$geo_coeff))
  
  new_df <- bind_rows(sub, full)
  stopifnot(nrow(new_df) == nrow(subdf))
  
  return(new_df)
}

# SMAPE: https://www.r-bloggers.com/2021/08/how-to-calculate-smape-in-r/
# a = actual, f = forecasted
smape <- function(a, f) {  return (1/length(a) * sum(2*abs(f-a) / (abs(a)+abs(f))*100))}
