MEGAPLOT <- function(moddf, stat_name, minv = NULL, maxv = NULL, aggfunc = "mean", colpal = "plasma", direction = -1, divergent = FALSE, na.rm=TRUE){
  
  moddf <- moddf[moddf$sampstrat != "full",]
  moddf$stat <- moddf[,stat_name]
  
  
  if(aggfunc == "mean"){
    agg <- moddf %>%
      group_by(K, phi, m, H, r, nsamp, sampstrat) %>%
      summarize(stat = mean(stat, na.rm = na.rm), .groups = "keep")
  }
  
  
  if(aggfunc == "count_na"){
    agg <- moddf %>%
      group_by(K, phi, m, H, r, nsamp, sampstrat) %>%
      summarize(stat = sum(is.na(stat)), .groups = "keep")
  }
  
  if(aggfunc == "count_null"){
    agg <- moddf %>%
      group_by(K, phi, m, H, r, nsamp, sampstrat) %>%
      summarize(stat = sum(is.na(stat)), .groups = "keep")
  }
  
  
  
  if(aggfunc == "prop_na"){
    agg <- moddf %>%
      group_by(K, phi, m, H, r, nsamp, sampstrat) %>%
      summarize(stat = mean(is.na(stat)), .groups = "keep")
  }
  
  if(aggfunc == "var") {
    agg <- moddf %>%
      group_by(K, phi, m, H, r, nsamp, sampstrat) %>%
      summarize(stat = var(stat, na.rm=na.rm), .groups = "keep")
  }
  
  
  
  colnames(agg) <- c("K", "phi", "m", "H", "r", "nsamp", "sampstrat", "mean")
  
  
  #SEED & IT NOT INCLUDED
  params <- expand.grid(K = c(1, 2), 
                        phi = c(0.1, 0.5),
                        m = c(0.25, 1.0),
                        H = c(0.05 , 0.5),
                        r = c(0.3, 0.6))
  
  
  
  # define max and min for plotting
  if(is.null(maxv)){ maxv <- max(agg$mean, na.rm = TRUE)}
  if(is.null(minv)){ minv <- min(agg$mean, na.rm = TRUE)}
  
  plts <- list()
  for(i in 1:nrow(params)){
    tempdf <- merge(params[i,], agg)
    
    ptitle <- paramset <- paste0("K=",params[i,"K"],
                                 " phi=",params[i,"phi"],
                                 " m=",params[i,"m"],
                                 "\nH=",params[i,"H"],
                                 " r=",params[i,"r"])
    
    p <- ggplot(tempdf, aes(nsamp, sampstrat)) +
      ggtitle(ptitle) +
      geom_tile(aes(fill = mean)) + 
      geom_text(aes(label = round(mean, digits = 2), hjust = 0.5), size = 5) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), legend.position = "none",
            axis.title.x=element_blank(), axis.ticks.x=element_blank(),
            axis.title.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.x = element_text(color = "grey50", size = 18),
            axis.text.y = element_text(color = "gray50", size = 18),
            plot.title = element_text(size=20),
            plot.margin=unit(rep(0.4,4),"cm")) +
      coord_fixed()
    
    if(divergent){
      p <- p + scale_fill_gradient2(low = "#2066AC", mid="#F7F7F7", high = "#d79232", midpoint = 0, limits=c(minv, maxv))
    } else {
      p <- p + scale_fill_viridis(limits=c(minv, maxv), option = colpal) 
    }
    
    plts[[i]] <- p
  }
  
  
  bp <- do.call(grid.arrange, c(plts, nrow=4))

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
    ## add justified jitter from the {gghalves} package
    #gghalves::geom_half_point(
    # aes(col = get(var)),
    # ## draw jitter on the left
    # side = "l", 
    # ## control range of jitter
    # range_scale = 0, 
    # ## add some transparency
    # alpha = .1,
    # ## size of points
    # position = position_dodge(width = 0.9)
    #)+
    # color
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

#rewrite in tidy
mean_na <- function(x){
  y <- mean(x, na.rm = TRUE) 
  return(y)
}

sd_na <- function(x){
  y <- sd(x, na.rm = TRUE) 
  return(y)
}

prop_na <- function(x){
  y <- mean(is.na(x))
}

summary_hplot <- function(df, stat_name = "stat", na.rm = TRUE, colpal = "plasma", full=FALSE, sigdig=2, aggfunc = "mean", minv = NULL, maxv = NULL, direction = 1, divergent = FALSE){
  #create column with stat called "stat"
  df$stat <- df[,stat_name]
  
  #remove full data from data frame
  sampstrat <- unique(df$sampstrat)
  nsamp <- unique(df$nsamp)
  if(!full){sampstratsub <- sampstrat[-which(sampstrat=="full")]}
  if(!full){nsampsub <- nsamp[-which(nsamp == 1000 | nsamp == 2000)]}
  
  #create dataframe to store results from summaries
  resdf <- data.frame()
  
  #loop through each number of samples and sampling strategy and summarize the stat
  for(n in nsampsub){
    for(s in sampstratsub){
      #subset the sampling number and strategy
      subdf <- df[df$sampstrat == s & df$nsamp == n, ]
      #create empty dataframe to store results
      meandf <- data.frame()
      #loop through each parameter
      for(p in c("K", "m", "phi", "H", "r")){
        #calculate summary stat for each level of parameter given a function 
        if(aggfunc == "mean"){aggdf <- aggregate(subdf$stat, list(subdf[,p]), mean_na)}
        if(aggfunc == "sd"){aggdf <- aggregate(subdf$stat, list(subdf[,p]), sd_na)}
        if(aggfunc == "prop_na"){aggdf <- aggregate(subdf$stat, list(subdf[,p]), prop_na)}
        if(aggfunc == "max"){aggdf <- aggregate(subdf$stat, list(subdf[,p]), max)}
        if(aggfunc == "min"){aggdf <- aggregate(subdf$stat, list(subdf[,p]), min)}
        if(aggfunc == "var"){aggdf <- aggregate(subdf$stat, list(subdf[,p]), var)}
        
        #create new column for each parameter named p (will be overwritten, probably a cleaner way to do this)
        aggdf[,1] <- paste(p,"=", aggdf[,1])
        aggdf$param <- p
        meandf <- rbind.data.frame(meandf, aggdf)
      }
      
      colnames(meandf) <- c("group", "mean", "param")
      
      #save results
      tempdf <- data.frame(meandf)
      tempdf$nsamp <- n
      tempdf$sampstrat <- s
      resdf <- rbind.data.frame(resdf, tempdf)
    }
  }
  
  # remove row names
  row.names(resdf) <- NULL
  
  # convert nsamp for plotting
  resdf$nsamp <- factor(resdf$nsamp, levels = nsampsub[order(nsampsub)])
  
  # define max and min for plotting
  if(is.null(maxv)){ maxv <- max(resdf$mean, na.rm = TRUE)}
  if(is.null(minv)){ minv <- min(resdf$mean, na.rm = TRUE)}
  
  ## plot data
  plts <- list()
  # For each param value make a plot
  for(i in 1:length(unique(resdf$group))){
    tempdf <- resdf[resdf$group == unique(resdf$group)[i],]
    
    p <- ggplot(tempdf, aes(nsamp, sampstrat)) +
      ggtitle(unique(tempdf$group)) +
      geom_tile(aes(fill = mean)) + 
      geom_text(aes(label = signif(mean, digits = sigdig), hjust = 0.5)) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), legend.position = "none",
            axis.title.x=element_blank(), axis.ticks.x=element_blank(),
            axis.title.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.x = element_text(color = "grey50", size = 14),
            axis.text.y = element_text(color = "gray50", size = 14), 
            plot.margin=unit(rep(0.4,4),"cm")) +
      coord_fixed()
    
    if(divergent){
      p <- p + scale_fill_gradient2(low = "#2066AC", mid="#F7F7F7", high = "#d79232", midpoint = 0, limits=c(minv, maxv))
    } else {
      p <- p + scale_fill_viridis(limits=c(minv, maxv), option = colpal, direction = direction)
    }
    
    plts[[i]] <- p
  }
  
  # This part just plots the pairs of parameters together
  plts2 <- list()
  o <- 1
  for(i in seq(1, length(plts), 2)){
    plts2[[o]] <- do.call(arrangeGrob, c(plts[c(i,i+1)], nrow=2))
    o <- o+1
  }
  
  # Make final plot
  plt <- do.call(grid.arrange, c(plts2, nrow = 1))
}


#Function to determine how sampling strategy affects the relationships between the sim vars/sampstrat
sampstrat_mod <- function(df, alpha = 0.05, padj = "fdr", sampstrat, full = FALSE){
  resdf <- data.frame()
  
  if(!full){sampstrat <- sampstrat[-which(sampstrat=="full")]}
  
  for(s in unique(sampstrat)){
    subdf <- df[df$sampstrat == s, ]
    
    #Need to figure out if you need to subset this anova/remove n.s. variables at this stage or later
    #mixed effect model
    fullmod <- lmer(stat ~ K + m + phi + H + r + (1 | seed), subdf)
    
    #anova to get p-values
    aovmod <- anova(fullmod)
    
    #fixed effects (minus intercept)
    efmod <- fixef(fullmod)[-1]
    
    #save results
    tempdf <- data.frame(fixef = efmod, p = aovmod$`Pr(>F)`)
    tempdf$param <- row.names(aovmod)
    tempdf$sampstrat <- s
    resdf <- rbind.data.frame(resdf, tempdf)
  }
  
  
  row.names(resdf) <- NULL
  resdf$padj <- p.adjust(resdf$p, padj)
  
  tempdf <- resdf
  tempdf[tempdf$padj > alpha, "fixef"] <- NA
  
  #change sampstrat and npts to ordered factor
  tempdf$sampstrat <- factor(tempdf$sampstrat, ordered=TRUE, levels = unique(sampstrat))
  
  p <- ggplot(tempdf, aes(param, sampstrat)) +
    geom_tile(aes(fill = fixef)) + 
    geom_text(aes(label = signif(fixef,  digits = 2), hjust = 0.5)) +
    scale_fill_gradient2(low = "#2066AC", mid="#F7F7F7", high = "#B2182B", midpoint = 0, limits=c(min(resdf$fixef),max(resdf$fixef))) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          #legend.position = "none",
          axis.title.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.ticks.y=element_blank(),
          axis.text.x = element_text(color = "grey50", size = 12),
          axis.text.y = element_text(color = "gray50", size = 12), 
          plot.margin=unit(rep(0.4,4),"cm")) +
    coord_fixed()
  
  print(p)
}


#Function to determine how nsamp affects the relationships between the sim vars/nsamp
nsamp_mod <- function(df, alpha = 0.05, padj = "fdr", nsamp, full = FALSE){
  resdf <- data.frame()
  
  if(!full){nsamp <- nsamp[-which(nsamp==2000)]}
  
  for(n in unique(nsamp)){
    subdf <- df[df$nsamp == n, ]
    
    #Need to figure out if you need to subset this anova/remove n.s. variables at this stage or later
    #mixed effect model
    fullmod <- lmer(stat ~ K + m + phi + H + r + (1 | seed), subdf)
    
    #anova to get p-values
    aovmod <- anova(fullmod)
    
    #fixed effects (minus intercept)
    efmod <- fixef(fullmod)[-1]
    
    #save results
    tempdf <- data.frame(fixef = efmod, p = aovmod$`Pr(>F)`)
    tempdf$param <- row.names(aovmod)
    tempdf$nsamp <- n
    resdf <- rbind.data.frame(resdf, tempdf)
  }
  
  
  row.names(resdf) <- NULL
  resdf$padj <- p.adjust(resdf$p, padj)
  
  tempdf <- resdf
  tempdf[tempdf$padj > alpha, "fixef"] <- NA
  
  #change sampstrat and npts to ordered factor
  tempdf$nsamp <- factor(tempdf$nsamp, ordered=TRUE, levels = nsamp)
  
  p <- ggplot(tempdf, aes(param, nsamp)) +
    geom_tile(aes(fill = fixef)) + 
    geom_text(aes(label = signif(fixef,  digits = 2), hjust = 0.5)) +
    scale_fill_gradient2(low = "#2066AC", mid="#F7F7F7", high = "#B2182B", midpoint = 0, limits=c(min(resdf$fixef),max(resdf$fixef))) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          #legend.position = "none",
          axis.title.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.ticks.y=element_blank(),
          axis.text.x = element_text(color = "grey50", size = 12),
          axis.text.y = element_text(color = "gray50", size = 12), 
          plot.margin=unit(rep(0.4,4),"cm")) +
    coord_fixed()
  
  print(p)
}

#Function to determine how sampling strategy affects the relationships between the sim vars/nsamp
sampstrat_nsamp_mod_params <- function(df, alpha = 0.05, padj = "fdr", sampstrat, nsamp, full = FALSE){
  
  if(!full){sampstratsub <- sampstrat[-which(sampstrat=="full")]}
  if(!full){nsampsub <- nsamp[-which(nsamp==2000)]}
  
  resdf <- data.frame()
  
  for(n in nsampsub){
    for(s in sampstratsub){
      subdf <- df[df$sampstrat == s & df$nsamp == n, ]
      
      #Need to figure out if you need to subset this anova/remove n.s. variables at this stage or later
      #mixed effect model
      fullmod <- lmer(stat ~ K + m + phi + H + r + (1 | seed), subdf)
      
      #anova to get p-values
      aovmod <- anova(fullmod)
      
      #fixed effects (minus intercept)
      efmod <- fixef(fullmod)[-1]
      
      #save results
      tempdf <- data.frame(fixef = fixef(fullmod)[-1], p = aovmod$`Pr(>F)`)
      tempdf$param <- row.names(aovmod)
      tempdf$nsamp <- n
      tempdf$sampstrat <- s
      resdf <- rbind.data.frame(resdf, tempdf)
    }
  }
  
  row.names(resdf) <- NULL
  resdf$nsamp <- as.factor(resdf$nsamp)
  resdf$padj <- p.adjust(resdf$p, padj)
  
  ## plot data
  plts <- list()
  for(i in 1:length(unique(resdf$param))){
    tempdf <- resdf[resdf$param == unique(resdf$param)[i],]
    tempdf[tempdf$padj > alpha, "fixef"] <- NA
    tempdf$lab <- NA
    tempdf$lab[tempdf$padj < 0.10] <- "plain"
    tempdf$lab[tempdf$padj <= 0.05] <-"bold"
    
    p <- ggplot(tempdf, aes(nsamp, sampstrat)) +
      ggtitle(unique(tempdf$param)) +
      geom_tile(aes(fill = fixef)) + 
      geom_text(aes(label = signif(fixef, digits = 1), hjust = 0.5, fontface = lab, size = 10)) +
      scale_fill_gradient2(low = "#2066AC", mid="#F7F7F7", high = "#B2182B", midpoint = 0, limits=c(min(resdf$fixef),max(resdf$fixef))) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), legend.position = "none",
            axis.title.x=element_blank(), axis.ticks.x=element_blank(),
            axis.title.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.x = element_text(color = "grey50", size = 12),
            axis.text.y = element_text(color = "gray50", size = 12), 
            plot.margin=unit(rep(0.4,4),"cm")) +
      coord_fixed()
    
    plts[[i]] <- p
  }
  
  
  do.call(grid.arrange, c(plts, nrow=1))
}

#Function to determine which sampling strategy performs best given a number of points
nsamp_sampstrat_mods <- function(df, alpha = 0.05, padj = "fdr", sampstrat, nsamp){
  
  sampstratsub <- sampstrat[-which(sampstrat=="full")]
  nsampsub <- nsamp[-which(nsamp==2000)]
  
  plts <- list()
  for(i in 1:length(nsampsub)){
    n <- nsampsub[i]
    subdf <- df[df$nsamp == n, ]
    
    #Need to figure out if you need to subset this anova/remove n.s. variables at this stage or later
    #mixed effect model
    fullmod <- lmer(stat ~ sampstrat + K + m + phi + H + r + (1 | seed), subdf)
    
    pw <- emmeans(fullmod, list(pairwise ~ sampstrat), adjust = "tukey")
    
    resdf <- data.frame(pw$`pairwise differences of sampstrat`)
    resdf$p1 <- gsub("\\ -.*", "", resdf$X1)
    resdf$p2 <- gsub(".*- ", "", resdf$X1)
    
    resdf$estimate[resdf$p.value > 0.05] <- NA
    
    resdf$p1 <- factor(resdf$p1, ordered=TRUE, levels = as.character(c("trans","rand","grid","envgeo")))
    resdf$p2 <- factor(resdf$p2, ordered=TRUE, levels = as.character(c("trans","rand","grid","envgeo")))
    
    plt <- ggplot(resdf, aes(p1,p2)) +
      ggtitle(n) +
      geom_tile(aes(fill = estimate)) + 
      geom_text(aes(label = signif(estimate,  digits = 2), hjust = 0.5)) +
      scale_fill_gradient2(low = "#2066AC", mid="#F7F7F7", high = "#B2182B", midpoint = 0, limits=c(min(resdf$estimate),max(resdf$estimate))) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            #legend.position = "none",
            axis.title.x=element_blank(), axis.ticks.x=element_blank(),
            axis.title.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.x = element_text(color = "grey50", size = 12),
            axis.text.y = element_text(color = "gray50", size = 12), 
            legend.position = "none",
            plot.margin=unit(rep(0.4,4),"cm")) +
      coord_fixed()
    
    print(plt)
    
    plts[[i]] <- plt
  }
  
  do.call(grid.arrange, c(plts, nrow=3))
  
}

#Function to determine which sampling strategy performs best given a number of points
sampstrat_nsamp_mods <- function(df, alpha = 0.05, padj = "fdr", sampstrat, nsamp){
  
  sampstratsub <- sampstrat[-which(sampstrat=="full")]
  nsampsub <- nsamp[-which(nsamp==2000)]
  
  plts <- list()
  for(i in 1:length(sampstratsub)){
    s <- sampstratsub[i]
    subdf <- df[df$sampstrat == s, ]
    
    #Need to figure out if you need to subset this anova/remove n.s. variables at this stage or later
    #mixed effect model
    fullmod <- lmer(stat ~ nsamp + K + m + phi + H + r + (1 | seed), subdf)
    
    pw <- emmeans(fullmod, list(pairwise ~ nsamp), adjust = "tukey")
    
    resdf <- data.frame(pw$`pairwise differences of nsamp`)
    resdf$p1 <- gsub("\\ -.*", "", resdf$X1)
    resdf$p2 <- gsub(".*- ", "", resdf$X1)
    
    resdf$estimate[resdf$p.value > 0.05] <- NA
    
    resdf$p1 <- factor(resdf$p1, ordered=TRUE, levels = as.character(c(36, 81, 144, 225, 324)))
    resdf$p2 <- factor(resdf$p2, ordered=TRUE, levels = as.character(c(36, 81, 144, 225, 324)))
    
    plt <- ggplot(resdf, aes(p1,p2)) +
      ggtitle(s) +
      geom_tile(aes(fill = estimate)) + 
      geom_text(aes(label = signif(estimate,  digits = 2), hjust = 0.5)) +
      scale_fill_gradient2(low = "#2066AC", mid="#F7F7F7", high = "#B2182B", midpoint = 0, limits=c(min(resdf$estimate),max(resdf$estimate))) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            #legend.position = "none",
            axis.title.x=element_blank(), axis.ticks.x=element_blank(),
            axis.title.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.x = element_text(color = "grey50"),
            axis.text.y = element_text(color = "gray50"), 
            text=element_text(size=18),
            plot.margin=unit(rep(0.4,4),"cm")) +
      coord_fixed()
    
    print(plt)
    
    plts[[i]] <- plt
  }
  
  do.call(grid.arrange, c(plts, nrow=2))
  
}

var_to_fact <- function(df){
  vars <- c("K", "phi", "m", "seed", "H", "r", "nsamp", "sampstrat", "it")
  df[, vars] <- data.frame(lapply(df[, vars], as.factor))
  return(df)
}

format_mmrr <- function(path){
  df <- read.csv(path)
  
  df <- var_to_fact(df)
  
  df$comboenv_err <- (df$env1_err + df$env2_err)/2
  df$comboenv_coeff <- (df$env1_coeff + df$env2_coeff)/2
  
  df$tpIBE <- ((df$env1_p < 0.05) + (df$env2_p < 0.05))/2
  df$tpIBD <- as.numeric(df$geo_p < 0.05)
  
  df$env_ae <- abs(df$comboenv_err)
  df$geo_ae <- abs(df$geo_err)
  df$ratio_ae <- abs(df$ratio)
  
  return(df)
}

run_lmer <- function(df, stat, table_main = ""){
  # mixed effect model
  moddf <- df[df$sampstrat != "full",]
  moddf$stat <- moddf[, stat]
  
  # na.action needs to be set to na.fail for dredge to work
  fullmod <- lmerTest::lmer(stat ~ nsamp + sampstrat + K + m + phi + H + r + (1 | seed), 
                            moddf, na.action = "na.fail", subset = NULL, weights = NULL, offset = NULL)
  
  # print anova result
  aov <- anova(fullmod)
  print(pretty_anova(aov))
  
  # make HTML table
  #tab_model(fullmod, 
            #p.val = "kr", 
            #show.intercept = FALSE, 
            #dv.labels = table_main,
            #show.re.var = FALSE,
            #show.icc = FALSE,
            #show.r2 = FALSE,
            #show.ngroups = FALSE,
            #show.obs = FALSE,
            #show.reflvl = TRUE,
            #prefix.labels = "varname",
            #digits = 3)
  #
}

pretty_anova <- function(aov){
  stopifnot(class(anova(fullmod))[1] == "anova")
  
  aov_df <- data.frame(Variable = rownames(aov), aov)
  
  aov_df$Pr..F. <- signif(aov_df$Pr..F., 2)
  
  aov_tb <- aov_df %>%
    gt::gt() %>%
    cols_label(
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
      decimals = 2,
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
    opt_footnote_marks(marks = c("***", "**", "*"))
  
  return(aov_tb)
}
