
MEGAPLOT <- function(moddf, stat, minv = 0, maxv = max(moddf[,stat]), option = "plasma"){
  meanagg <- aggregate(moddf[,stat], list(moddf$K, moddf$phi, moddf$m, moddf$H, moddf$r, moddf$nsamp, moddf$sampstrat), mean)
  colnames(meanagg) <- c("K", "phi", "m", "H", "r", "nsamp", "sampstrat", "mean")
  
  
  params <- expand.grid(K = c(2, 4), 
                        phi = c(0.1, 0.5),
                        m = c(0.25, 1.0))
  
  plts <- list()
  for(i in 1:nrow(params)){
    tempdf <- merge(params[i,], meanagg)
    
    ptitle <- paramset <- paste0("K=",params[i,"K"],
                                 " phi=",params[i,"phi"],
                                 " m=",params[i,"m"],
                                 "\nH=",params[i,"H"],
                                 " r=",params[i,"r"])
    
    p <- ggplot(tempdf, aes(nsamp, sampstrat)) +
      ggtitle(ptitle) +
      geom_tile(aes(fill = mean)) + 
      geom_text(aes(label = round(mean, digits = 2), hjust = 0.5), size = 5) +
      scale_fill_viridis(limits=c(minv, maxv), option = option) +
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
    
    plts[[i]] <- p
  }
  
  
  bp <- do.call(grid.arrange, c(plts, nrow=12))
}


summary_vplot <- function(df, allplots = TRUE, colpal = "plasma"){
  (pK <- ggplot(df, aes(fill=K, y=stat, x=factor(nsamp))) + 
     #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
     scale_fill_viridis(discrete=T, option=colpal) +
     xlab("nsamp") +
     theme_bw()+
     geom_boxplot(width = 0.4, position = position_dodge(width = 0.9))+
     #geom_point(aes(col=K),position = position_dodge(width = 0.9), alpha=0.1)+
     scale_colour_viridis(discrete=T, option=colpal) +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=21)))
  
  (pphi <- ggplot(df, aes(fill=phi, y=stat, x=factor(nsamp))) + 
      #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
      scale_fill_viridis(discrete=T, option=colpal) +
      xlab("nsamp") +
      theme_bw()+
      geom_boxplot(width = 0.4, position = position_dodge(width = 0.9))+
      #geom_point(aes(col=phi),position = position_dodge(width = 0.9), alpha=0.1)+
      scale_colour_viridis(discrete=T, option=colpal) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=21)))
  
  (pm <- ggplot(df, aes(fill=m, y=stat, x=factor(nsamp))) + 
      #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
      scale_fill_viridis(discrete=T, option=colpal) +
      xlab("nsamp") +
      theme_bw()+
      geom_boxplot(width = 0.4, position = position_dodge(width = 0.9))+
      scale_colour_viridis(discrete=T, option=colpal) +
      #geom_point(aes(col=m),position = position_dodge(width = 0.9), alpha=0.1)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=21)))
  
  (ps <- ggplot(df, aes(fill=sampstrat, y=stat, x=factor(nsamp))) + 
      #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
      scale_fill_viridis(discrete=T, option=colpal) +
      xlab("nsamp") +
      theme_bw()+
      geom_boxplot(width = 0.4, position = position_dodge(width = 0.9))+
      #geom_point(aes(fill=sampstrat, col=sampstrat),position = position_dodge(width = 0.9), alpha=0.1)+
      scale_colour_viridis(discrete=T, option=col) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=21)))
  
  if(allplots){
    print(pK)
    print(pphi)
    print(pm)
    print(ps)
  }

  plt <- grid.arrange(pK, pphi, pm, ps, nrow = 2)
  return(plt)
}


summary_hplot <- function(df, colpal = "plasma", full=FALSE){
  
  if(!full){sampstratsub <- sampstrat[-which(sampstrat=="full")]}
  if(!full){nsampsub <- nsamp[-which(nsamp==2000)]}
  
  resdf <- data.frame()
  
  for(n in nsampsub){
    for(s in sampstratsub){
      subdf <- df[df$sampstrat == s & df$nsamp == n, ]
      
      meandf <- data.frame()
      for(p in c("K", "m", "phi")){
        aggdf <- aggregate(subdf$stat, list(subdf[,p]), mean)
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
  
  row.names(resdf) <- NULL
  resdf$nsamp <- as.factor(resdf$nsamp)
  
  ## plot data
  plts <- list()
  for(i in 1:length(unique(resdf$group))){
    tempdf <- resdf[resdf$group == unique(resdf$group)[i],]
    p <- ggplot(tempdf, aes(nsamp, sampstrat)) +
      ggtitle(unique(tempdf$group)) +
      geom_tile(aes(fill = mean)) + 
      geom_text(aes(label = signif(mean, digits = 2), hjust = 0.5)) +
      scale_fill_viridis(limits=c(min(resdf$mean),max(resdf$mean)), option = colpal) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), legend.position = "none",
            axis.title.x=element_blank(), axis.ticks.x=element_blank(),
            axis.title.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.x = element_text(color = "grey50", size = 14),
            axis.text.y = element_text(color = "gray50", size = 14), 
            plot.margin=unit(rep(0.4,4),"cm")) +
      coord_fixed()
    
    plts[[i]] <- p
  }
  
  
  plts2 <- list()
  o <- 1
  for(i in seq(1, length(plts), 2)){
    plts2[[o]] <- do.call(grid.arrange, c(plts[c(i,i+1)], nrow=2))
    o <- o+1
  }
  
  plt <- do.call(grid.arrange, c(plts2, nrow=1))
  return(plt)
}


#Function to determine how sampling strategy affects the relationships between the sim vars/sampstrat
sampstrat_mod <- function(df, alpha = 0.05, padj = "fdr", sampstrat, full = FALSE){
  resdf <- data.frame()
  
  if(!full){sampstrat <- sampstrat[-which(sampstrat=="full")]}
  
  for(s in unique(sampstrat)){
    subdf <- df[df$sampstrat == s, ]
    
    #Need to figure out if you need to subset this anova/remove n.s. variables at this stage or later
    #mixed effect model
    fullmod <- lm(stat ~ K + m + phi, subdf)
    
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
    fullmod <- lm(stat ~ K + m + phi, subdf)
    
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
      fullmod <- lm(stat ~ K + m + phi, subdf)
      
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
    fullmod <- lm(stat ~ sampstrat + K + m + phi, subdf)
    
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
    fullmod <- lm(stat ~ nsamp + K + m + phi, subdf)
    
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
