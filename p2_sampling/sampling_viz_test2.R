#code to create some plots of example subsamples
source("general_functions.R")
pdf("gsd_plots.pdf")
for(i in sample(1:nrow(params),5)){
    gsd_df <- get_data(i, params = params, "gsd")
    plot(gsd_df$x, gsd_df$y, xlim = c(0,ldim), ylim = c(0,ldim),
             xlab = "", ylab = "", pch = 19, col = "#756fd6", cex = 0.1)
  }
dev.off()


