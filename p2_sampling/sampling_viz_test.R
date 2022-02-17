#code to create some plots of example subsamples
source("general_functions.R")


pdf("sampling_plots.pdf")
par(mfrow=c(length(npts), length(sampstrats)), pty="s", mar = rep(2, 4))
for(i in sample(1:nrow(params), 10)){
    gsd_df <- get_data(i, params = params, "gsd")
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params = params, sampstrat, nsamp)
        subgsd_df <- gsd_df[subIDs,]
        plot(subgsd_df$x, subgsd_df$y, xlim = c(0,ldim), ylim = c(0,ldim), main = paste0(nsamp," / ", sampstrat),
             xlab = "", ylab = "", pch = 19, col = "#756fd6")
      }
    }
}
dev.off()

source("site_functions.R")
pdf("site_sampling_plots.pdf")
par(mfrow=c(length(nsites), length(sampstrats)), pty="s", mar = rep(2, 4))
for(i in sample(1:nrow(params),5)){
    gsd_df <- get_data(i, params = params, "gsd")
    for(nsite in nsites){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <-  get_samples(params[i,], params, sampstrat, nsite)
        subgsd_df <- gsd_df[subIDs,]
        plot(subgsd_df$x, subgsd_df$y, xlim = c(0,ldim), ylim = c(0,ldim), main = paste0(nsite," / ", sampstrat),
             xlab = "", ylab = "", pch = 19, col = "#756fd6")
      }
    }
  }
dev.off()


