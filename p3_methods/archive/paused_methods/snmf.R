# run tess for all K values
tess3_obj <- tess3(X = gen[1:300,1:100], coord = as.matrix(coords[1:300,]), K = 1:30, ploidy = 2) 

# plot CV results and mark the K-value automatically selected  
plot(tess3_obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# get best K value
K <-  bestK(tess3_obj, 1:30)


res <- adegenet::find.clusters(gen[1:200,],  pca.select = "percVar", perc.pca = 90, choose.n.clust = FALSE, criterion = "diffNgroup", max.n.clust = 199)
K <- max(as.numeric(res$grp))

