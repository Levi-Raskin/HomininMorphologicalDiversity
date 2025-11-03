library(geomorph)
landmarks <- read.delim("/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/simulatedShapes/twoDimensionSimTreeIndex6Alpha1.000000Beta1.000000LM25Rate0.300000Dataset0nodeShapes.tsv", header = F)
cn <- colnames(landmarks)
cn[1] <- "label"
cn[2] <- "lm"
colnames(landmarks) <- cn

numLandmarks = 25
dimension = 2

pcaMat <- matrix(data = NA, nrow = length(unique(landmarks$label)), ncol = dimension * numLandmarks)
rownames(pcaMat) = unique(landmarks$label)
for(i in unique(landmarks$label)){
  lm <- as.matrix(landmarks[which(landmarks[,1] == i), ])
  lm <- c(lm[, 3:ncol(lm)])
  for(j in 1:length(lm)){
    pcaMat[which(rownames(pcaMat)== i), j] = as.numeric(lm[j]) 
  }
}
### gpagen expects an array containing n x d slices where n is num landmarks and d is dimension of ladnmarks


subset_mat <- pcaMat[!grepl("^Tip", rownames(pcaMat)), , drop = FALSE]
tempFullMat <- pcaMat
disparityPlotMat <- matrix(data = NA, nrow = 160-16, ncol = 2)
for (i in 1:(160-16)) {
  tip_rows <- grep("^Tip", rownames(tempFullMat), value = TRUE)
  sampled_row <- sample(tip_rows, 1)
  subset_mat <- rbind(subset_mat, tempFullMat[sampled_row, , drop = FALSE])
  tempFullMat <- tempFullMat[setdiff(rownames(tempFullMat), sampled_row), , drop = FALSE]
  
  garray <- geomorph::arrayspecs(subset_mat, numLandmarks, k = dimension)
  gpa <- gpagen(garray, print.progress = FALSE)
  gdf <- geomorph.data.frame(gpa)
  
  
  
  disparityPlotMat[i, ] <- c(i + 16, as.numeric(morphol.disparity(coords ~ 1, groups = NULL, data = gdf, print.progress = FALSE)) )
}

plot(disparityPlotMat)
