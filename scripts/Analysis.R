library(geomorph)
library(geometry)
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



# Disparity ---------------------------------------------------------------
subset_mat <- pcaMat[!grepl("^Tip", rownames(pcaMat)), , drop = FALSE]
tempFullMat <- pcaMat
disparityPlotMat <- matrix(data = NA, nrow = 160-16, ncol = 2)
for (i in 1:(160-16)) {
  tip_rows <- grep("^Tip", rownames(tempFullMat), value = TRUE)
  sampled_row <- sample(tip_rows, 1)
  subset_mat <- rbind(subset_mat, tempFullMat[sampled_row, , drop = FALSE])
  tempFullMat <- tempFullMat[setdiff(rownames(tempFullMat), sampled_row), , drop = FALSE]
  
  garray <- geomorph::arrayspecs(subset_mat, numLandmarks, k = dimension)
  gpa <- gpagen(garray, print.progress = FALSE, ProcD = TRUE)
  gdf <- geomorph.data.frame(gpa)

  disparityPlotMat[i, ] <- c(i + 16, as.numeric(morphol.disparity(coords ~ 1, groups = NULL, data = gdf, print.progress = FALSE)) )
}

colnames(disparityPlotMat) <- c("num taxa", "disparity")
plot(disparityPlotMat)

# Chull volume ------------------------------------------------------------

garray <- geomorph::arrayspecs(pcaMat, numLandmarks, k = dimension)
gpa <- gpagen(garray, print.progress = FALSE, ProcD = TRUE)
fullSpecimenPCA <- geomorph::gm.prcomp(gpa$coords)
fullSpecimenPCs <- fullSpecimenPCA$x[,1:3] #just axes 1:3

full_hull <- convhulln(fullSpecimenPCs, options = "FA")
full_vol <- full_hull$vol

subset_mat <- fullSpecimenPCs[!grepl("^Tip", rownames(fullSpecimenPCs)), , drop = FALSE]
tempFullMat <- fullSpecimenPCs
chullPlotMat <- matrix(data = NA, nrow = 160-16, ncol = 2)
for (i in 1:(160-16)) {
  tip_rows <- grep("^Tip", rownames(tempFullMat), value = TRUE)
  sampled_row <- sample(tip_rows, 1)
  subset_mat <- rbind(subset_mat, tempFullMat[sampled_row, , drop = FALSE])
  tempFullMat <- tempFullMat[setdiff(rownames(tempFullMat), sampled_row), , drop = FALSE]
  
  subset_hull <- convhulln(subset_mat, options = "FA")
  subset_vol <- subset_hull$vol

    chullPlotMat[i, ] <- c(i + 16, subset_vol)
}

colnames(chullPlotMat) <- c("num taxa", "subset hull volume")
plot(chullPlotMat)



