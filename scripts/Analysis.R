library(geomorph)
library(geometry)
library(parallel)

set.seed(5)

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

### total morpho volume ###
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

### volume gain per taxon ###
garray <- geomorph::arrayspecs(pcaMat, numLandmarks, k = dimension)
gpa <- gpagen(garray, print.progress = FALSE, ProcD = TRUE)
fullSpecimenPCA <- geomorph::gm.prcomp(gpa$coords)
fullSpecimenPCs <- fullSpecimenPCA$x[, 1:3] # just axes 1:3

full_hull <- convhulln(fullSpecimenPCs, options = "FA")
full_vol <- full_hull$vol

subset_mat <- fullSpecimenPCs[!grepl("^Tip", rownames(fullSpecimenPCs)), , drop = FALSE]
tempFullMat <- fullSpecimenPCs
chullPlotMat <- matrix(NA, nrow = 160 - 16, ncol = 2)

subset_hull <- convhulln(subset_mat, options = "FA")
prev_vol <- subset_hull$vol  # store previous total volume

for (i in 1:(160 - 16)) {
  tip_rows <- grep("^Tip", rownames(tempFullMat), value = TRUE)
  sampled_row <- sample(tip_rows, 1)
  subset_mat <- rbind(subset_mat, tempFullMat[sampled_row, , drop = FALSE])
  tempFullMat <- tempFullMat[setdiff(rownames(tempFullMat), sampled_row), , drop = FALSE]
  
  subset_hull <- convhulln(subset_mat, options = "FA")
  subset_vol <- subset_hull$vol
  
  vol_gain <- subset_vol - prev_vol  # difference from last total volume
  chullPlotMat[i, ] <- c(i + 16, vol_gain)
  
  prev_vol <- subset_vol  # update for next iteration
}

colnames(chullPlotMat) <- c("num taxa", "volume gained with each add. taxon")
plot(chullPlotMat)


# Full sim data analysis --------------------------------------------------

numRepetitions <- 100 #number of disparty/morpho variance vs. num taxa subsamplings we perform
outputDisparity <- "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/analysisResults/disparity/"
outputMorphospace <- "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/analysisResults/morphospace/"

# files <- list.files("/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/simulatedShapes", full.names = T)
# saveRDS(files, "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/simulatedShapesToAnalyze.rds")
files <- readRDS("/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/simulatedShapesToAnalyze.rds")

total_files <- length(files)
start_all <- Sys.time()

cat("Starting analysis of", total_files, "files at", format(start_all), "\n\n")

while (length(files) > 0) {
  f <- files[1]
  start_file <- Sys.time()
  f_idx <- total_files - length(files) + 1
  
  cat(sprintf("[%d/%d] Processing file: %s\n", f_idx, total_files, basename(f)))
  
  ### regex and string parsing ###
  word_to_num <- c(
    one=1, two=2, three=3, four=4, five=5,
    six=6, seven=7, eight=8, nine=9, ten=10
  )
  rx <- paste0(
    "(one|two|three|four|five|six|seven|eight|nine|ten)Dimension",
    "SimTreeIndex(\\d+)",
    "Alpha([-+0-9\\.eE]+)",
    "Beta([-+0-9\\.eE]+)",
    "LM(\\d+)",
    "Rate([-+0-9\\.eE]+)",
    "Dataset(\\d+)"
  )
  m <- regexec(rx, f)
  got <- regmatches(f, m)
  mat <- do.call(rbind, lapply(got, function(x) if (length(x)) x[-1] else rep(NA_character_, 7)))
  colnames(mat) <- c("dimension_word","tree_index","alpha","beta","lm","rate","dataset")

  dimension = unname(word_to_num[mat[,"dimension_word"]])
  treeIndex = as.integer(mat[,"tree_index"])
  alpha = as.numeric(mat[,"alpha"])
  beta = as.numeric(mat[,"beta"])
  numLandmarks = as.integer(mat[,"lm"])
  rate = as.numeric(mat[,"rate"])
  dataset = as.integer(mat[,"dataset"])
  
  ### only analyze the first 99 of these for paleoanthro abstract ###
  if (treeIndex > 99) {
    cat(sprintf("[%d/%d] Skipping file %s (TreeIdx > 99)\n",
                f_idx, total_files, basename(f)))
    files <- files[-1]
    saveRDS(files, "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/simulatedShapesToAnalyze.rds")
    next
  }
  
  ### data read in and set up
  landmarks <- read.delim(f, header = F)
  cn <- colnames(landmarks)
  cn[1] <- "label"
  cn[2] <- "lm"
  colnames(landmarks) <- cn
  
  pcaMat <- matrix(data = NA, nrow = length(unique(landmarks$label)), ncol = dimension * numLandmarks)
  rownames(pcaMat) = unique(landmarks$label)
  for(i in unique(landmarks$label)){
    lm <- as.matrix(landmarks[which(landmarks[,1] == i), ])
    lm <- c(lm[, 3:ncol(lm)])
    for(j in 1:length(lm)){
      pcaMat[which(rownames(pcaMat)== i), j] = as.numeric(lm[j]) 
    }
  }
  
  ### disparity through taxa ###
  ## disparity through taxa analysis is super expensive; commenting out ###
  # start_disp <- Sys.time()
  # dispFunc <- function(x){
  #   subset_mat <- pcaMat[!grepl("^Tip", rownames(pcaMat)), , drop = FALSE]
  #   tempFullMat <- pcaMat
  #   disparityPlotMat <- matrix(data = NA, nrow = 160-16, ncol = 2)
  #   for (i in 1:(160-16)) {
  #     tip_rows <- grep("^Tip", rownames(tempFullMat), value = TRUE)
  #     sampled_row <- sample(tip_rows, 1)
  #     subset_mat <- rbind(subset_mat, tempFullMat[sampled_row, , drop = FALSE])
  #     tempFullMat <- tempFullMat[setdiff(rownames(tempFullMat), sampled_row), , drop = FALSE]
  #     
  #     garray <- geomorph::arrayspecs(subset_mat, numLandmarks, k = dimension)
  #     gpa <- gpagen(garray, print.progress = FALSE, ProcD = TRUE)
  #     gdf <- geomorph.data.frame(gpa)
  #     
  #     disparityPlotMat[i, ] <- c(i + 16, as.numeric(morphol.disparity(coords ~ 1, groups = NULL, data = gdf, print.progress = FALSE)) )
  #   }
  #   return(disparityPlotMat)
  # }
  # dispTaxa <- parallel::mclapply(1:numRepetitions, dispFunc, mc.cores = parallel::detectCores() - 1)
  # saveRDS(dispTaxa, file = paste(
  #   outputDisparity,
  #   "disparityAnalysis",
  #   dimension, "Dimensions",
  #   treeIndex, "TreeIdx",
  #   alpha, "Alpha",
  #   beta, "Beta",
  #   numLandmarks, "Landmarks",
  #   rate, "Rate",
  #   dataset, "Dataset.rds",
  #   sep = ""
  # ))
  # end_disp <- Sys.time()
  # cat("  Disparity completed in", round(difftime(end_disp, start_disp, units="mins"), 2), "minutes\n")
  
  
  ### morphospace volume over taxa
  start_morpho <- Sys.time()
  morphoFunc <- function(x){
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
    return(chullPlotMat)
  }
  morphoTaxa <- parallel::mclapply(1:numRepetitions, morphoFunc, mc.cores = parallel::detectCores() - 1)
  saveRDS(morphoTaxa, file = paste(
    outputMorphospace,
    "morphospaceAnalysis",
    dimension, "Dimensions",
    treeIndex, "TreeIdx",
    alpha, "Alpha",
    beta, "Beta",
    numLandmarks, "Landmarks",
    rate, "Rate",
    dataset, "Dataset.rds",
    sep = ""
  ))
  end_morpho <- Sys.time()
  cat("  Morphospace completed in", round(difftime(end_morpho, start_morpho, units="mins"), 2), "minutes\n")
  
  ### remove file from files so that we avoid reanalysis if R crashes, etc. ###
  files <- files[-1]
  saveRDS(files, "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/simulatedShapesToAnalyze.rds")
  
  end_file <- Sys.time()
  elapsed_file <- as.numeric(difftime(end_file, start_file, units="mins"))
  elapsed_total <- as.numeric(difftime(end_file, start_all, units="mins"))
  avg_per_file <- elapsed_total / f_idx
  remaining <- (total_files - f_idx) * avg_per_file
  hours_remaining <- floor(remaining / 60)
  minutes_remaining <- round(remaining %% 60)
  
  cat(sprintf(
    "Finished file %d/%d in %.2f min | Avg per file: %.2f min | Time remaining: %d hr %d min\n\n",
    f_idx, total_files, elapsed_file, avg_per_file, hours_remaining, minutes_remaining
  ))

    
}

