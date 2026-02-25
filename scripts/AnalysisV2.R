library(dplyr)
library(geomorph)
library(geometry)
library(parallel)

set.seed(5)

# Setup -------------------------------------------------------------------

input <- "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/simulatedShapes/"

outputDisparity <- "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/analysisResults/disparity/"
outputMorphospace <- "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/analysisResults/morphospace/"

#parameters

#constatnt parameters
numReps = 100 #number of resims per condition
numTrees = 1000 #number of trees we analyzed
lddmmSigma = 1.0
betaDistAlpha = 1.0
betaDistBeta = 1.0

#variable parameters
alphas = c(0.1, 0.2, 0.3, 0.4)
numAdditionalTaxa = c(3, 5, 8, 11, 15, 20, 25, 50, 100, 144)
#numLandmarks = c(10, 25, 50)
numLandmarks = c(10)

#analysis functions

morphospaceVolume <- function(lmDF){
  #expecting lmDF with columns:
  #taxa name
  #landmark idx
  #x
  #y
  #z (if 3D)
  
  numLandmarks <- max(lmDF[,2])
  
  if(ncol(lmDF) == 4){
    #two dimensions
    dimension = 2
  }else if (ncol(lmDF) == 5){
    #three dimensions
    dimension = 3
  }else{
    stop("lmDF must have 4 columns (2D) or 5 columns (3D), got ", ncol(lmDF), " columns")
  }
  
  
  
  #rearrange
  taxa = unique(lmDF[,1])
  pcaMat <- matrix(data = NA, nrow = length(taxa), ncol = dimension * numLandmarks)
  rownames(pcaMat) = taxa
  for (i in taxa) {
    pcaMat[i, ] <- as.numeric(as.matrix(lmDF[lmDF[, 1] == i, 3:ncol(lmDF)]))
  }
  
  #PCA
  garray <- geomorph::arrayspecs(pcaMat, numLandmarks, k = dimension)
  gpa <- gpagen(garray, print.progress = FALSE, ProcD = TRUE)
  fullSpecimenPCA <- geomorph::gm.prcomp(gpa$coords)
  fullSpecimenPCs <- fullSpecimenPCA$x[,1:3] #just axes 1:3
  
  #chull volume including additional tips
  full_hull <- convhulln(fullSpecimenPCs, options = "FA")
  full_vol <- full_hull$vol
  
  #chull volume just including hominins
  subset_mat <- fullSpecimenPCs[!grepl("^Tip", rownames(fullSpecimenPCs)), , drop = FALSE]
  subset_hull <- convhulln(subset_mat, options = "FA")
  subset_vol <- subset_hull$vol
  
  resVec <- c(subset_vol, full_vol)
  names(resVec) <- c("homininSubsetChullVolume", "fullChullVolume")
  return(resVec)
}

# Morphological volume ----------------------------------------------------

#ETA tracker
totalIterations <- length(numLandmarks) * length(alphas) * length(numAdditionalTaxa) * numTrees
completedIterations <- 0
startTime <- Sys.time()

processTree <- function(
    numLand,
    lddmmAlpha,
    numAddTax, 
    treeIdx,
    rep
){
  
}

for(nlm in numLandmarks){
  for(a in alphas){
    for(nt in numAdditionalTaxa){
      
      #create output directory
      outdir <- paste0(
          outputMorphospace,
          "numLandmarks", nlm, "/",
          "numAdditionalTaxa", nt ,  "/" ,
          "alpha" , a , "0" ,  "/"
        )
      if(dir.exists(outdir) == FALSE){
        dir.create(
          outdir,
          recursive = T
        )
      }
      
      for(treeIdx in 0:(numTrees-1)){
          #create data storage obj
          resultsMatrix2D <- matrix(data = NA, nrow = numReps, ncol = 2)
          resultsMatrix3D <- matrix(data = NA, nrow = numReps, ncol = 2)
          
          colnames(resultsMatrix2D) <- c("homininSubsetChullVolume", "fullChullVolume")
          colnames(resultsMatrix3D) <- c("homininSubsetChullVolume", "fullChullVolume")
          
          for(rep in 0:(numReps-1)){
            #2D analysis
            shapeInputFile2D <- 
              paste0(
                input,
                "numLandmarks", nlm, "/",
                "numAdditionalTaxa", nt ,  "/" ,
                "alpha" , a , "0" ,  "/" ,
                "2D_TreeIndex" , treeIdx,
                "_TreeAlpha1.00_TreeBeta1.00_" , 
                "NumAdditionalTaxa", nt ,
                "_numLandmarks" , nlm ,
                "_lddmmAlpha" , a , "0" ,
                "_Sigma1.00",
                "_rep" , rep ,
                "nodeShapes.tsv.gz"
              )
            
            con <- gzfile(shapeInputFile2D)
            readFile <- read.delim(con, header = F)
            
            resultsMatrix2D[rep + 1, ] <- morphospaceVolume(readFile)
            
            #3D analysis
            shapeInputFile3D <- 
              paste0(
                input,
                "numLandmarks", nlm, "/",
                "numAdditionalTaxa", nt ,  "/" ,
                "alpha" , a , "0" ,  "/" ,
                "3D_TreeIndex" , treeIdx,
                "_TreeAlpha1.00_TreeBeta1.00_" , 
                "NumAdditionalTaxa", nt ,
                "_numLandmarks" , nlm ,
                "_lddmmAlpha" , a , "0" ,
                "_Sigma1.00",
                "_rep" , rep ,
                "nodeShapes.tsv.gz"
              )
            
            con <- gzfile(shapeInputFile3D)
            readFile <- read.delim(con, header = F)
            
            resultsMatrix3D[rep + 1, ] <- morphospaceVolume(readFile)
            rm(readFile)
          } 
          
          #write data storage obj to memory
          saveRDS(
            resultsMatrix2D,
            paste0(
              outdir,
              "2D_TreeIndex" , treeIdx,
              "_TreeAlpha1.00_TreeBeta1.00_" , 
              "NumAdditionalTaxa", nt ,
              "_numLandmarks" , nlm ,
              "_lddmmAlpha" , a , "0" ,
              "_Sigma1.00.rds"
            )
          )
          saveRDS(
            resultsMatrix3D,
            paste0(
              outdir,
              "3D_TreeIndex" , treeIdx,
              "_TreeAlpha1.00_TreeBeta1.00_" , 
              "NumAdditionalTaxa", nt ,
              "_numLandmarks" , nlm ,
              "_lddmmAlpha" , a , "0" ,
              "_Sigma1.00.rds"
            )
          )
          
          
          #ETA tracking
          completedIterations <- completedIterations + 1
          elapsed <- as.numeric(difftime(Sys.time(), startTime, units = "secs"))
          avgSecsPerIter <- elapsed / completedIterations
          remainingIters <- totalIterations - completedIterations
          etaSecs <- avgSecsPerIter * remainingIters
          etaTime <- Sys.time() + etaSecs
          
          message(sprintf(
            "[%d/%d] nlm=%s a=%s nt=%s treeIdx=%d | Elapsed: %.1f min | ETA: %s",
            completedIterations, totalIterations,
            nlm, a, nt, treeIdx,
            elapsed / 60,
            format(etaTime, "%Y-%m-%d %H:%M:%S")
          ))
        }
      
      gc()
    }
  }
}


# Full sim data analysis --------------------------------------------------


#files <- list.files("/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/simulatedShapes", full.names = T)
#saveRDS(files, "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/simulatedShapesToAnalyze.rds")
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
  start_disp <- Sys.time()
  dispFunc <- function(x){
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
    return(disparityPlotMat)
  }
  dispTaxa <- parallel::mclapply(1:numRepetitions, dispFunc, mc.cores = parallel::detectCores() - 1)
  saveRDS(dispTaxa, file = paste(
    outputDisparity,
    "disparityAnalysis",
    dimension, "Dimensions",
    treeIndex, "TreeIdx",
    alpha, "Alpha",
    beta, "Beta",
    numLandmarks, "Landmarks",
    rate, "Rate",
    dataset, "Dataset.rds",
    sep = ""
  ))
  end_disp <- Sys.time()
  cat("  Disparity completed in", round(difftime(end_disp, start_disp, units="mins"), 2), "minutes\n")

  
  ### morphospace volume over taxa
  # start_morpho <- Sys.time()
  # morphoFunc <- function(x){
  #   garray <- geomorph::arrayspecs(pcaMat, numLandmarks, k = dimension)
  #   gpa <- gpagen(garray, print.progress = FALSE, ProcD = TRUE)
  #   fullSpecimenPCA <- geomorph::gm.prcomp(gpa$coords)
  #   fullSpecimenPCs <- fullSpecimenPCA$x[,1:3] #just axes 1:3
  #   
  #   full_hull <- convhulln(fullSpecimenPCs, options = "FA")
  #   full_vol <- full_hull$vol
  #   
  #   subset_mat <- fullSpecimenPCs[!grepl("^Tip", rownames(fullSpecimenPCs)), , drop = FALSE]
  #   tempFullMat <- fullSpecimenPCs
  #   chullPlotMat <- matrix(data = NA, nrow = 160-16, ncol = 2)
  #   for (i in 1:(160-16)) {
  #     tip_rows <- grep("^Tip", rownames(tempFullMat), value = TRUE)
  #     sampled_row <- sample(tip_rows, 1)
  #     subset_mat <- rbind(subset_mat, tempFullMat[sampled_row, , drop = FALSE])
  #     tempFullMat <- tempFullMat[setdiff(rownames(tempFullMat), sampled_row), , drop = FALSE]
  #     
  #     subset_hull <- convhulln(subset_mat, options = "FA")
  #     subset_vol <- subset_hull$vol
  #     
  #     chullPlotMat[i, ] <- c(i + 16, subset_vol)
  #   }
  #   return(chullPlotMat)
  # }
  # morphoTaxa <- parallel::mclapply(1:numRepetitions, morphoFunc, mc.cores = parallel::detectCores() - 1)
  # saveRDS(morphoTaxa, file = paste(
  #   outputMorphospace,
  #   "morphospaceAnalysis",
  #   dimension, "Dimensions",
  #   treeIndex, "TreeIdx",
  #   alpha, "Alpha",
  #   beta, "Beta",
  #   numLandmarks, "Landmarks",
  #   rate, "Rate",
  #   dataset, "Dataset.rds",
  #   sep = ""
  # ))
  # end_morpho <- Sys.time()
  # cat("  Morphospace completed in", round(difftime(end_morpho, start_morpho, units="mins"), 2), "minutes\n")
  # 
  # ### remove file from files so that we avoid reanalysis if R crashes, etc. ###
  # files <- files[-1]
  # saveRDS(files, "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/simulatedShapesToAnalyze.rds")
  # 
  # end_file <- Sys.time()
  # elapsed_file <- as.numeric(difftime(end_file, start_file, units="mins"))
  # elapsed_total <- as.numeric(difftime(end_file, start_all, units="mins"))
  # avg_per_file <- elapsed_total / f_idx
  # remaining <- (total_files - f_idx) * avg_per_file
  # hours_remaining <- floor(remaining / 60)
  # minutes_remaining <- round(remaining %% 60)
  # 
  # cat(sprintf(
  #   "Finished file %d/%d in %.2f min | Avg per file: %.2f min | Time remaining: %d hr %d min\n\n",
  #   f_idx, total_files, elapsed_file, avg_per_file, hours_remaining, minutes_remaining
  # ))

    
}


# Hominin morpho volume baseline ------------------------------------------

files <- list.files("/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/simulatedShapes", full.names = T)

processHominins <- function(file){
    f <- file

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
    if (treeIndex > 99 || alpha != 1 || beta != 1 || numLandmarks != 10) {

    }else{
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
      
      garray <- geomorph::arrayspecs(pcaMat, numLandmarks, k = dimension)
      gpa <- gpagen(garray, print.progress = FALSE, ProcD = TRUE)
      fullSpecimenPCA <- geomorph::gm.prcomp(gpa$coords)
      fullSpecimenPCs <- fullSpecimenPCA$x[,1:3] #just axes 1:3
      
      subset_mat <- fullSpecimenPCs[!grepl("^Tip", rownames(fullSpecimenPCs)), , drop = FALSE]
      subset_hull <- convhulln(subset_mat, options = "FA")
      subset_vol <- subset_hull$vol
      
      row <- c(
        dimension,
        treeIndex,
        alpha,
        beta,
        numLandmarks,
        rate,
        dataset,
        subset_vol
      )
      
      names(row) <- c( "dimension", "treeIndex", "alpha", "beta", "numLandmarks", "rate", "dataset", "chullVol" )
      return(row)
    }
}

results <- parallel::mclapply(files,
                              processHominins, 
                              mc.cores = 6)
homininbaseline <- dplyr::bind_rows(results)
saveRDS(homininbaseline, "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/homininOnlyChullVolPaleoanthroSoc.rds")

# Morpho volume analysis -------------------------------------------------

### How much does each additional taxon contribute to the morphospace? ###


## For paleoanthro abstract: just 10 landmarks on uniform tree

## set up ##
outputMorphospace <- "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/analysisResults/morphospace/"
morphofiles <- list.files(outputMorphospace, full.names = T)
homininBaseline <- readRDS("/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/homininOnlyChullVolPaleoanthroSoc.rds")

## results objects ##
analysisFunc <- function(filePath){
  f <- filePath
  
  # --- extract metadata ---
  dimensions = as.numeric(sub(".*morphospaceAnalysis([0-9]+)Dimensions.*", "\\1", f))
  tree_index = as.numeric(sub(".*Dimensions([0-9]+)TreeIdx.*", "\\1", f))
  alpha      = as.numeric(sub(".*TreeIdx([0-9]+)Alpha.*", "\\1", f))
  beta       = as.numeric(sub(".*Alpha([0-9.]+)Beta.*", "\\1", f))
  landmarks  = as.numeric(sub(".*Beta([0-9.]+)Landmarks.*", "\\1", f))
  rate       = as.numeric(sub(".*Landmarks([0-9.]+)Rate.*", "\\1", f))
  dataset    = as.numeric(sub(".*Rate([0-9.]+)Dataset.*", "\\1", f))
  
  # skip unwanted files
  if (alpha != 1 || beta != 1 || landmarks != 10) {
    
  }else{
    homininChull <- dplyr::filter(homininBaseline,
                                  dimension==dimensions,
                                  treeIndex==tree_index,
                                  alpha==alpha,
                                  beta==beta,
                                  numLandmarks==landmarks)
    
    homininChull = homininChull[which(homininChull$rate == rate),]
    homininChull = homininChull[which(homininChull$dataset == dataset),]
    if(nrow(homininChull) != 1){
      stop("did not narrow down baseline chull")
    }
    homininChull = homininChull$chullVol
    
    # --- load file ---
    file <- readRDS(f)
    rawDF <- cbind("ntax" = file[[1]][,1])
    rawDF <- cbind(rawDF, "dimensions" = rep(dimensions, 144))
    rawDF <- cbind(rawDF, "treeIdx" = rep(tree_index, 144))
    rawDF <- cbind(rawDF, "rate" = rep(rate, 144))
    
    for(i in seq_along(file)){
      rawDF <- cbind(rawDF, file[[i]][,2])
    }
    
    #add in hominin baseline vol
    baselineHominin <- c(
      16,
      dimensions,
      tree_index,
      rate,
      rep(homininChull, 100)
    )
    rawDF <- rbind(baselineHominin, rawDF)
    
    ### normalize every column by the total volume of the morphospace ###
    for(i in 5:ncol(rawDF)){
      rawDF[,i] <- rawDF[,i] / max(rawDF[,i])
    }
    
    ### percent morphospace sampled ###
    
    # 1:25, 50, 100, and 144
    morphospaceRow <- c()
    for(i in c(1:25, 50, 100, 144)){
      dat <- rawDF[which(rawDF[,1] == i+16), 5:104]
      dat <- rawDF[1, 5:104] / dat
      morphospaceRow <- c(morphospaceRow, mean(dat))
    }
    names(morphospaceRow) <- c(1:25, 50, 100, 144)
    
    ### average contribution made by each taxon ###
    prawDF <- rawDF[,5:104]
    volGains <- c()
    for(i in 2:145){
      volGains <- c(
        volGains, 
        prawDF[i, ] - prawDF[i-1, ]
      )
    }
    averageVolumeGain <- mean(volGains)
    percent0VolumeGain <- sum(volGains == 0) / length(volGains)
    
    return(c(morphospaceRow, 
             "avPercentVolGain" = averageVolumeGain,
             "percent0VolumeGain" = percent0VolumeGain)) 
  }
}

res <- parallel::mclapply(morphofiles,
                          analysisFunc,
                          mc.cores = 10)
res2 <- dplyr::bind_rows(res)
saveRDS(res2, "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/homininMorphospaceResultsUpdate1to25incl.rds")

summary(res2)


# Analysis after resampling -----------------------------------------------

res = readRDS("/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/homininMorphospaceResultsUpdate1to25incl.rds")

