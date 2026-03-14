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

processData2D <- function(
    numLand,
    lddmmAlpha,
    numAddTax, 
    tIdx,
    r
){
  #2D analysis
  shapeInputFile2D <- 
    paste0(
      input,
      "numLandmarks", numLand, "/",
      "numAdditionalTaxa", numAddTax ,  "/" ,
      "alpha" , lddmmAlpha , "0" ,  "/" ,
      "2D_TreeIndex" , tIdx,
      "_TreeAlpha1.00_TreeBeta1.00_" , 
      "NumAdditionalTaxa", numAddTax ,
      "_numLandmarks" , numLand ,
      "_lddmmAlpha" , lddmmAlpha , "0" ,
      "_Sigma1.00",
      "_rep" , r ,
      "nodeShapes.tsv.gz"
    )
  
  con <- gzfile(shapeInputFile2D)
  readFile <- read.delim(con, header = F)
  
  return(morphospaceVolume(readFile))
}

processData3D <- function(
    numLand,
    lddmmAlpha,
    numAddTax, 
    tIdx,
    r
){
  #3D analysis
  shapeInputFile3D <- 
    paste0(
      input,
      "numLandmarks", nlm, "/",
      "numAdditionalTaxa", nt ,  "/" ,
      "alpha" , a , "0" ,  "/" ,
      "3D_TreeIndex" , tIdx,
      "_TreeAlpha1.00_TreeBeta1.00_" , 
      "NumAdditionalTaxa", nt ,
      "_numLandmarks" , nlm ,
      "_lddmmAlpha" , a , "0" ,
      "_Sigma1.00",
      "_rep" , r ,
      "nodeShapes.tsv.gz"
    )
  
  con <- gzfile(shapeInputFile3D)
  readFile <- read.delim(con, header = F)
  
  return(morphospaceVolume(readFile))
}

# alphas = c(0.1, 0.2, 0.3, 0.4)
# numAdditionalTaxa = c(3, 5, 8, 11, 15, 20, 25, 50, 100, 144)
numLandmarks = c(10, 25, 50)

totalIter <- 0
for(nlm in numLandmarks){
  if(nlm == 10){ alphas = c(0.3, 0.4) } else { alphas = c(0.1, 0.2, 0.3, 0.4) }
  for(a in alphas){
    if(nlm == 10 && a == 0.3){ numAdditionalTaxa = c(144) } else { numAdditionalTaxa = c(3, 5, 8, 11, 15, 20, 25, 50, 100, 144) }
    totalIter <- totalIter + length(numAdditionalTaxa) * numTrees
  }
}

completedIter <- 0
startTime <- Sys.time()

for(nlm in numLandmarks){
  
  if(nlm == 10){
    alphas = c(0.3, 0.4)
  }else{
    alphas = c(0.1, 0.2, 0.3, 0.4)
  }
  
  for(a in alphas){
    if(nlm == 10 && a == 0.3){
      numAdditionalTaxa = c(144)
    }else{
      numAdditionalTaxa = c(3, 5, 8, 11, 15, 20, 25, 50, 100, 144)
    }
    
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
          resultsMatrix2D <- mclapply(
            X = 0:(numReps-1),
            FUN = function(rep){
                processData2D(
                  nlm,
                  a,
                  nt, 
                  treeIdx,
                  rep
                )
              }
            ,
            mc.cores = 12
          )
          resultsMatrix2D <- as.data.frame(do.call(rbind, resultsMatrix2D))

          resultsMatrix3D <- mclapply(
            X = 0:(numReps-1),
            FUN = function(rep){
              processData3D(
                nlm,
                a,
                nt, 
                treeIdx,
                rep
              )
            }
            ,
            mc.cores = 12
          )
          
          resultsMatrix3D <- as.data.frame(do.call(rbind, resultsMatrix3D))
          
         
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
          
          
          # ETA tracking
          completedIter <- completedIter + 1
          elapsed <- as.numeric(difftime(Sys.time(), startTime, units = "secs"))
          avgSecPerIter <- elapsed / completedIter
          remainingIter <- totalIter - completedIter
          etaSecs <- avgSecPerIter * remainingIter
          etaTime <- Sys.time() + etaSecs
          
          message(sprintf(
            "[%d/%d | %.1f%%] nlm=%s | alpha=%.1f | nt=%d | tree=%d | Elapsed: %s | ETA: %s (~%s remaining)",
            completedIter, totalIter,
            100 * completedIter / totalIter,
            nlm, a, nt, treeIdx,
            format(round(as.difftime(elapsed, units="secs")), format="%H:%M:%S"),
            format(etaTime, "%H:%M:%S"),
            format(round(as.difftime(etaSecs, units="secs")), format="%H:%M:%S")
          ))
          
          
          
        }
      
      gc()
    }
    
  }
}