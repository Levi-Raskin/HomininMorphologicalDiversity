library(ggplot2)
library(ggtree)
library(phytools)

output <- "/Users/levir/Documents/GitHub/HomininTaxicDiversity/figures/"
outputMorphospace <- "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/analysisResults/morphospace/"


# num add tax, alpha = 0.1 ------------------------------------------------

numReps = 100 #number of resims per condition
numTrees = 1000 #number of trees we analyzed
lddmmSigma = 1.0
betaDistAlpha = 1.0
betaDistBeta = 1.0

#variable parameters
# alphas = c(0.1, 0.2, 0.3, 0.4)
alphas = c(0.1)
numAdditionalTaxa = c(3, 5, 8, 11, 15, 20, 25, 50, 100, 144)
#numLandmarks = c(10, 25, 50)
numLandmarks = c(10)

twoDimensionalPlotDat <- matrix(data = NA,
                  nrow = 0,
                  ncol = 4)
colnames(twoDimensionalPlotDat) <- c(
  "numLandmarks",
  "alpha",
  "numAdditionalTaxa",
  "percentRecovered"
)

threeDimensionalPlotDat <- matrix(data = NA,
                                nrow = 0,
                                ncol = 4)
colnames(threeDimensionalPlotDat) <- c(
  "numLandmarks",
  "alpha",
  "numAdditionalTaxa",
  "percentRecovered"
)

for(nlm in numLandmarks){
  for(a in alphas){
    for(nt in numAdditionalTaxa){
      outdir <- paste0(
        outputMorphospace,
        "numLandmarks", nlm, "/",
        "numAdditionalTaxa", nt ,  "/" ,
        "alpha" , a , "0" ,  "/"
      )
      
      numLandmarks <- rep(nlm, 100)
      lddmmAlpha <- rep(a, 100)
      numAddTax <- rep(nt, 100)
      
      for(treeIdx in 0:(numTrees-1)){
        
        twoD <- readRDS(
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
        
        temp <- cbind(
          numLandmarks,
          lddmmAlpha,
          numAddTax,
          twoD[,1] / twoD[,2]
        )
        
        twoDimensionalPlotDat <- rbind(twoDimensionalPlotDat, temp)
        
        if(any(is.na(twoD))){
          print(paste0(
            "NA HERE: ",
            outdir,
            "2D_TreeIndex" , treeIdx,
            "_TreeAlpha1.00_TreeBeta1.00_" , 
            "NumAdditionalTaxa", nt ,
            "_numLandmarks" , nlm ,
            "_lddmmAlpha" , a , "0" ,
            "_Sigma1.00.rds"
          ))
        }
        
        threeD <- readRDS(
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
        
        if(any(is.na(threeD))){
          print(paste0(
            "NA HERE: ",
            outdir,
            "2D_TreeIndex" , treeIdx,
            "_TreeAlpha1.00_TreeBeta1.00_" , 
            "NumAdditionalTaxa", nt ,
            "_numLandmarks" , nlm ,
            "_lddmmAlpha" , a , "0" ,
            "_Sigma1.00.rds"
          ))
        }
        
        temp <- cbind(
          numLandmarks,
          lddmmAlpha,
          numAddTax,
          threeD[,1] / threeD[,2]
        )
        
        threeDimensionalPlotDat <- rbind(threeDimensionalPlotDat, temp)
        
      }
      
      print(paste(
        "a: ", a,
        "nt: ", nt,
        "nlm: ", nlm
      ))
    }
  }
}

twoDimensionalPlotDat <- data.frame(twoDimensionalPlotDat)
threeDimensionalPlotDat <- data.frame(threeDimensionalPlotDat)

