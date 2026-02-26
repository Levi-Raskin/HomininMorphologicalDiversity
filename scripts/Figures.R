library(dplyr)
library(gghalves)
library(ggnewscale)
library(ggplot2)
library(ggtree)
library(phytools)
library(RColorBrewer)

output <- "/Users/levir/Documents/GitHub/HomininTaxicDiversity/figures/"
outputMorphospace <- "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/analysisResults/morphospace/"


# num add tax, alpha = 0.1 ------------------------------------------------

numReps = 100 #number of resims per condition
numTrees = 1000 #number of trees we analyzed
lddmmSigma = 1.0
betaDistAlpha = 1.0
betaDistBeta = 1.0

#variable parameters
#alphas = c(0.1, 0.2, 0.3, 0.4)
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
        
        twoD = dplyr::mutate(twoD, perc = twoD$homininSubsetChullVolume /twoD$fullChullVolume )
        
        temp <- cbind(
          numLandmarks,
          lddmmAlpha,
          numAddTax,
          twoD$perc
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
        threeD$homininSubsetChullVolume <- as.numeric(threeD$homininSubsetChullVolume)
        threeD$fullChullVolume <- as.numeric(threeD$fullChullVolume)
        threeD = dplyr::mutate(threeD, perc = threeD$homininSubsetChullVolume /threeD$fullChullVolume  )
        
        temp <- cbind(
          numLandmarks,
          lddmmAlpha,
          numAddTax,
          threeD$perc
        )
        
        threeDimensionalPlotDat <- rbind(threeDimensionalPlotDat, temp)
        
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

twoDimensionalPlotDat$numAdditionalTaxa <- as.factor(twoDimensionalPlotDat$numAdditionalTaxa)
threeDimensionalPlotDat$numAdditionalTaxa <- as.factor(threeDimensionalPlotDat$numAdditionalTaxa)

saveRDS(object = twoDimensionalPlotDat,
        file = paste0(output, "twoDimensional_LM10_rate0.1.rds"))
saveRDS(object = threeDimensionalPlotDat,
        file = paste0(output, "threeDimensional_LM10_rate0.1.rds"))

nullFunc = function(i){
  return (17 / (17 + as.numeric(as.character(i))))
}

null_dat <- data.frame(
  numAdditionalTaxa = unique(twoDimensionalPlotDat$numAdditionalTaxa)
) %>%
  mutate(null_val = nullFunc(numAdditionalTaxa))


# two dimensions
mean_dat <- twoDimensionalPlotDat %>%
  group_by(numAdditionalTaxa) %>%
  summarise(
    mean_val = mean(percentRecovered),
    lower_95 = quantile(percentRecovered, 0.025),
    upper_95 = quantile(percentRecovered, 0.975)
  )

p1 <- ggplot() +
  geom_half_violin(data = twoDimensionalPlotDat,
                   aes(x = numAdditionalTaxa, y = percentRecovered), 
                   fill = "#3b7ea1", 
                   side = "l") +
  geom_point(data = mean_dat,
             aes(x = numAdditionalTaxa, y = mean_val),
             color = "black", size = 2.5, shape = 18) +
  geom_text(data = mean_dat,
            aes(x = numAdditionalTaxa, y = mean_val,
                label = paste0("Mean: ", round(mean_val, 2))),
            hjust = -0.1, color = "black", size = 3, fontface = "bold",
            family = "Georgia") +
  geom_text(data = mean_dat,
            aes(x = numAdditionalTaxa, y = lower_95,
                label = paste0("2.5%: ", round(lower_95, 2))),
            hjust = -0.1, color = "black", size = 2.8,
            family = "Georgia") +
  geom_text(data = mean_dat,
            aes(x = numAdditionalTaxa, y = upper_95,
                label = paste0("97.5%: ", round(upper_95, 2))),
            hjust = -0.1, color = "black", size = 2.8,
            family = "Georgia") +
  geom_text(data = null_dat,
            aes(x = numAdditionalTaxa, y = null_val,
                label = paste0("Ex: ", round(null_val, 2))),
            hjust = -0.1, color = "black", size = 2.8,
            family = "Georgia") +
  xlab("Number of additional taxa") +
  ylab("Percent recovered disparity") +
  theme_minimal(base_family = "Georgia") +
  scale_x_discrete(expand = expansion(add = c(0, 1.0)))
p1

ggsave(
  paste0(output, "2Dalpha0.1.svg"),
  p1,
  height = 3.5,
  width = 11
)

#three dimensions
mean_dat <- threeDimensionalPlotDat %>%
  group_by(numAdditionalTaxa) %>%
  summarise(
    mean_val = mean(percentRecovered, na.rm = T),
    lower_95 = quantile(percentRecovered, 0.025, na.rm = T),
    upper_95 = quantile(percentRecovered, 0.975, na.rm = T)
  )

p1 <- ggplot() +
  geom_half_violin(data = threeDimensionalPlotDat,
                   aes(x = numAdditionalTaxa, y = percentRecovered), 
                   fill = "#3b7ea1", 
                   side = "l") +
  geom_point(data = mean_dat,
             aes(x = numAdditionalTaxa, y = mean_val),
             color = "black", size = 2.5, shape = 18) +
  geom_text(data = mean_dat,
            aes(x = numAdditionalTaxa, y = mean_val,
                label = paste0("Mean: ", round(mean_val, 2))),
            hjust = -0.1, color = "black", size = 3, fontface = "bold",
            family = "Georgia") +
  geom_text(data = mean_dat,
            aes(x = numAdditionalTaxa, y = lower_95,
                label = paste0("2.5%: ", round(lower_95, 2))),
            hjust = -0.1, color = "black", size = 2.8,
            family = "Georgia") +
  geom_text(data = mean_dat,
            aes(x = numAdditionalTaxa, y = upper_95,
                label = paste0("97.5%: ", round(upper_95, 2))),
            hjust = -0.1, color = "black", size = 2.8,
            family = "Georgia") +
  geom_text(data = null_dat,
            aes(x = numAdditionalTaxa, y = null_val,
                label = paste0("Ex: ", round(null_val, 2))),
            hjust = -0.1, color = "black", size = 2.8,
            family = "Georgia") +
  xlab("Number of additional taxa") +
  ylab("Percent recovered disparity") +
  theme_minimal(base_family = "Georgia") +
  scale_x_discrete(expand = expansion(add = c(0, 1.0)))
p1

ggsave(
  paste0(output, "3Dalpha0.1.svg"),
  p1,
  height = 3.5,
  width = 11
)


#twoD and threeD
mean_dat2d <- twoDimensionalPlotDat %>%
  group_by(numAdditionalTaxa) %>%
  summarise(
    mean_val = mean(percentRecovered, na.rm = T)
  )
mean_dat3d <- threeDimensionalPlotDat %>%
  group_by(numAdditionalTaxa) %>%
  summarise(
    mean_val = mean(percentRecovered, na.rm = T)
  )

p1 <- ggplot() +
  geom_half_violin(data = threeDimensionalPlotDat,
                   aes(x = numAdditionalTaxa, y = percentRecovered), 
                   fill = "#215834", 
                   side = "r") +
  geom_half_violin(data = twoDimensionalPlotDat,
                   aes(x = numAdditionalTaxa, y = percentRecovered), 
                   fill = "#3b7ea1", 
                   side = "l") +
  geom_point(data = mean_dat3d,
             aes(x = numAdditionalTaxa, y = mean_val),
             color = "black", size = 2.5, shape = 18) +
  geom_text(data = mean_dat3d,
            aes(x = numAdditionalTaxa, y = mean_val,
                label = paste0("3D: ", round(mean_val, 2))),
            hjust = -0.1, color = "black", size = 3, fontface = "bold",
            family = "Georgia") +
  geom_point(data = mean_dat2d,
             aes(x = numAdditionalTaxa, y = mean_val),
             color = "black", size = 2.5, shape = 18) +
  geom_text(data = mean_dat2d,
            aes(x = numAdditionalTaxa, y = mean_val,
                label = paste0("2D: ", round(mean_val, 2))),
            hjust = 1.1, color = "black", size = 3, fontface = "bold",
            family = "Georgia") +
  xlab("Number of additional taxa") +
  ylab("Percent recovered disparity") +
  theme_minimal(base_family = "Georgia") +
  scale_x_discrete(expand = expansion(add = c(0.75, 0.75)))
p1

ggsave(
  paste0(output, "twoDthreeDCompared.svg"),
  p1,
  height = 3.5,
  width = 11
)

for(i in unique(twoDimensionalPlotDat$numAdditionalTaxa)){
  print(wilcox.test(filter(twoDimensionalPlotDat, twoDimensionalPlotDat$numAdditionalTaxa ==i)$percentRecovered, filter(threeDimensionalPlotDat, threeDimensionalPlotDat$numAdditionalTaxa ==i)$percentRecovered))
}



# num add tax, alpha = 0.2 ------------------------------------------------

numReps = 100 #number of resims per condition
numTrees = 1000 #number of trees we analyzed
lddmmSigma = 1.0
betaDistAlpha = 1.0
betaDistBeta = 1.0

#variable parameters
#alphas = c(0.1, 0.2, 0.3, 0.4)
alphas = c(0.2)
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
        
        twoD = dplyr::mutate(twoD, perc = twoD$homininSubsetChullVolume /twoD$fullChullVolume )
        
        temp <- cbind(
          numLandmarks,
          lddmmAlpha,
          numAddTax,
          twoD$perc
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
        threeD$homininSubsetChullVolume <- as.numeric(threeD$homininSubsetChullVolume)
        threeD$fullChullVolume <- as.numeric(threeD$fullChullVolume)
        threeD = dplyr::mutate(threeD, perc = threeD$homininSubsetChullVolume /threeD$fullChullVolume  )
        
        temp <- cbind(
          numLandmarks,
          lddmmAlpha,
          numAddTax,
          threeD$perc
        )
        
        threeDimensionalPlotDat <- rbind(threeDimensionalPlotDat, temp)
        
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

twoDimensionalPlotDat$numAdditionalTaxa <- as.factor(twoDimensionalPlotDat$numAdditionalTaxa)
threeDimensionalPlotDat$numAdditionalTaxa <- as.factor(threeDimensionalPlotDat$numAdditionalTaxa)

saveRDS(object = twoDimensionalPlotDat,
        file = paste0(output, "twoDimensional_LM10_rate0.2.rds"))
saveRDS(object = threeDimensionalPlotDat,
        file = paste0(output, "threeDimensional_LM10_rate0.2.rds"))

nullFunc = function(i){
  return (17 / (17 + as.numeric(as.character(i))))
}

null_dat <- data.frame(
  numAdditionalTaxa = unique(twoDimensionalPlotDat$numAdditionalTaxa)
) %>%
  mutate(null_val = nullFunc(numAdditionalTaxa))


# two dimensions
mean_dat <- twoDimensionalPlotDat %>%
  group_by(numAdditionalTaxa) %>%
  summarise(
    mean_val = mean(percentRecovered),
    lower_95 = quantile(percentRecovered, 0.025),
    upper_95 = quantile(percentRecovered, 0.975)
  )

p1 <- ggplot() +
  geom_half_violin(data = twoDimensionalPlotDat,
                   aes(x = numAdditionalTaxa, y = percentRecovered), 
                   fill = "#3b7ea1", 
                   side = "l") +
  geom_point(data = mean_dat,
             aes(x = numAdditionalTaxa, y = mean_val),
             color = "black", size = 2.5, shape = 18) +
  geom_text(data = mean_dat,
            aes(x = numAdditionalTaxa, y = mean_val,
                label = paste0("Mean: ", round(mean_val, 2))),
            hjust = -0.1, color = "black", size = 3, fontface = "bold",
            family = "Georgia") +
  geom_text(data = mean_dat,
            aes(x = numAdditionalTaxa, y = lower_95,
                label = paste0("2.5%: ", round(lower_95, 2))),
            hjust = -0.1, color = "black", size = 2.8,
            family = "Georgia") +
  geom_text(data = mean_dat,
            aes(x = numAdditionalTaxa, y = upper_95,
                label = paste0("97.5%: ", round(upper_95, 2))),
            hjust = -0.1, color = "black", size = 2.8,
            family = "Georgia") +
  geom_text(data = null_dat,
            aes(x = numAdditionalTaxa, y = null_val,
                label = paste0("Ex: ", round(null_val, 2))),
            hjust = -0.1, color = "black", size = 2.8,
            family = "Georgia") +
  xlab("Number of additional taxa") +
  ylab("Percent recovered disparity") +
  theme_minimal(base_family = "Georgia") +
  scale_x_discrete(expand = expansion(add = c(0, 1.0)))
p1

ggsave(
  paste0(output, "2Dalpha0.1.svg"),
  p1,
  height = 3.5,
  width = 11
)

#three dimensions
mean_dat <- threeDimensionalPlotDat %>%
  group_by(numAdditionalTaxa) %>%
  summarise(
    mean_val = mean(percentRecovered, na.rm = T),
    lower_95 = quantile(percentRecovered, 0.025, na.rm = T),
    upper_95 = quantile(percentRecovered, 0.975, na.rm = T)
  )

p1 <- ggplot() +
  geom_half_violin(data = threeDimensionalPlotDat,
                   aes(x = numAdditionalTaxa, y = percentRecovered), 
                   fill = "#3b7ea1", 
                   side = "l") +
  geom_point(data = mean_dat,
             aes(x = numAdditionalTaxa, y = mean_val),
             color = "black", size = 2.5, shape = 18) +
  geom_text(data = mean_dat,
            aes(x = numAdditionalTaxa, y = mean_val,
                label = paste0("Mean: ", round(mean_val, 2))),
            hjust = -0.1, color = "black", size = 3, fontface = "bold",
            family = "Georgia") +
  geom_text(data = mean_dat,
            aes(x = numAdditionalTaxa, y = lower_95,
                label = paste0("2.5%: ", round(lower_95, 2))),
            hjust = -0.1, color = "black", size = 2.8,
            family = "Georgia") +
  geom_text(data = mean_dat,
            aes(x = numAdditionalTaxa, y = upper_95,
                label = paste0("97.5%: ", round(upper_95, 2))),
            hjust = -0.1, color = "black", size = 2.8,
            family = "Georgia") +
  geom_text(data = null_dat,
            aes(x = numAdditionalTaxa, y = null_val,
                label = paste0("Ex: ", round(null_val, 2))),
            hjust = -0.1, color = "black", size = 2.8,
            family = "Georgia") +
  xlab("Number of additional taxa") +
  ylab("Percent recovered disparity") +
  theme_minimal(base_family = "Georgia") +
  scale_x_discrete(expand = expansion(add = c(0, 1.0)))
p1

ggsave(
  paste0(output, "3Dalpha0.1.svg"),
  p1,
  height = 3.5,
  width = 11
)


#twoD and threeD
mean_dat2d <- twoDimensionalPlotDat %>%
  group_by(numAdditionalTaxa) %>%
  summarise(
    mean_val = mean(percentRecovered, na.rm = T)
  )
mean_dat3d <- threeDimensionalPlotDat %>%
  group_by(numAdditionalTaxa) %>%
  summarise(
    mean_val = mean(percentRecovered, na.rm = T)
  )

p1 <- ggplot() +
  geom_half_violin(data = threeDimensionalPlotDat,
                   aes(x = numAdditionalTaxa, y = percentRecovered), 
                   fill = "#215834", 
                   side = "r") +
  geom_half_violin(data = twoDimensionalPlotDat,
                   aes(x = numAdditionalTaxa, y = percentRecovered), 
                   fill = "#3b7ea1", 
                   side = "l") +
  geom_point(data = mean_dat3d,
             aes(x = numAdditionalTaxa, y = mean_val),
             color = "black", size = 2.5, shape = 18) +
  geom_text(data = mean_dat3d,
            aes(x = numAdditionalTaxa, y = mean_val,
                label = paste0("3D: ", round(mean_val, 2))),
            hjust = -0.1, color = "black", size = 3, fontface = "bold",
            family = "Georgia") +
  geom_point(data = mean_dat2d,
             aes(x = numAdditionalTaxa, y = mean_val),
             color = "black", size = 2.5, shape = 18) +
  geom_text(data = mean_dat2d,
            aes(x = numAdditionalTaxa, y = mean_val,
                label = paste0("2D: ", round(mean_val, 2))),
            hjust = 1.1, color = "black", size = 3, fontface = "bold",
            family = "Georgia") +
  xlab("Number of additional taxa") +
  ylab("Percent recovered disparity") +
  theme_minimal(base_family = "Georgia") +
  scale_x_discrete(expand = expansion(add = c(0.75, 0.75)))
p1

ggsave(
  paste0(output, "twoDthreeDCompared.svg"),
  p1,
  height = 3.5,
  width = 11
)

for(i in unique(twoDimensionalPlotDat$numAdditionalTaxa)){
  print(wilcox.test(filter(twoDimensionalPlotDat, twoDimensionalPlotDat$numAdditionalTaxa ==i)$percentRecovered, filter(threeDimensionalPlotDat, threeDimensionalPlotDat$numAdditionalTaxa ==i)$percentRecovered))
}
