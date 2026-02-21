library(ggplot2)
library(ggtree)
library(phytools)

output <- "/Users/levir/Documents/GitHub/HomininTaxicDiversity/figures/"

### LDDMM Stacks ###
treeList <- read.delim(
  gzfile("/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/simulatedTrees/numLandmarks10/numAdditionalTaxa5/alpha0.10/TreeIndex0_TreeAlpha1.00_TreeBeta1.00_NumAdditionalTaxa5.tsv.gz"),
  header = F
)
tree <- read.tree(text =treeList[1,1])

tree$tip.label <- gsub("_", tree$tip.label, replacement = " ")
p1 <- ggtree(tree)+
  geom_tiplab(fontface = 4)+
  xlim(NA, 5)
p1 <- ggtree::rotate(p1, 32)
p1
ggsave(paste(output, "tree1LDDMM.svg", sep = ""), p1)

#LDDMMM stacks
shapeInputFile2D <- 
  paste0(
    "/Users/levir/Documents/GitHub/HomininTaxicDiversity/results/simulatedShapes/",
    "numLandmarks", 10, "/",
    "numAdditionalTaxa", 3 ,  "/" ,
    "alpha" , 0.1 , "0" ,  "/" ,
    "2D_TreeIndex" , 0,
    "_TreeAlpha1.00_TreeBeta1.00_" , 
    "NumAdditionalTaxa", 3 ,
    "_numLandmarks" , 10 ,
    "_lddmmAlpha" , 0.1 , "0" ,
    "_Sigma1.00",
    "_rep" , 0 ,
    "nodeShapes.tsv.gz"
  )
lddmmRes1Alpha0.2 <- read.delim(gzfile(shapeInputFile2D), header = FALSE)
colnames(lddmmRes1Alpha0.2) <- c("taxon", "id", "x", "y")
lddmmRes1Alpha0.2$taxon <- gsub("_", lddmmRes1Alpha0.2$taxon, replacement = " ")

tip_order <- p1$data %>%
  filter(isTip) %>%
  arrange(-y) %>%
  pull(label)

# Generate plots with fixed axis limits
x_range <- range(lddmmRes1Alpha0.2$x)

# Make individual plots with custom y range
stackFunc <- function(tip) {
  dat <- filter(lddmmRes1Alpha0.2, taxon == tip)
  
  #center dat
  dat$x <- dat$x - mean(dat$x)
  dat$y <- dat$y - mean(dat$y)
  
  y_pad <- 0.5 # Add vertical padding
  y_min <- min(dat$y) - y_pad
  y_max <- max(dat$y) + y_pad
  
  x_pad <- 0.5  # Add vertical padding
  x_min <- min(dat$x) - x_pad
  x_max <- max(dat$x) + x_pad
  
  connections_sim <- dat %>%
    arrange(id) %>%
    mutate(lm_next = lead(id),
           x_next = lead(x),
           y_next = lead(y))
  
  # For the last point, set next to the first point (to close the loop)
  connections_sim[nrow(connections_sim), c("id_next", "x_next", "y_next")] <- 
    c(dat$id[1], dat$x[1], dat$y[1])
  
  ggplot(dat, aes(x = x, y = y)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_segment(data = connections_sim,
                 aes(x = x, y = y, xend = x_next, yend = y_next),
                 color = "black", size = 1)+
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    theme_minimal() +
    ggtitle(tip) +
    theme(
      plot.title = element_text(hjust = 0, size = 9),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
}
plots <- lapply(tip_order, stackFunc)

# Stack vertically
stack_plot <-patchwork::wrap_plots(plots, ncol = 1)
stack_plot
ggsave(paste(output, "fig1LDDMMAlpha0.1.svg", sep = ""), 
       stack_plot, 
       width = 2,
       height = 42)


lddmmRes1Alpha0.2 <- read.delim("resultsGit/twoDimensionSimTreeIndex0LM10Alpha0.200000Dataset1nodeShapes.tsv", header = FALSE)
colnames(lddmmRes1Alpha0.2) <- c("taxon", "id", "x", "y")
lddmmRes1Alpha0.2$taxon <- gsub("_", lddmmRes1Alpha0.2$taxon, replacement = " ")

tip_order <- p1$data %>%
  filter(isTip) %>%
  arrange(-y) %>%
  pull(label)

# Generate plots with fixed axis limits
x_range <- range(lddmmRes1Alpha0.2$x)

# Make individual plots with custom y range
plots <- lapply(tip_order, stackFunc)

# Stack vertically
stack_plot <-patchwork::wrap_plots(plots, ncol = 1)
stack_plot
ggsave(paste(output, "fig1LDDMMAlpha0.2.svg", sep = ""), 
       stack_plot, 
       width = 2,
       height = 42)

lddmmRes1Alpha0.2 <- read.delim("resultsGit/twoDimensionSimTreeIndex0LM10Alpha0.300000Dataset1nodeShapes.tsv", header = FALSE)
colnames(lddmmRes1Alpha0.2) <- c("taxon", "id", "x", "y")
lddmmRes1Alpha0.2$taxon <- gsub("_", lddmmRes1Alpha0.2$taxon, replacement = " ")

tip_order <- p1$data %>%
  filter(isTip) %>%
  arrange(-y) %>%
  pull(label)

# Generate plots with fixed axis limits
x_range <- range(lddmmRes1Alpha0.2$x)

# Make individual plots with custom y range
plots <- lapply(tip_order, stackFunc)

# Stack vertically
stack_plot <-patchwork::wrap_plots(plots, ncol = 1)
stack_plot
ggsave(paste(output, "fig1LDDMMAlpha0.3.svg", sep = ""), 
       stack_plot, 
       width = 2,
       height = 42)

lddmmRes1Alpha0.2 <- read.delim("resultsGit/twoDimensionSimTreeIndex0LM10Alpha0.400000Dataset1nodeShapes.tsv", header = FALSE)
colnames(lddmmRes1Alpha0.2) <- c("taxon", "id", "x", "y")
lddmmRes1Alpha0.2$taxon <- gsub("_", lddmmRes1Alpha0.2$taxon, replacement = " ")

tip_order <- p1$data %>%
  filter(isTip) %>%
  arrange(-y) %>%
  pull(label)

# Generate plots with fixed axis limits
x_range <- range(lddmmRes1Alpha0.2$x)

# Make individual plots with custom y range
plots <- lapply(tip_order, stackFunc)

# Stack vertically
stack_plot <-patchwork::wrap_plots(plots, ncol = 1)
stack_plot
ggsave(paste(output, "fig1LDDMMAlpha0.4.svg", sep = ""), 
       stack_plot, 
       width = 2,
       height = 42)