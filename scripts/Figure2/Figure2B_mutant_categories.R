# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

########################################
## Figure 2B - Mutant category splits ##
########################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "ggplot2")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Pre-processed DiMSum data ##
##################################

## Load
load("../../data/DiMSum/PYL1-ABI1/PYL1-ABI1_preprocessed.RData")

## remove synonymous variants
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){y <- x[-which(x[,"Nham_aa"] == 0 & is.na(x[,"WT"]) == T),]; return(y)})
PYL1.ABI1 <- lapply(PYL1.ABI1, function(x){y <- x[-which(PYL1.ABI1[[1]][,"Nham_aa"] == 0 & is.na(PYL1.ABI1[[1]][,"WT"]) == F)[-1],]; return(y)})

## summarise in one matrix
PYL1.summary <- matrix(NA, nrow = 177*21, ncol = 12)
colnames(PYL1.summary) <- names(PYL1.ABI1)
ids <- c()
for(i in 1:177){
  tmp.pos <- rep(str_split_fixed(PYL1.ABI1$`2500`[which(PYL1.ABI1$`2500`$WT)[1],"aa_seq"], "", 177)[i], 21)
  ids <- c(ids, paste0(tmp.pos, i + 32))
}
rownames(PYL1.summary) <- paste0(ids, rep(c("G", "A", "V", "L", "M",
                                            "I", "F", "Y", "W", "K",
                                            "R", "H", "D", "E", "S",
                                            "T", "C", "N", "Q", "P",
                                            "*"),177))
WT <- str_split_fixed(PYL1.ABI1$`2500`[which(PYL1.ABI1$`2500`$WT)[1],"aa_seq"],"",177)[1,]
for(j in 1:length(PYL1.ABI1)){
  
  tmp.vars <- str_split_fixed(PYL1.ABI1[[j]]$aa_seq,"",177)
  
  for(i in 1:nrow(tmp.vars)){
    
    ### only take the "best" WT
    if(all(tmp.vars[i,] == WT)){
      
      PYL1.summary[paste0(WT, 33:209, WT),j] <- PYL1.ABI1[[j]][i,"gr_normalised_WTscaled"]
      next
      
    }
    
    pos.mut <- which(tmp.vars[i,] != WT)
    
    ### skip if higher-order mutant
    if(length(pos.mut) > 1){
      
      next
      
    }
    
    ### otherwise obtain the corresponding fitness
    tmp.mut <- paste0(WT[pos.mut], pos.mut + 32, tmp.vars[i,pos.mut])
    PYL1.summary[tmp.mut,j] <- PYL1.ABI1[[j]][i,"gr_normalised_WTscaled"]
    
  }
  
}

## remove WT repetition 
WTs.pos <- match(paste0(WT, 33:209, WT), rownames(PYL1.summary))
PYL1.summary.nonWT <- PYL1.summary[-WTs.pos[-1],]
PYL1.summary.nonWT <- PYL1.summary.nonWT[c(16,1:15,17:nrow(PYL1.summary.nonWT)),]
rownames(PYL1.summary.nonWT)[1] <- "WT"

PYL1.summary.nonWT <- PYL1.summary.nonWT[-sort(unique(do.call(c,apply(PYL1.summary.nonWT[,1:12], 2, function(x){which(is.na(x) == T)})))),1:12]
colnames(PYL1.summary.nonWT) <- round(as.numeric(colnames(PYL1.summary.nonWT)),2)


## 2. Doughnut distribution of mutants ##
#########################################

## distribution, plainly inferred from unbiased clustering
mut.clust <- cutree(hclust(dist(PYL1.summary.nonWT)), k = 6)
mut.clust[which(mut.clust == 1 | mut.clust == 3)] <- "Hill shape"
mut.clust[which(mut.clust == 4 | mut.clust == 6)] <- "Constitutive binding"
mut.clust[which(mut.clust == 2)] <- "Minimal binding"
mut.clust[which(mut.clust == 5)] <- "Inverted Hill shape"

## after manual inspection, overwrite four inverted Hill shape muts
mut.clust["S119Q"] <- "Inverted Hill shape"
mut.clust["S119R"] <- "Inverted Hill shape"
mut.clust["S119Y"] <- "Inverted Hill shape"
mut.clust["H142P"] <- "Inverted Hill shape"

## rescue those ~56 incomplete measurements which could not be h-clustered?

# Load the required package

# Prepare the data frame
doughnut.in <- data.frame(category = paste0(names(sort(table(mut.clust), decreasing = T)), " (", round(100*sort(table(mut.clust), decreasing = T)/sum(sort(table(mut.clust), decreasing = T)), 1), "%)"),
                          value = as.numeric(sort(table(mut.clust), decreasing = T)))
doughnut.in$category <- factor(doughnut.in$category, 
                               levels = rev(doughnut.in$category))

# Create a basic doughnut plot
custom_colors <- c("grey80", "grey60", "grey40", "grey20")

pdf('../../results/Figure2/Figure2B_mutant_distribution_doughnut.pdf', width = 10, height = 10)

ggplot(doughnut.in, aes(x = 2, y = value, fill = category)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_polar(theta = "y") +
  xlim(0.1, 2.9) +
  scale_fill_manual(values = custom_colors) +
  theme_void() +
  theme(legend.position = "none")

dev.off()


## 3. Version ##
################

# sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sonoma 14.6.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Europe/Madrid
# tzcode source: internal
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ggplot2_3.5.1 scales_1.3.0  stringr_1.5.1
# 
# loaded via a namespace (and not attached):
# [1] labeling_0.4.3    R6_2.6.1          tidyselect_1.2.1  farver_2.1.2      magrittr_2.0.3    gtable_0.3.6     
# [7] glue_1.8.0        tibble_3.2.1      pkgconfig_2.0.3   generics_0.1.3    dplyr_1.1.4       lifecycle_1.0.4  
# [13] cli_3.6.4         grid_4.4.1        vctrs_0.6.5       withr_3.0.2       compiler_4.4.1    rstudioapi_0.17.1
# [19] tools_4.4.1       pillar_1.10.1     munsell_0.5.1     colorspace_2.1-1  rlang_1.1.5       stringi_1.8.4  