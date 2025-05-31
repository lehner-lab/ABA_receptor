# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

##############################################
## Figure 2A - residue conservation display ##
##############################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "scales")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input ConSurf-DB output ##
################################

## import
consurf <- read.table("../../data/PYL1_consurf_grades.txt", sep = "\t", 
                      col.names = c("POS", "SEQ", "ATOM",
                                    "SCORE", "COLOR",	"CONFIDENCE INTERVAL",
                                    "B/E",	"F/S",	"MSA DATA", "RESIDUE VARIETY"),
                      skip = 28, strip.white = T)

consurf <- consurf[-c(1:13,191:202),]
consurf[,"POS"] <- 33:209


## 2. Plot ##
#############

pdf("../../results/Figure2/Figure2A_conservation_layer.pdf", width = 25, height = 2.5)

plot(x = 33:209,
     y = -consurf$SCORE, 
     pch = 16, 
     type = "s",
     lwd = 2,
     bty = "n",
     xaxt = "n",
     xlab = "",
     yaxt = "n",
     ylab = "",
     xlim = c(33,210))

lines(x = c(209, 210), y = c(-consurf$SCORE[177], -consurf$SCORE[177]), lwd = 2)

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
# [1] scales_1.3.0  stringr_1.5.1
# 
# loaded via a namespace (and not attached):
# [1] colorspace_2.1-1  compiler_4.4.1    R6_2.6.1          magrittr_2.0.3    cli_3.6.4         tools_4.4.1      
# [7] glue_1.8.0        rstudioapi_0.17.1 stringi_1.8.4     lifecycle_1.0.4   munsell_0.5.1     rlang_1.1.5    