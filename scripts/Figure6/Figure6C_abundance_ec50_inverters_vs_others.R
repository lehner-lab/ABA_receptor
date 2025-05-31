# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

########################################################
## Figure 6C - abundance/EC50 of inverters vs. others ##
########################################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "reshape", "ggplot2", "ggtext", "ggrepel",
              "wesanderson", "beeswarm", "readxl")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Identify hypersensitive and "inverter" variants from hierarchical clustering ##
#####################################################################################

## input hclust results from Table S3
hclust.out <- as.matrix(read_xlsx("../../data/Supplementary Tables/Table-S3.xlsx"))
colnames(hclust.out) <- hclust.out[2,]
hypersensitive.hclust <- hclust.out[which(hclust.out[,"Variant label by hierarchical clustering"] == "Constitutive binding"
                                        | hclust.out[,"Variant label by hierarchical clustering"] == "Inverted Hill shape"),]

hyper.binders <- hypersensitive.hclust[which(hypersensitive.hclust[,"Variant label by hierarchical clustering"] == "Constitutive binding"), 1]
inverted.binders <- hypersensitive.hclust[which(hypersensitive.hclust[,"Variant label by hierarchical clustering"] == "Inverted Hill shape"), 1]


## 2. Fetch raw measurements ##
###############################

## abundance
abundance.out <- as.matrix(read_xlsx("../../data/Supplementary Tables/Table-S2.xlsx", sheet = 2))
rownames(abundance.out) <- abundance.out[,1]
abundance.out <- abundance.out[,-1]
colnames(abundance.out) <- c(0, 250)
abundance.out <- abundance.out[-c(1:3),]

## EC50
ec50.out <- as.matrix(read_xlsx("../../data/Supplementary Tables/Table-S3.xlsx", sheet = 1))
rownames(ec50.out) <- ec50.out[,1]
ec50.out <- ec50.out[,-1]
colnames(ec50.out) <- ec50.out[2,]
ec50.out <- ec50.out[-c(1:2),]


## 3. Plot aPCA beeswarms ##
############################

abundance.comp <- list("Hypersensitive" = as.numeric(abundance.out[match(hyper.binders, rownames(abundance.out)),"0"]),
                       "Inverters" = as.numeric(abundance.out[match(inverted.binders, rownames(abundance.out)),"0"]),
                       "Others" = as.numeric(abundance.out[-sort(unique(c(which(rownames(abundance.out) %in% hyper.binders),
                                                                          which(rownames(abundance.out) %in% inverted.binders)))),"0"]))

wilcox.test(x = abundance.comp$Hypersensitive, y = abundance.comp$Others) # p < 2.2e-16
wilcox.test(x = abundance.comp$Inverters, y = abundance.comp$Others) # p = 2.927e-05

## set up colouring
beeswarm.cols <- c("G" = wes_palette("Darjeeling2", 6, type = "continuous")[6], "P" = wes_palette("Darjeeling2", 6, type = "continuous")[6],
                   "W" = wes_palette("Darjeeling2", 6, type = "continuous")[5], "Y" = wes_palette("Darjeeling2", 6, type = "continuous")[5], 
                   "F" = wes_palette("Darjeeling2", 6, type = "continuous")[5],
                   "A" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "V" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "L" = wes_palette("Darjeeling2", 6, type = "continuous")[4], 
                   "I" = wes_palette("Darjeeling2", 6, type = "continuous")[4], "M" = wes_palette("Darjeeling2", 6, type = "continuous")[4], 
                   "S" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "T" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "C" = wes_palette("Darjeeling2", 6, type = "continuous")[3],
                   "Q" = wes_palette("Darjeeling2", 6, type = "continuous")[3], "N" = wes_palette("Darjeeling2", 6, type = "continuous")[3],
                   "K" = wes_palette("Darjeeling2", 6, type = "continuous")[2], "R" = wes_palette("Darjeeling2", 6, type = "continuous")[2], "H" = wes_palette("Darjeeling2", 6, type = "continuous")[2],
                   "D" = wes_palette("Darjeeling2", 6, type = "continuous")[1], "E" = wes_palette("Darjeeling2", 6, type = "continuous")[1])

labels.cols <- c(names(abundance.out[match(hyper.binders, rownames(abundance.out)),"0"]),
                 names(abundance.out[match(inverted.binders, rownames(abundance.out)),"0"]),
                 names(abundance.out[-sort(unique(c(which(rownames(abundance.out) %in% hyper.binders),
                                                    which(rownames(abundance.out) %in% inverted.binders)))),"0"]))
labels.cols <- substr(labels.cols, nchar(labels.cols), nchar(labels.cols))
labels.cols <- beeswarm.cols[match(labels.cols, names(beeswarm.cols))]
labels.cols[which(is.na(labels.cols) == T)] <- "black"
names(labels.cols[which(is.na(names(labels.cols) == T))]) <- "*"

pdf("../../results/Figure6/Figure6C_aPCA_inverters.pdf", width = 14, height = 8)
par(mar = c(4,13,1,0))

boxplot(abundance.comp, 
        cex = 1, 
        ylim = c(0, 130), 
        outline = F, 
        horizontal = F,
        medcol = 'white',
        whiskcol = 'white',
        boxcol = 'white',
        outcol = 'white',
        border = "white",
        col = "white",
        frame = F,
        xlab = '',
        xaxt = 'n',
        yaxt = 'n')

## mark WT
abline(h = 100, lwd = 3)

## set up label positions (via beeswarm)
labels.pos1 <- beeswarm(abundance.comp$Hypersensitive, method = "swarm", do.plot = FALSE, cex = 1.5)
labels.pos2 <- beeswarm(abundance.comp$Inverters, method = "compactswarm", do.plot = FALSE, cex = 8)
labels.pos3 <- beeswarm(abundance.comp$Others, method = "swarm", do.plot = FALSE, cex = 0.35)

## add mutations to plot
points(x = labels.pos1$x, 
       y = labels.pos1$y,
       col = "black",
       pch = 16,
       cex = 1.5)

labels.out <- rownames(abundance.out[match(inverted.binders, rownames(abundance.out)),])
labels.out <- substr(labels.out, nchar(labels.out), nchar(labels.out))
text(x = labels.pos2$x + 1, 
     y = labels.pos2$y, 
     labels = rownames(abundance.out[match(inverted.binders, rownames(abundance.out)),]),
     col = beeswarm.cols[match(labels.out, names(beeswarm.cols))],
     cex = 1.75)

points(x = labels.pos3$x + 2, 
       y = labels.pos3$y,
       col = "black", 
       pch = 16,
       cex = 0.35)

## add X-axis
axis(1, at = 1, labels = names(abundance.comp)[1], las = 1, col = 'black', 
     cex.axis = 3, lwd = 3, col.ticks = F, padj = 0.5, tick = F)
axis(1, at = 2, labels = names(abundance.comp)[2], las = 1, col = 'black', 
     cex.axis = 3, lwd = 3, col.ticks = F, padj = 0.5, tick = F)
axis(1, at = 3, labels = names(abundance.comp)[3], las = 1, col = 'black', 
     cex.axis = 3, lwd = 3, col.ticks = F, padj = 0.5, tick = F)

## add Y-axis
axis(2, at = seq(f = 0, t = 120, length.out = 7), 
     labels = seq(f = 0, t = 120, length.out = 7), las = 2,
     cex.axis = 2)

## Y label
mtext(text = "Abundance", side = 2, line = 6, cex = 4)

dev.off()


## 4. Plot EC50 beeswarms ##
############################

ec50.comp <- list("Inverters" = as.numeric(ec50.out[match(inverted.binders, rownames(ec50.out)),"EC_50 [mean]"]),
                  "Others" = as.numeric(ec50.out[-sort(unique(c(which(rownames(ec50.out) %in% hyper.binders),
                                                                which(rownames(ec50.out) %in% inverted.binders)))),"EC_50 [mean]"]))

pdf("../../results/Figure6/Figure6C_EC50_inverters.pdf", width = 9.5, height = 8)
par(mar = c(4,13,1,0))

boxplot(ec50.comp, 
        cex = 1, 
        log = "y", 
        ylim = c(0.01, 5000), 
        outline = F, 
        horizontal = F,
        medcol = 'white',
        whiskcol = 'white',
        boxcol = 'white',
        outcol = 'white',
        border = "white",
        col = "white",
        frame = F,
        xlab = '',
        xaxt = 'n',
        yaxt = 'n')

## mark WT
abline(h = as.numeric(ec50.out["WT", "EC_50 [mean]"]), 
       lwd = 3)

## set up label positions (via beeswarm)
labels.pos1 <- beeswarm(ec50.comp$Inverters, method = "compactswarm", do.plot = FALSE, cex = 5.7)
labels.pos2 <- beeswarm(ec50.comp$Others, method = "swarm", do.plot = FALSE, cex = 0.35)

## add mutations to plot
labels.out <- names(ec50.out[match(inverted.binders, rownames(ec50.out)),"EC_50 [mean]"])
labels.out <- substr(labels.out, nchar(labels.out), nchar(labels.out))
text(x = labels.pos1$x, 
     y = labels.pos1$y, 
     labels = names(ec50.out[match(inverted.binders, rownames(ec50.out)),"EC_50 [mean]"]),
     col = beeswarm.cols[match(labels.out, names(beeswarm.cols))],
     cex = 1.2)

points(x = labels.pos2$x + 1, 
       y = labels.pos2$y,
       col = "black",
       pch = 16,
       cex = 0.35)

## add X-axis
axis(1, at = 1, labels = names(ec50.comp)[1], las = 1, col = 'black', 
     cex.axis = 3, lwd = 3, col.ticks = F, padj = 0.5, tick = F)
axis(1, at = 2, labels = names(ec50.comp)[2], las = 1, col = 'black', 
     cex.axis = 3, lwd = 3, col.ticks = F, padj = 0.5, tick = F)

## add Y-axis
axis(2, at = c(0.01, 0.1, 1, 10, 100, 1000),
     labels = c("0.01", "0.1", "1", "10", "100", "1,000"),
     las = 2,
     cex.axis = 2)

## Y labela
mtext(text = bquote(EC[50] ~ "((+)-ABA, µM)"), side = 2, line = 6, cex = 4)

dev.off()


## 5. Version ##
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
# [1] readxl_1.4.3      beeswarm_0.4.0    wesanderson_0.3.7 ggrepel_0.9.6     ggtext_0.1.2      ggplot2_3.5.1    
# [7] reshape_0.8.9     scales_1.3.0      stringr_1.5.1    
# 
# loaded via a namespace (and not attached):
# [1] vctrs_0.6.5       cli_3.6.4         rlang_1.1.5       stringi_1.8.4     generics_0.1.3    glue_1.8.0       
# [7] colorspace_2.1-1  plyr_1.8.9        gridtext_0.1.5    cellranger_1.1.0  grid_4.4.1        munsell_0.5.1    
# [13] tibble_3.2.1      lifecycle_1.0.4   compiler_4.4.1    dplyr_1.1.4       Rcpp_1.0.14       pkgconfig_2.0.3  
# [19] rstudioapi_0.17.1 R6_2.6.1          tidyselect_1.2.1  pillar_1.10.1     magrittr_2.0.3    withr_3.0.2      
# [25] tools_4.4.1       gtable_0.3.6      xml2_1.3.6    