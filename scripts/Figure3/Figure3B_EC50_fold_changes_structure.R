# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###############################################################
## Figure 3B - min. EC50 fold changes on PYL1-ABI1 structure ##
###############################################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "bio3d")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Collect EC50 fold-changes within residues ##
##################################################

## input
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")

## filter: remove stops
Hill.parameters.filtered <- parameters.Hill[-grep("[*]", rownames(parameters.Hill)),]

## residue-level: summarise EC50s and EC50 fold-changes
ec50_per_res <- matrix(NA, nrow = 177, ncol = 21)
rownames(ec50_per_res) <- paste0(strsplit("TQDEFTQLSQSIAEFHTYQLGNGRCSSLLAQRIHAPPETVWSVVRRFDRPQIYKHFIKSCNVSEDFEMRVGCTRDVNVISGLPANTSRERLDLLDDDRRVTGFSITGGEHRLRNYKSVTTVHRFEKEEEEERIWTVVLESYVVDVPEGNSEEDTRLFADTVIRLNLQKLASITEAMN", "")[[1]],
                                 33:209)
colnames(ec50_per_res) <- c(c("G", "A", "V", "L", "M", "I", "F", 
                              "Y", "W", "K", "R", "H", "D", "E", 
                              "S", "T", "C", "N", "Q", "P","*"))
for(i in 2:nrow(parameters.Hill)){
  tmp.id <- rownames(parameters.Hill)[i]
  aa.id <- substr(tmp.id, 1, nchar(tmp.id) - 1)
  mut.id <- substr(tmp.id, nchar(tmp.id), nchar(tmp.id))
  ec50_per_res[aa.id,mut.id] <- parameters.Hill[i,"EC50"]
}
ec50_per_res_FC <- ec50_per_res / parameters.Hill["WT","EC50"]
ec50_per_res_FC_log10 <- log10(ec50_per_res_FC)
ec50_per_res_FC_log10.min <- apply(ec50_per_res_FC_log10[,-21], 1, min, na.rm = T) + 4.5


## 2. Project on complex structure ##
#####################################

PYL1_ABI1.pdb <- read.pdb(file = "../../data/PDB/original/Yin_2009_ABA-PYL1-ABI1/pdb3kdj.ent")
PYL1_ABI1_min_EC50_fc <- PYL1_ABI1.pdb

### overwrite B-factor values of PYL1
for (i in 1:length(ec50_per_res_FC_log10.min)){
  
  PYL1_ABI1_min_EC50_fc$atom[which(PYL1_ABI1_min_EC50_fc$atom[,"resno"] == i + 32),"b"] <- ec50_per_res_FC_log10.min[i]
  
}
PYL1_ABI1_min_EC50_fc$atom[which(PYL1_ABI1_min_EC50_fc$atom$chain == "A" & c(PYL1_ABI1_min_EC50_fc$atom[,"resno"] == 31 | PYL1_ABI1_min_EC50_fc$atom[,"resno"] == 32)),"b"] <- 100

### overwrite the B-factor values of the corresponding ABI1 (chain B)
PYL1_ABI1_min_EC50_fc$atom[which(PYL1_ABI1_min_EC50_fc$atom$chain == "B"),"b"] <- 100

### export
write.pdb(PYL1_ABI1_min_EC50_fc, file = paste0("../../data/PDB/modified/Figure3B_min_EC50_FC_new.ent"))


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
# [1] bio3d_2.4-5   scales_1.3.0  stringr_1.5.1
# 
# loaded via a namespace (and not attached):
# [1] R6_2.6.1          magrittr_2.0.3    glue_1.8.0        parallel_4.4.1    lifecycle_1.0.4   cli_3.6.4        
# [7] grid_4.4.1        compiler_4.4.1    rstudioapi_0.17.1 tools_4.4.1       munsell_0.5.1     Rcpp_1.0.14      
# [13] colorspace_2.1-1  rlang_1.1.5       stringi_1.8.4  