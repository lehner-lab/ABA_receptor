# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###################################
## Figure 3A - EC50 Volcano plot ##
###################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "scales", "reshape", "ggplot2", "ggtext", "ggrepel")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Identify sites with significantly lower EC50 ##
#####################################################

## input
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")

## remove stops
Hill.parameters.filtered <- parameters.Hill[-grep("[*]", rownames(parameters.Hill)),]

## calculate log10-fold change between WT and variant
ec50_logFC <- log10(Hill.parameters.filtered[-1,"EC50"] / Hill.parameters.filtered["WT","EC50"])

## fetch EC50s and values
all_ec50 <- log10(Hill.parameters.filtered[-1,"EC50"])
all_ec50_se <- Hill.parameters.filtered[-1,"EC50 SE"] / c(Hill.parameters.filtered[-1,"EC50"] * log(10))
ref_ec50 <- log10(Hill.parameters.filtered["WT","EC50"])
ref_se <- Hill.parameters.filtered["WT","EC50 SE"] / c(Hill.parameters.filtered["WT","EC50"] * log(10))

## calculate Z-scores
z_scores <- (all_ec50 - ref_ec50) / sqrt(all_ec50_se^2 + ref_se^2)

## calculate p-values (BH adjusted)
p.ec50 <- 2 * (1 - pnorm(abs(z_scores)))
p.adj.ec50 <- p.adjust(p.ec50, method = "BH")

## summarise data
data.out <- as.data.frame(cbind("logFC" = ec50_logFC,
                                "p" = p.ec50,
                                "p.adj" = p.adj.ec50))
data.out$label <- names(ec50_logFC)
data.out$p.adj[data.out$p.adj == 0] <- 1e-18


## 2. Display ##
################

## highlight selected hits
labs.out <- subset(data.out, `p.adj` < 0.1 & `logFC` < 0)
labs.out <- labs.out[order(labs.out$logFC),]
labs.out <- data.out[grep("H87P|H87V|H87R|H87Q|H87C|H87A
                           |D107V|V108I|N109F|V110R|N117K|T118W|S119W|L125I
                           |T138I|T138V|E141D|K158P|R164D|R164N|I165E|W166M|V168M|E171T|S182T
                           |D185E|L188T|A190V|V193I|L196T|Q199R|K200S|S203K", rownames(data.out)),]

pdf("../../results/Figure3/Figure3A_Volcano_EC50.pdf", width = 7, height = 8)

out.3A <- ggplot(data.out, aes(x = logFC, y = -log10(p.adj))) +
  geom_point(data = subset(data.out, !(p.adj < 0.1 & logFC < 0)),
             aes(x = logFC, y = -log10(p.adj), color = squish(logFC, range = c(-4.5, 4.5))),
             shape = 16,
             size = 1,
             alpha = 0.3) +
  geom_point(data = subset(data.out, p.adj < 0.1 & logFC < 0),
             aes(x = logFC, y = -log10(p.adj), fill = squish(logFC, range = c(-4.5, 4.5))),
             shape = 21,
             color = "black",
             size = 3,
             alpha = 1) +
  scale_colour_gradient2(low = "blue", mid = "gray90", high = "red", 
                         midpoint = 0, 
                         limits = c(-2,2),
                         oob = squish) +
  scale_fill_gradient2(low = "blue", mid = "gray90", high = "red", 
                       midpoint = 0, 
                       limits = c(-2,2),
                       oob = squish) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black", alpha = 0.2) +
  geom_text_repel(data = labs.out,
                  aes(label = label),
                  size = 4, 
                  family = "Helvetica",
                  max.overlaps = 1000,
                  force = 2,
                  force_pull = 0.5,
                  box.padding = 0.85,
                  point.padding = 0,
                  segment.size = 0.3,
                  segment.color = "darkgrey",
                  color = "black") +
  scale_x_continuous(breaks = seq(f = -5, t = 5, by = 1), 
                     limits = c(-2, 1),
                     labels = c("A", "B", "C", "0.01", "0.1", "1", "10", "100", "I", "J", "K")) +
  scale_y_continuous(breaks = seq(f = 0, t = 18, by = 3), 
                     limits = c(0, 18), oob = scales::squish) +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        legend.position = "none",
        title = element_text(size = 40),
        axis.text = element_text(size = 18),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 25, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 25, vjust = 1),
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(1, 1, 1, 1),"cm")) +
  labs(x = bquote(EC[50] * " fold change"),
       y = bquote(-log[10] * "(adjusted P)"))

print(out.3A)

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
# [1] ggrepel_0.9.6 ggtext_0.1.2  ggplot2_3.5.1 reshape_0.8.9 scales_1.3.0  stringr_1.5.1
# 
# loaded via a namespace (and not attached):
# [1] crayon_1.5.3      vctrs_0.6.5       cli_3.6.4         rlang_1.1.5       stringi_1.8.4     generics_0.1.3   
# [7] labeling_0.4.3    glue_1.8.0        colorspace_2.1-1  plyr_1.8.9        gridtext_0.1.5    grid_4.4.1       
# [13] munsell_0.5.1     tibble_3.2.1      lifecycle_1.0.4   compiler_4.4.1    dplyr_1.1.4       Rcpp_1.0.14      
# [19] pkgconfig_2.0.3   rstudioapi_0.17.1 farver_2.1.2      R6_2.6.1          tidyselect_1.2.1  pillar_1.10.1    
# [25] magrittr_2.0.3    withr_3.0.2       tools_4.4.1       gtable_0.3.6      xml2_1.3.6 