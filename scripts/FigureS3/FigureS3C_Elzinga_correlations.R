# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

##################################################################################
## Supplementary Figure 3C - correlation of Elzinga PYR1 mutants vs. PYL1 EC50s ##
##################################################################################


## 0. Environment ##
####################

# Libraries
packages <- c("stringr", "readxl", "growthrates", "drc", 
              "scales", "reshape", "ggplot2", "ggtext")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Identify sites with significantly lower EC50 ##
#####################################################

## input
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")


## 2. Load Elzinga et al. single mutations (ACS ChemBio 2019) ##
################################################################

## note that the exact system here is:
## PYR1 (not PYL1)
## HAB1 (not ABI1)
## Y2M (not DHFR complementation)
pyr1.muts.elzinga <- matrix(NA, nrow = 9, ncol = 2)
rownames(pyr1.muts.elzinga) <- c("WT", "A160C", "A160V", "F61L",
                                 "F61M", "I110C", "I110S", "V81I",
                                 "E141L")
colnames(pyr1.muts.elzinga) <- c("EC50", "SD")
pyr1.muts.elzinga[,"EC50"] <- c(0.864, 0.167, 0.0982, 0.613,
                                2.98, 0.873, 1.15, 0.242,
                                0.238)
pyr1.muts.elzinga[,"SD"] <- c(0.061, 0.022, 0.011, 0.079,
                              1.1, 0.17, 0.28, 0.242,
                              0.238)

## translate the mutants to PYL1 positions
rownames(pyr1.muts.elzinga)[-1] <- paste0(substr(rownames(pyr1.muts.elzinga)[-1], 1, 1),
                                          as.numeric(substr(rownames(pyr1.muts.elzinga)[-1], 2, nchar(rownames(pyr1.muts.elzinga)[-1]) - 1)) + 30,
                                          substr(rownames(pyr1.muts.elzinga)[-1], nchar(rownames(pyr1.muts.elzinga)[-1]), nchar(rownames(pyr1.muts.elzinga)[-1])))
rownames(pyr1.muts.elzinga)[4] <- "F88L"
rownames(pyr1.muts.elzinga)[5] <- "F88M"
rownames(pyr1.muts.elzinga)[6] <- "I137C"
rownames(pyr1.muts.elzinga)[7] <- "I137S"
rownames(pyr1.muts.elzinga)[8] <- "V108I"


## 3. Combine and plot ##
#########################

data.out <- pyr1.muts.elzinga
data.out <- cbind(data.out,
                  parameters.Hill[match(rownames(data.out), rownames(parameters.Hill)),c(4,8)])
colnames(data.out) <- c("Elzinga EC50", "Elzinga EC50_SD", "DMS EC50", "DMS EC50_SE")
data.out <- as.data.frame(data.out)

## correlation
r.EC50 <- cor(x = log(data.out$`Elzinga EC50`),  y = log(data.out$`DMS EC50`), method = "pearson") ## 0.62
p.EC50 <- summary(lm(log10(data.out$`Elzinga EC50`) ~ log10(data.out$`DMS EC50`)))$coefficients[2,4] ## P = 0.07292539


pdf('../../results/FigureS3/FigureS3C_Elzinga_correlation.pdf', height = 15, width = 15)

ggplot(data.out, aes(x = `Elzinga EC50`, y = `DMS EC50`)) +
  scale_x_log10(breaks = c(0.1, 0.1, 1, 10), 
                labels = c(0.1, 0.1, 1, 10),
                limits = c(0.0000000001,1000000)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10), 
                labels = c(0.01, 0.1, 1, 10),
                limits = c(0.0000000001,1000000)) +
  coord_cartesian(xlim = c(0.1, 75), ylim = c(0.1, 75), expand = T) +
  geom_smooth(data = data.out,
              mapping = aes(x = `Elzinga EC50`, 
                            y = `DMS EC50`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = NA) +
  geom_point(aes(x = `Elzinga EC50`, y = `DMS EC50`), pch = 16, 
             col = "black", size = 9) +
  geom_errorbar(aes(ymin = `DMS EC50` - `DMS EC50_SE`, ymax = `DMS EC50` + `DMS EC50_SE`),
                color = "black", linewidth = 1) +
  geom_errorbarh(aes(xmin = `Elzinga EC50` - `Elzinga EC50_SD`, xmax = `Elzinga EC50` + `Elzinga EC50_SD`),
                 color = "black", linewidth = 1) +
  annotate("text",
           x = 0.1,
           y = 60,
           label = bquote(italic(r) == .(format(r.EC50, digits = 2, nsmall = 2)) ~ 
                            ", " ~ italic(P) == .(format(p.EC50, digits = 2, nsmall = 2))),
           hjust = 0, size = 15, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 40),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 2),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = expression(atop("PYR1-HAB1 "*EC[50]* " (µM)", "(Elzinga et al., ACS Chem. Biol. 2019)")),
       y = expression(atop("PYL1-ABI1 "*EC[50]* " (µM)", "(this study)")))

dev.off()


## 4. Version ##
################

