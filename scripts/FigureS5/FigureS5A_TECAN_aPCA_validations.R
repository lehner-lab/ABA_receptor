# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###########################################################
## Supplementary Figure 5A - PYL1-PYL1 TECAN validations ##
###########################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "readxl", "growthrates", "scales", 
              "reshape", "ggplot2", "ggtext")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input PYL1-PYL1 data ##
#############################

## import

load("../../data/DiMSum/PYL1-PYL1/0uM_ABA/PYL1-PYL1_0uM_ABA_preprocessed.RData")
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.0uM.ABA[,"WT"]) == T),]
PYL1.PYL1.0uM.ABA <- PYL1.PYL1.0uM.ABA[-which(PYL1.PYL1.0uM.ABA[,"Nham_aa"] == 0 & is.na(PYL1.PYL1.0uM.ABA[,"WT"]) == F)[-1],]


## 2. Import TECAN data ##
##########################

## raw data import
PYL1.PYL1.TECAN <- as.matrix(read_xlsx('../../data/TECAN/PYL1_validations/pMS-PYL1-PYL1_mutant_validations_TECAN.xlsx', sheet = 1))[22:407,]
rownames(PYL1.PYL1.TECAN) <- PYL1.PYL1.TECAN[,1]
colnames(PYL1.PYL1.TECAN) <- PYL1.PYL1.TECAN[1,]
PYL1.PYL1.TECAN <- PYL1.PYL1.TECAN[,-1]
PYL1.PYL1.TECAN <- PYL1.PYL1.TECAN[-c(1:2),]
colnames(PYL1.PYL1.TECAN) <- as.numeric(colnames(PYL1.PYL1.TECAN))/3600 ## convert to hours
class(PYL1.PYL1.TECAN) <- "numeric"
PYL1.PYL1.TECAN <- as.data.frame(PYL1.PYL1.TECAN)

## assign mutants
setup <- as.matrix(read_xlsx('../../data/TECAN/PYL1_validations/pMS-PYL1-PYL1_mutant_validations_TECAN.xlsx', sheet = 1))[2:18,1:25]
colnames(setup) <- setup[1,]
setup <- setup[-1,]
rownames(setup) <- setup[,1]
setup <- setup[,-1]
mut.names <- unique(c(setup["A",],setup["I",]))
mut.names <- mut.names[!is.na(mut.names)]
PYL1.PYL1.TECAN.muts <- vector(mode = "list", length = length(mut.names))
names(PYL1.PYL1.TECAN.muts) <- mut.names
for(i in 1:length(PYL1.PYL1.TECAN.muts)){
  PYL1.PYL1.TECAN.muts[[i]] <- which(setup == mut.names[i], arr.ind = TRUE)
  PYL1.PYL1.TECAN.muts[[i]] <- paste0(rownames(PYL1.PYL1.TECAN.muts[[i]]), PYL1.PYL1.TECAN.muts[[i]][,2])
  PYL1.PYL1.TECAN.muts[[i]] <- PYL1.PYL1.TECAN[PYL1.PYL1.TECAN.muts[[i]],]
}

## calculate growth rates
PYL1.PYL1.TECAN.muts.rates <- PYL1.PYL1.TECAN.muts
for (i in 1:length(PYL1.PYL1.TECAN.muts.rates)){
  
  print(i)
  
  tmp.out <- vector(mode = 'list', length = nrow(PYL1.PYL1.TECAN.muts.rates[[i]]))
  
  for (k in 1:length(tmp.out)){
    
    tmp.out[[k]] <- rbind(as.numeric(colnames(PYL1.PYL1.TECAN.muts.rates[[i]])),
                          as.numeric(PYL1.PYL1.TECAN.muts.rates[[i]][k,]))
    tmp.out[[k]] <- tmp.out[[k]][,!is.na(tmp.out[[k]][2,])]
    tmp.out[[k]] <- fit_easylinear(time = tmp.out[[k]][1,],
                                   y = tmp.out[[k]][2,],
                                   h = 15)
    tmp.out[[k]] <- tmp.out[[k]]@par[['mumax']]
    
  }
  
  PYL1.PYL1.TECAN.muts.rates[[i]] <- do.call(c, tmp.out)
  
}

## rescale to WT
PYL1.PYL1.TECAN.muts.rates <- lapply(PYL1.PYL1.TECAN.muts.rates, function(x){y <- 100*x/mean(PYL1.PYL1.TECAN.muts.rates$WT); return(y)})

## clean-up
rm(i,k,tmp.out)


## 3. Plot ##
#############

## summarise data
out.all <- matrix(NA, nrow = length(PYL1.PYL1.TECAN.muts.rates) - 1, ncol = 4)
colnames(out.all) <- c("bulk_mean", "bulk_sd", "tecan_mean", "tecan_sd")
rownames(out.all) <- names(PYL1.PYL1.TECAN.muts.rates)[-1]

### add bulk library mean/sd
for(i in 1:nrow(out.all)){
  tmp.name <- rownames(out.all)[i]
  tmp.id <- match(tmp.name, paste0(PYL1.PYL1.0uM.ABA[,"WT_AA"],
                                   PYL1.PYL1.0uM.ABA[,"Pos"],
                                   PYL1.PYL1.0uM.ABA[,"Mut"]))
  out.all[i,"bulk_mean"] <- PYL1.PYL1.0uM.ABA[tmp.id,"gr_normalised_WTscaled"]
  out.all[i,"bulk_sd"] <- PYL1.PYL1.0uM.ABA[tmp.id,"gr_sigma_normalised_WTscaled"]
}

### add TECAN mean/sd
out.all[,"tecan_mean"] <- sapply(PYL1.PYL1.TECAN.muts.rates, mean)[-1]
out.all[,"tecan_sd"] <- sapply(PYL1.PYL1.TECAN.muts.rates, sd)[-1]
out.all <- as.data.frame(out.all)

## to calculate Pearson's coefficients
r.out <- cor(x = out.all[,"bulk_mean"],
         y = out.all[,"tecan_mean"],
         method = "pearson", use = "complete.obs")
p.out <- summary(lm(out.all[,"bulk_mean"] ~ out.all[,"tecan_mean"]))$coefficients[2,4] ## P = 0.07292539

pdf("../../results/FigureS5/FigureS5A_abundancePCA_validation.pdf", height = 15, width = 18)

out.S5A <- ggplot(out.all, aes(x = `tecan_mean`, y = `bulk_mean`)) +
  scale_x_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-50, 150)) +
  scale_y_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-50, 150)) +
  coord_cartesian(xlim = c(-5, 130), ylim = c(-5, 130)) +
  geom_point(data = out.all,
             mapping = aes(x = `tecan_mean`, y = `bulk_mean`),
             color = "black", size = 6) +
  geom_errorbarh(aes(xmin = `tecan_mean` - `tecan_sd`, xmax = `tecan_mean` + `tecan_sd`),
                 color = "black") +
  geom_errorbar(aes(ymin = `bulk_mean` - `bulk_sd`, ymax = `bulk_mean` + `bulk_sd`),
                color = "black") +
  geom_smooth(data = out.all,
              mapping = aes(x = `tecan_mean`, y = `bulk_mean`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = NA) +
  annotate("text",
           x = 0,
           y = 120,
           label = bquote(italic(r) == .(format(r.out, digits = 2, nsmall = 2)) ~ 
                            ", " ~ italic(P) == .(format(p.out, digits = 2, nsmall = 2))),
           hjust = 0, size = 20, color = "black") +
  theme_classic(base_size = 50) +
  theme(plot.title = element_markdown(),
        plot.subtitle = element_markdown(size = 20),
        title = element_text(size = 40),
        axis.text = element_text(size = 40),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = "Abundance (individual mutants)",
       y = "Abundance (library sequencing)")

print(out.S5A)

dev.off()


## 4. Version ##
################

