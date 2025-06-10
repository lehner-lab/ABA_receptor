# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

###########################################################
## Supplementary Figure 3B - PYL1-ABI1 TECAN validations ##
###########################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "readxl", "growthrates", "drc", 
              "scales", "reshape", "ggplot2", "ggtext", "cowplot")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Input PYL1-ABI1 data ##
#############################

## input Hill parameter distributions from dose-response curve fits
load("../../data/DRCs/PYL1-ABI1_parameters_Hill.RData")
parameters.Hill <- parameters.Hill[!is.na(parameters.Hill[,1]),]


## 2. Import TECAN data ##
##########################

## raw data import
PYL1.ABI1.TECAN <- as.matrix(read_xlsx('../../data/TECAN/PYL1_validations/pMS-GluePCA1_PYL1_mutant_validations_dose_response_TECAN.xlsx', sheet = 1))[20:405,1:263]
rownames(PYL1.ABI1.TECAN) <- PYL1.ABI1.TECAN[,1]
colnames(PYL1.ABI1.TECAN) <- PYL1.ABI1.TECAN[1,]
PYL1.ABI1.TECAN <- PYL1.ABI1.TECAN[,-1]
PYL1.ABI1.TECAN <- PYL1.ABI1.TECAN[-c(1:2),]
colnames(PYL1.ABI1.TECAN) <- as.numeric(colnames(PYL1.ABI1.TECAN))/3600 ## convert to hours
class(PYL1.ABI1.TECAN) <- "numeric"
PYL1.ABI1.TECAN <- as.data.frame(PYL1.ABI1.TECAN)

## assign mutants
setup <- as.matrix(read_xlsx('../../data/TECAN/PYL1_validations/pMS-GluePCA1_PYL1_mutant_validations_dose_response_TECAN.xlsx', sheet = 1))[1:16,1:25]
rownames(setup) <- setup[,1]
setup <- setup[,-1]
colnames(setup) <- 1:24
mut.names <- unique(c(setup[,1],setup[,13]))
mut.names <- mut.names[!is.na(mut.names)]
PYL1.ABI1.TECAN.muts <- vector(mode = "list", length = length(mut.names))
names(PYL1.ABI1.TECAN.muts) <- mut.names
for(i in 1:length(PYL1.ABI1.TECAN.muts)){
  PYL1.ABI1.TECAN.muts[[i]] <- which(setup == mut.names[i], arr.ind = TRUE)
  PYL1.ABI1.TECAN.muts[[i]] <- paste0(rownames(PYL1.ABI1.TECAN.muts[[i]]), PYL1.ABI1.TECAN.muts[[i]][,2])
  PYL1.ABI1.TECAN.muts[[i]] <- PYL1.ABI1.TECAN[PYL1.ABI1.TECAN.muts[[i]],]
}

## calculate growth rates
PYL1.ABI1.TECAN.muts.rates <- PYL1.ABI1.TECAN.muts
for (i in 1:length(PYL1.ABI1.TECAN.muts.rates)){
  
  print(i)
  
  tmp.out <- vector(mode = 'list', length = nrow(PYL1.ABI1.TECAN.muts.rates[[i]]))
  
  for (k in 1:length(tmp.out)){
    
    tmp.out[[k]] <- rbind(as.numeric(colnames(PYL1.ABI1.TECAN.muts.rates[[i]])),
                          as.numeric(PYL1.ABI1.TECAN.muts.rates[[i]][k,]))
    tmp.out[[k]] <- tmp.out[[k]][,!is.na(tmp.out[[k]][2,])]
    tmp.out[[k]] <- fit_easylinear(time = tmp.out[[k]][1,],
                                   y = tmp.out[[k]][2,],
                                   h = 15)
    tmp.out[[k]] <- tmp.out[[k]]@par[['mumax']]
    
  }
  
  PYL1.ABI1.TECAN.muts.rates[[i]] <- do.call(c, tmp.out)
  
}

## Dosages
dosages <- rep(NA, 12)
dosages[1] <- 2500
dosages[12] <- 0
for(i in 2:11){
  dosages[i] <- 2500 * (1/3.5)^c(i-1)
}

## fit WT curve(s)
WT.PYL1.drc <- cbind(PYL1.ABI1.TECAN.muts.rates$WT[1:12],rev(dosages))
class(WT.PYL1.drc) <- "numeric"
WT.PYL1.drc <- as.data.frame(WT.PYL1.drc)
colnames(WT.PYL1.drc) <- c("GR", "concentration")

WT.PYL1.drc <- drm(WT.PYL1.drc$GR ~ WT.PYL1.drc$concentration,
                   fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                   type = 'continuous')
WT.PYL1.drc.par <- WT.PYL1.drc$fit$par
names(WT.PYL1.drc.par) <- c("Hill", "B[0]", "B[inf]", "EC50")
WT.PYL1.drc.par <- WT.PYL1.drc.par[c(2:4,1)]
WT.PYL1.drc.par[4] <- -WT.PYL1.drc.par[4]

## rescale all curves to WT Bmax
PYL1.ABI1.TECAN.muts.rates <- lapply(PYL1.ABI1.TECAN.muts.rates, function(x){y <- 100*x/WT.PYL1.drc.par["B[inf]"]; return(y)})

## clean-up
rm(i,k,tmp.out)


## 3. Fit all TECAN curves' dose response profiles ##
#####################################################

parameters.Hill.TECAN <- matrix(NA, nrow = length(PYL1.ABI1.TECAN.muts.rates), ncol = 16)
colnames(parameters.Hill.TECAN) <- c("Hill", "B[0]", "B[inf]", "EC50", 
                                     "Hill SE", "B[0] SE", "B[inf] SE", "EC50 SE", 
                                     "Hill P", "B[0] P", "B[inf] P", "EC50 P", 
                                     "Data points", "AIC", "Residual var", "R^2")
rownames(parameters.Hill.TECAN) <- names(PYL1.ABI1.TECAN.muts.rates)

## Calculate and plot all the dose-response curves with the 4-parametric Hill model
#pdf("../results/FigureS3B_all_curves.pdf", height = 15, width = 18)
#pdf("../results/FigureS3B_inverter_curves.pdf", height = 15, width = 18)
for(i in 1:nrow(parameters.Hill.TECAN)){
  
  print(i)
  
  ### DRC modelling, using the R DRM package
  
  ### output parameter translations:
  ### b: Hill coefficient, i.e. steepness of the curve (n)
  ### c: basal binding fitness (B[0])
  ### d: saturated binding fitness (B[inf])
  ### e: curve inflection point (EC50)
  ### model: binding(ABA conc.) = c + ((d-c)/(1 + (e/x)^b))
  
  ### fetch genotype's data
  tmp.PYL1.drc.in <- cbind(PYL1.ABI1.TECAN.muts.rates[[i]][1:12],rev(dosages))
  class(tmp.PYL1.drc.in) <- "numeric"
  tmp.PYL1.drc.in <- as.data.frame(tmp.PYL1.drc.in)
  colnames(tmp.PYL1.drc.in) <- c("GR", "concentration")
    
  #if(!rownames(parameters.Hill.TECAN)[i] %in% c("V110H", "A116H", "V193H", "V193W", "A190E", "H87P")){
  if(!rownames(parameters.Hill.TECAN)[i] %in% c("V110H", "A116H", "V193H", "V193W", "A116R", "A190E", "L125I")){
  #if(rownames(parameters.Hill.TECAN)[i] %in% c("V110Y", "V110W", "A116W", "H142W")){
    
    ### curve fit
    tmp.PYL1.drc <- drm(tmp.PYL1.drc.in$GR ~ tmp.PYL1.drc.in$concentration,
                        fct = LL.4(fixed = c(NA, NA, NA, NA), names = c("Hill", "B[0]", "B[inf]", "EC50")),
                        type = 'continuous')
    
    ### parameter take-over
    parameters.Hill.TECAN[i,"Data points"] <- sum(is.na(tmp.PYL1.drc.in$GR) == F)
    parameters.Hill.TECAN[i,1:4] <- summary(tmp.PYL1.drc)$coefficients[,1]
    parameters.Hill.TECAN[i,5:8] <- summary(tmp.PYL1.drc)$coefficients[,2]
    parameters.Hill.TECAN[i,9:12] <- summary(tmp.PYL1.drc)$coefficients[,4]
    if(parameters.Hill.TECAN[i,1] < 0){
      parameters.Hill.TECAN[i,1] <- -parameters.Hill.TECAN[i,1] ## invert Hill parameter output
    }else{
      parameters.Hill.TECAN[i,2:3] <- parameters.Hill.TECAN[i,c(3,2)]
      parameters.Hill.TECAN[i,10:11] <- parameters.Hill.TECAN[i,c(11,10)]
    }
    parameters.Hill.TECAN[i,14:15] <- as.numeric(mselect(tmp.PYL1.drc, icfct = AIC)[c(2,4)])
    
    ### calculation of non-linear R^2
    predicted_values <- predict(tmp.PYL1.drc)
    rss <- sum((tmp.PYL1.drc.in$GR - predicted_values)^2)
    tss <- sum((tmp.PYL1.drc.in$GR - mean(tmp.PYL1.drc.in$GR))^2)
    tmp.r_squared.Hill <- 1 - (rss/tss)
    parameters.Hill.TECAN[i,16] <- tmp.r_squared.Hill
    rm(predicted_values, rss, tss, tmp.r_squared.Hill)
   
    ### predict the full, smoothened curve (using 1000 data points) and confidence interval
    tmp.PYL1.drc.Hill.predict.newdata <- expand.grid(conc = c(0, exp(seq(log(0.001), log(5000), length = 999))))
    tmp.PYL1.drc.Hill.predict <- predict(tmp.PYL1.drc, 
                                         newdata = tmp.PYL1.drc.Hill.predict.newdata, 
                                         interval = "confidence")
    
    #### new data with predictions
    tmp.PYL1.drc.Hill.predict.newdata$p <- tmp.PYL1.drc.Hill.predict[,1]
    tmp.PYL1.drc.Hill.predict.newdata$pmin <- tmp.PYL1.drc.Hill.predict[,2]
    tmp.PYL1.drc.Hill.predict.newdata$pmax <- tmp.PYL1.drc.Hill.predict[,3]
    
    ### plot curve
    tmp.PYL1.drc.in$concentration[1] <- 9.062741e-03/3.5/3.5/3.5 ## "0-conc." positioning for log scale
    # print(ggplot(tmp.PYL1.drc.in, aes(x = concentration, y = GR)) +
    #         geom_ribbon(data = tmp.PYL1.drc.Hill.predict.newdata,
    #                     aes(x = conc, y = p, ymin = pmin, ymax = pmax),
    #                     alpha = 0.2, fill = "grey50") +
    #         geom_point(data = tmp.PYL1.drc.in, aes(x = concentration, y = GR),
    #                    color = "black", size = 10) +
    #         geom_line(data = tmp.PYL1.drc.Hill.predict.newdata, aes(x = conc, y = p), linewidth = 1.5) +
    #         scale_x_log10(breaks = c(9.062741e-03/3.5/3.5/3.5, 0.01, 0.1, 1, 10, 100, 1000),
    #                       labels = c(0, 0.01, 0.1, 1, 10, 100, "1,000"),
    #                       limits = c(9.062741e-03/3.5/3.5/3.5, 5000)) +
    #         scale_y_continuous(breaks = seq(from = 0, to = 140, length.out = 8), limits = c(-20000, 20000)) +
    #         coord_cartesian(ylim = c(-5, 150)) +
    #         annotate("text",
    #                  x = 9.062741e-03/3.5/3.5/3.5,
    #                  y = 140,
    #                  label = bquote(italic(R)^2 == .(format(parameters.Hill.TECAN[i,"R^2"], digits = 2))),  ##
    #                  hjust = 0, size = 25, color = "black") +
    #         theme_classic(base_size = 50) +
    #         theme(plot.title = element_markdown(),
    #               plot.subtitle = element_markdown(size = 20),
    #               title = element_text(size = 40),
    #               axis.text = element_text(size = 40),
    #               axis.line.x = element_line(size = 1, color = 'black'),
    #               axis.line.y = element_line(size = 1, color = 'black'),
    #               axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
    #               axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),
    #               legend.position = "none",
    #               text = element_text(family="Helvetica"),
    #               plot.margin = unit(c(2, 2, 2, 2),"cm")) +
    #         labs(x = "(+)-ABA conc. (µM)",
    #              y = "Binding (library sequencing)",
    #              title = rownames(parameters.Hill.TECAN)[i],
    #              subtitle = paste0("Hill parameters: B[0] = ", round(parameters.Hill.TECAN[i,"B[0]"],3), if(!is.na(parameters.Hill.TECAN[i,"B[0] P"])){if(parameters.Hill.TECAN[i,"B[0] P"] < 0.05){"*"}},
    #                                ", B[inf] = ", round(parameters.Hill.TECAN[i,"B[inf]"],3), if(!is.na(parameters.Hill.TECAN[i,"B[inf] P"])){if(parameters.Hill.TECAN[i,"B[inf] P"] < 0.05){"*"}},
    #                                ", EC50 = ", round(parameters.Hill.TECAN[i,"EC50"],3), if(!is.na(parameters.Hill.TECAN[i,"EC50 P"])){if(parameters.Hill.TECAN[i,"EC50 P"] < 0.05){"*"}},
    #                                ", Hill = ", round(parameters.Hill.TECAN[i,"Hill"],3), if(!is.na(parameters.Hill.TECAN[i,"Hill P"])){if(parameters.Hill.TECAN[i,"Hill P"] < 0.05){"*"}})))

    rm(tmp.PYL1.drc, tmp.PYL1.drc.Hill.predict, tmp.PYL1.drc.Hill.predict.newdata, tmp.PYL1.drc.in)
    
  }

}
#dev.off()


## 4. Plot ##
#############

## summarise data
out.all <- matrix(NA, nrow = 22, ncol = 16)
colnames(out.all) <- c("bulk_B0", "bulk_B0_SE", "bulk_Binf", "bulk_Binf_SE", 
                       "bulk_EC50", "bulk_EC50_SE", "bulk_Hill", "bulk_Hill_SE",
                       "TECAN_B0", "TECAN_B0_SE", "TECAN_Binf", "TECAN_Binf_SE", 
                       "TECAN_EC50", "TECAN_EC50_SE", "TECAN_Hill", "TECAN_Hill_SE")
rownames(out.all) <- c("Q34Y", "Q34I", "E36R", "T38K", "Q39F", 
                       "H48C", "D80M", "P82G", "Y85P", "H87A",
                       "H87P", "A116R", "T118I", "L125I", "V132P", 
                       "I137G", "L144A", "E161A", "E171K", "A190E",
                       "R195M", "A202T")
for(i in 1:nrow(out.all)){
  tmp.name <- rownames(out.all)[i]
  tmp.id1 <- match(tmp.name, rownames(parameters.Hill))
  out1 <- parameters.Hill[tmp.id1,c(2,6,3,7,4,8,1,5)]
  tmp.id2 <- match(tmp.name, rownames(parameters.Hill.TECAN))
  out2 <- parameters.Hill.TECAN[tmp.id2,c(2,6,3,7,4,8,1,5)]
  out.all[i,] <- c(out1, out2)
}
out.all <- as.data.frame(out.all)

## raw linear regression
p.B0 <- summary(lm(out.all$bulk_B0 ~ out.all$TECAN_B0))$coefficients[2,4] ## P = 3.548e-05
p.Binf <- summary(lm(out.all$bulk_Binf ~ out.all$TECAN_Binf))$coefficients[2,4] ## P = 0.001079
p.EC50 <- summary(lm(log10(out.all$bulk_EC50) ~ log10(out.all$TECAN_EC50)))$coefficients[2,4] ## P = 0.002154
p.Hill <- summary(lm(log10(out.all$bulk_Hill) ~ log10(out.all$TECAN_Hill)))$coefficients[2,4] ## P = 0.6181

## B0
r.B0 <- cor(x = out.all$bulk_B0, y = out.all$TECAN_B0, 
            method = "pearson", use = "complete.obs")

out.S3B_B0 <- ggplot(out.all, aes(x = `TECAN_B0`, y = `bulk_B0`)) +
  scale_x_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  scale_y_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  coord_cartesian(xlim = c(-5, 130), ylim = c(-5, 130)) +
  geom_point(data = out.all,
             mapping = aes(x = `TECAN_B0`, y = `bulk_B0`),
             color = "black", size = 6) +
  geom_errorbarh(aes(xmin = pmax(`TECAN_B0` - `TECAN_B0_SE`, 0), 
                     xmax = `TECAN_B0` + `TECAN_B0_SE`), color = "black") +
  geom_errorbar(aes(ymin = pmax(`bulk_B0` - `bulk_B0_SE`, 0),
                    ymax = `bulk_B0` + `bulk_B0_SE`), color = "black") +
  geom_smooth(data = out.all,
              mapping = aes(x = `TECAN_B0`, y = `bulk_B0`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = NA) +
  annotate("text",
           x = 0,
           y = 120,
           label = bquote(italic(r) == .(format(r.B0, digits = 2, nsmall = 2)) ~ 
                            ", " ~ italic(P) == .(format(p.B0, digits = 2, nsmall = 2))),
           hjust = 0, size = 15, color = "black") +
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
  labs(x = bquote(B[0] ~ "(individual mutants)"),
       y = bquote(B[0] ~ "(library sequencing)"))

## Binf
r.Binf <- cor(x = out.all$bulk_Binf, y = out.all$TECAN_Binf, 
              method = "pearson", use = "complete.obs")

out.S3B_Binf <- ggplot(out.all, aes(x = `TECAN_Binf`, y = `bulk_Binf`)) +
  scale_x_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  scale_y_continuous(breaks = seq(f = 0, t = 120, length.out = 7), limits = c(-10000, 10000)) +
  coord_cartesian(xlim = c(-5, 130), ylim = c(-5, 130)) +
  geom_point(data = out.all,
             mapping = aes(x = `TECAN_Binf`, y = `bulk_Binf`),
             color = "black", size = 6) +
  geom_errorbarh(aes(xmin = pmax(`TECAN_Binf` - `TECAN_Binf_SE`, 0),
                     xmax = `TECAN_Binf` + `TECAN_Binf_SE`), color = "black") +
  geom_errorbar(aes(ymin = pmax(`bulk_Binf` - `bulk_Binf_SE`, 0),
                    ymax = `bulk_Binf` + `bulk_Binf_SE`), color = "black") +
  geom_smooth(data = out.all,
              mapping = aes(x = `TECAN_Binf`, y = `bulk_Binf`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = NA) +
  annotate("text",
           x = 0,
           y = 120,
           label = bquote(italic(r) == .(format(r.Binf, digits = 2, nsmall = 2)) ~ 
                            ", " ~ italic(P) == .(format(p.Binf, digits = 2, nsmall = 2))),
           hjust = 0, size = 15, color = "black") +
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
  labs(x = bquote(B[infinity] ~ "(individual mutants)"),
       y = bquote(B[infinity] ~ "(library sequencing)"))

## EC50
r.EC50 <- cor(x = log(out.all$bulk_EC50), y = log(out.all$TECAN_EC50), 
              method = "pearson", use = "complete.obs")

out.S3B_EC50 <- ggplot(out.all, aes(x = `TECAN_EC50`, y = `bulk_EC50`)) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000), 
                labels = c(0.01, 0.1, 1, 10, 100, "1,000"),
                limits = c(0.0000000001,1000000)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000), 
                labels = c(0.01, 0.1, 1, 10, 100, "1,000"),
                limits = c(0.0000000001,1000000)) +
  coord_cartesian(xlim = c(0.001, 10000), ylim = c(0.001, 10000), expand = T) +
  geom_point(data = out.all,
             mapping = aes(x = `TECAN_EC50`, y = `bulk_EC50`),
             color = "black", size = 6) +
  geom_errorbarh(aes(xmin = pmax(TECAN_EC50 - TECAN_EC50_SE, 0.001),
                     xmax = TECAN_EC50 + TECAN_EC50_SE), color = "black") +
  geom_errorbar(aes(ymin = pmax(`bulk_EC50` - `bulk_EC50_SE`, 0.001),
                    ymax = `bulk_EC50` + `bulk_EC50_SE`), color = "black") +
  geom_smooth(data = out.all,
              mapping = aes(x = `TECAN_EC50`, y = `bulk_EC50`),
              method = 'lm',
              color = "black", fullrange = T,
              linewidth = 1.5, alpha = 0.2, fill = NA) +
  annotate("text",
           x = 0.0018,
           y = 3025,
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
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm")) +
  labs(x = bquote(EC[50] ~ "(individual mutants)"),
       y = bquote(EC[50] ~ "(library sequencing)"))

## n
r.Hill <- cor(x = log(out.all$bulk_Hill), y = log(out.all$TECAN_Hill), 
              method = "pearson", use = "complete.obs")

out.S3B_Hill <- ggplot(out.all, aes(x = `TECAN_Hill`, y = `bulk_Hill`)) +
  scale_x_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5), 
                labels = c(0.1, 0.2, 0.5, 1, 2, 5),
                limits = c(0.0000000001,1000000)) +
  scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5), 
                labels = c(0.1, 0.2, 0.5, 1, 2, 5),
                limits = c(0.0000000001,1000000)) +
  coord_cartesian(xlim = c(0.1, 10), ylim = c(0.1, 10), expand = T) +
  geom_point(data = out.all,
             mapping = aes(x = `TECAN_Hill`, y = `bulk_Hill`),
             color = "black", size = 6) +
  geom_errorbarh(aes(xmin = pmax(TECAN_Hill - TECAN_Hill_SE, 0.1),
                     xmax = TECAN_Hill + TECAN_Hill_SE), color = "black") +
  geom_errorbar(aes(ymin = pmax(`bulk_Hill` - `bulk_Hill_SE`, 0.1),
                    ymax = `bulk_Hill` + `bulk_Hill_SE`), color = "black") +
  annotate("text",
           x = 0.1186,
           y = 7.1,
           label = bquote(italic(r) == .(format(r.Hill, digits = 2, nsmall = 2)) ~ 
                            ", " ~ italic(P) == .(format(p.Hill, digits = 4, nsmall = 4))),
           hjust = 0, size = 15, color = "black") +
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
  labs(x = bquote(n ~ "(individual mutants)"),
       y = bquote(n ~ "(library sequencing)"))

## combined plot
pdf("../../results/FigureS3/FigureS3B_Hill_parameter_validation.pdf", height = 15, width = 55)
plot_grid(out.S3B_B0, out.S3B_Binf, out.S3B_EC50, out.S3B_Hill, align = "hv", axis = "tblr", ncol = 4)
dev.off()


## 4. Version ##
################

