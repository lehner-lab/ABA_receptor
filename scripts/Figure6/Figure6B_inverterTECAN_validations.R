# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

####################################################################
## Figure 6B - TECAN validations of curves with inverted response ##
####################################################################


## 0. Environment ##
####################

## Libraries
packages <- c("stringr", "readxl", "growthrates", "drc", 
              "scales", "reshape", "ggplot2", "ggtext", "cowplot",
              "wesanderson")

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
colnames(setup) <- 0:24
rownames(setup) <- setup[,1]
setup <- setup[,-1]
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
PYL1.ABI1.TECAN.muts <- PYL1.ABI1.TECAN.muts[c("WT", "V110Y", "V110W", "A116W", "H142W")]
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
  
    rm(tmp.PYL1.drc, tmp.PYL1.drc.Hill.predict, tmp.PYL1.drc.Hill.predict.newdata, tmp.PYL1.drc.in)

}


## 4. Calculate the curves and confidence intervals ##
######################################################

## Define the four-parameter logistic function
four_pl_function <- function(x, Hill, B0, Binf, EC50) {
  B0 + (Binf - B0) / (1 + (EC50 / x)^c(Hill))
}

## Calculate the predicted response over the dose range
curves <- matrix(NA, nrow = nrow(parameters.Hill.TECAN) - 2, ncol = 1000)
rownames(curves) <- rownames(parameters.Hill.TECAN)[-c(1:2)]
colnames(curves) <- c(0, exp(seq(log(0.001), log(5000), length = 999)))
curves.min <- curves.max <- curves

dosages <- rep(NA, 12)
dosages[1] <- 2500
dosages[12] <- 0
for(i in 2:11){
  dosages[i] <- 2500 * (1/3.5)^c(i-1)
}

for(i in 1:nrow(curves)){
  
  print(i)
  
  ## fetch raw data
  tmp.id <- rownames(curves)[i]
  
  ## estimate confidence interval
  tmp.df <- as.data.frame(cbind("out" = PYL1.ABI1.TECAN.muts.rates[tmp.id][[1]],
                                "conc" = rev(dosages)))
  tmp.model <- drm(tmp.df$out ~ tmp.df$conc, 
                   fct = LL.4(fixed = c(NA, NA, NA, NA)),
                   type = 'continuous') 
  tmp.predict.newdata <- expand.grid(conc = c(0, exp(seq(log(0.001), log(5000), length = 999))))
  tmp.model.predict <- predict(tmp.model,
                               newdata = tmp.predict.newdata,
                               interval = "confidence")
  
  ## fetch curves and confidences interval
  curves[i,] <- tmp.model.predict[,"Prediction"]
  curves.min[i,] <- tmp.model.predict[,"Lower"]
  curves.max[i,] <- tmp.model.predict[,"Upper"]
  
}

## mean
curves <- rbind(curves,
                "conc" = c(0, exp(seq(log(0.001), log(5000), length = 999))))
curves <- t(curves)
curves <- as.data.frame(curves)
curves <- melt(curves, id.vars = "conc")

## max
curves.max <- rbind(curves.max,
                    "conc" = c(0, exp(seq(log(0.001), log(5000), length = 999))))
curves.max <- t(curves.max)
curves.max <- as.data.frame(curves.max)
curves.max <- melt(curves.max, id.vars = "conc")

## min
curves.min <- rbind(curves.min,
                    "conc" = c(0, exp(seq(log(0.001), log(5000), length = 999))))
curves.min <- t(curves.min)
curves.min <- as.data.frame(curves.min)
curves.min <- melt(curves.min, id.vars = "conc")

## clean up
rm(i, four_pl_function, tmp.df, tmp.model, tmp.model.predict, tmp.predict.newdata, tmp.id)


## 5. Plot ##
#############

### curves
inverter.plot <- function(id){
  
  ## extract curve fit and 95% confidence interval
  curve.tmp <- curves[grep(id, curves$variable),]
  curve.tmp$min.value <- curves.min[grep(id, curves.min$variable),"value"]
  curve.tmp$max.value <- curves.max[grep(id, curves.max$variable),"value"]
  curve.tmp$variable <- factor(curve.tmp$variable)
  
  ## extract data points
  points.tmp <- rbind(rev(PYL1.ABI1.TECAN.muts.rates[id][[1]]))
  points.tmp <- as.data.frame(t(points.tmp))
  dosages[12] <- 9.062741e-03/3.5/3.5/3.5
  points.tmp$conc <- dosages
  colnames(points.tmp)[1] <- "value"
  
  ## plot
  cols <- substr(id, nchar(id), nchar(id))
  if(cols %in% c("F", "Y", "W")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[5]
  }else if(cols %in% c("A", "V", "L", "I", "M")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[4]
  }else if(cols %in% c("Q", "S", "T", "N", "C")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[3]
  }else if(cols %in% c("K", "R", "H")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[2]
  }else if(cols %in% c("D", "E")){
    cols <- wes_palette("Darjeeling2", 6, type = "continuous")[1]
  }else if(cols %in% c("P", "G")){
    cols <- "black"
  }
  
  p <- ggplot(data = curve.tmp, aes(x = conc, y = variable)) +
    geom_ribbon(data = curve.tmp,
                aes(group = variable, x = conc, y = value, ymin = min.value, ymax = max.value),
                alpha = 0.1, fill = "grey50") +
    geom_line(data = curve.tmp,
              aes(group = variable, x = conc, y = value), 
              linewidth = 3, alpha = 1, colour = cols) +
    geom_point(data = points.tmp, aes(x = conc, y = value),
               size = 10, alpha = 1, colour = cols) +
    annotate("text", x = 0.0002113759, y = 30,
             label = bquote(italic(R)^2 == .(format(round(parameters.Hill.TECAN[grep(id, rownames(parameters.Hill.TECAN)), "R^2"], 2), digits = 2))),
    color = "black", hjust = 0, size = 17) +
    annotate("text", x = 0.0002113759, y = 10,
      label = bquote(Delta * "Binding = " * .(format(round(parameters.Hill.TECAN[grep(id, rownames(parameters.Hill.TECAN)), "B[inf]"] -
                                                           parameters.Hill.TECAN[grep(id, rownames(parameters.Hill.TECAN)),"B[0]"], 1), nsmall = 1)) * "%"),
      color = "black", hjust = 0, size = 17) +
    scale_x_log10(breaks = c(9.062741e-03/3.5/3.5/3.5, 0.01, 0.1, 1, 10, 100, 1000),
                  labels = c(0, 0.01, 0.1, 1, 10, 100, "1,000"),
                  limits = c(9.062741e-03/3.5/3.5/3.5, 5000)) +
    scale_y_continuous(breaks = seq(from = 0, to = 150, length.out = 6), limits = c(-10000, 10000)) +
    coord_cartesian(ylim = c(-5, 150)) +
    theme_classic(base_size = 50) +
    theme(plot.title = element_markdown(size = 60),
          axis.text = element_text(size = 25),
          axis.line.x = element_line(linewidth = 1, color = 'black'),
          axis.line.y = element_line(linewidth = 1, color = 'black'),
          axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
          axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = 3),
          legend.position = "none",
          text = element_text(family="Helvetica"),
          plot.margin = unit(c(1, 1, 2, 2),"cm")) +
    labs(title = id,
         x = "(+)-ABA conc. (µM)",
         y = "Binding (individual mutants)")
  return(p)
  
}
out.H142W <- inverter.plot("H142W")
out.V110W <- inverter.plot("V110W")
out.A116W <- inverter.plot("A116W")

### combine in one plot
pdf("../../results/Figure6/Figure6B_inverter_TECANs.pdf", width = 29, height = 12)

plot_grid(out.H142W, out.V110W, out.A116W, align = "hv", axis = "tblr", ncol = 3)

dev.off()


## 6. Version ##
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
# [1] wesanderson_0.3.7 cowplot_1.1.3     ggtext_0.1.2      ggplot2_3.5.1     reshape_0.8.9     scales_1.3.0     
# [7] drc_3.0-1         MASS_7.3-64       growthrates_0.8.4 deSolve_1.40      lattice_0.22-6    readxl_1.4.3     
# [13] stringr_1.5.1    
# 
# loaded via a namespace (and not attached):
# [1] generics_0.1.3    sandwich_3.1-1    xml2_1.3.6        gtools_3.9.5      stringi_1.8.4     FME_1.3.6.3      
# [7] magrittr_2.0.3    grid_4.4.1        mvtnorm_1.3-3     cellranger_1.1.0  plyr_1.8.9        Matrix_1.7-2     
# [13] Formula_1.2-5     survival_3.8-3    multcomp_1.4-28   TH.data_1.1-3     codetools_0.2-20  abind_1.4-8      
# [19] cli_3.6.4         rlang_1.1.5       commonmark_1.9.2  munsell_0.5.1     splines_4.4.1     withr_3.0.2      
# [25] plotrix_3.8-4     rootSolve_1.8.2.4 tools_4.4.1       parallel_4.4.1    coda_0.19-4.1     minpack.lm_1.2-4 
# [31] minqa_1.2.8       dplyr_1.1.4       colorspace_2.1-1  vctrs_0.6.5       R6_2.6.1          zoo_1.8-12       
# [37] lifecycle_1.0.4   car_3.1-3         pkgconfig_2.0.3   pillar_1.10.1     gtable_0.3.6      glue_1.8.0       
# [43] Rcpp_1.0.14       xfun_0.51         tidyselect_1.2.1  tibble_3.2.1      rstudioapi_0.17.1 farver_2.1.2     
# [49] carData_3.0-5     compiler_4.4.1    markdown_1.13     gridtext_0.1.5   
