# The genetic architecture of an allosteric hormone receptor
# Maximilian R. Stammnitz & Ben Lehner
# bioRxiv link: https://www.biorxiv.org/content/10.1101/2025.05.30.656975v1
# 31.05.2025
# © M.R.S. (maximilian.stammnitz@crg.eu)

#######################################################
## Supplementary Figure 1 - GluePCA TECAN benchmarks ##
#######################################################


## 0. Environment ##
####################

## Libraries
packages <- c("readxl", "growthrates", "drc", "scales")

## Install missing packages
install_if_missing <- function(pkg) {if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)}
lapply(packages, install_if_missing)

## Load libraries
lapply(packages, library, character.only = TRUE)


## 1. Import TECAN data ##
##########################

## PYL1-ABI1
pyl1.abi1 <- as.matrix(read_xlsx('../../data/TECAN/GluePCA_systems/PYL1-ABI1/GluePCA_benchmark_PYL1_ABI1_TECAN.xlsx', sheet = 1))
rownames(pyl1.abi1) <- pyl1.abi1[,1]
pyl1.abi1 <- pyl1.abi1[,-c(1:2)]
pyl1.abi1 <- pyl1.abi1[-1,]
colnames(pyl1.abi1) <- as.numeric(colnames(pyl1.abi1))/3600 ## convert to hours
class(pyl1.abi1) <- "numeric"
pyl1.abi1 <- as.data.frame(pyl1.abi1)

## GID1A-GAI ##
gid1a.gai <- as.matrix(read_xlsx('../../data/TECAN/GluePCA_systems/GID1A-GAI/GluePCA_benchmark_GID1A-GAI_TECAN.xlsx', sheet = 1))
rownames(gid1a.gai) <- gid1a.gai[,1]
gid1a.gai <- gid1a.gai[,-c(1:2)]
gid1a.gai <- gid1a.gai[-1,]
colnames(gid1a.gai) <- as.numeric(colnames(gid1a.gai))/3600 ## convert to hours
class(gid1a.gai) <- "numeric"
gid1a.gai <- as.data.frame(gid1a.gai)

## FKBP12-FRB ##
fkbp12.frb <- as.matrix(read_xlsx('../../data/TECAN/GluePCA_systems/FKBP12-FRB/GluePCA_benchmark_FKBP12_FRB_TECAN.xlsx', sheet = 1))
rownames(fkbp12.frb) <- fkbp12.frb[,1]
fkbp12.frb <- fkbp12.frb[,-c(1:2)]
fkbp12.frb <- fkbp12.frb[-1,]
colnames(fkbp12.frb) <- as.numeric(colnames(fkbp12.frb))/3600 ## convert to hours
class(fkbp12.frb) <- "numeric"
fkbp12.frb <- as.data.frame(fkbp12.frb)


## 2. Plot ##
#############

## PYL1-ABI1
pdf("../../results/FigureS1/FigureS1A_GluePCA_benchmark_PYL1-ABI1.pdf", width = 15, height = 8)
par(mar = c(6,8,3,0.5), mfrow = c(1,2))

### WT

#### without ABA
plot(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep1",], type = "p", pch = 16, cex = 0.6, bty = "n", 
     ylim = c(0.1,0.6), xaxt = "n", yaxt = "n", ylab = "", xlab = "", main = "WT", cex.main = 2.8, xlim = c(0,70))
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep1",])
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep2",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep2",])
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep3",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep3",])
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep4",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_0nM_rep4",])

#### with ABA
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep1",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep1",], col = "darkgreen")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep2",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep2",], col = "darkgreen")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep3",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep3",], col = "darkgreen")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep4",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["WT_380µM_rep4",], col = "darkgreen")

#### axes, etc
axis(1, at = c(0, 10, 20, 30, 40, 50, 60, 70), labels = c(0, 10, 20, 30, 40, 50, 60, 70), cex.axis = 2)
mtext("Time [hours]", side = 1, line = 4.5, cex = 2.5)
axis(2, at = seq(f = 0.1, t = 0.6, by = 0.1), labels = seq(f = 0.1, t = 0.6, by = 0.1), cex.axis = 2, las = 2)
mtext("OD", side = 2, line = 4.5, cex = 2.5)
legend("topleft", legend = c("380 µM (+)-ABA", "no ABA"), pch = 16, 
       col = c("darkgreen", "black"), bty = "n", cex = 1.4)


### PYL1-ABI1 (ABI1 W300A)

#### without ABA
plot(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep1",], type = "p", pch = 16, cex = 0.6, bty = "n", 
     ylim = c(0.1,0.6), xaxt = "n", yaxt = "n", ylab = "", xlab = "", col = "black", main = "ABI1 W300A", cex.main = 2.8, 
     xlim = c(0,70))
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep1",], col = "black")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep2",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep2",], col = "black")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep3",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep3",], col = "black")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep4",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_0nM_rep4",], col = "black")

#### with ABA
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep1",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep1",], col = "darkgreen")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep2",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep2",], col = "darkgreen")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep3",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep3",], col = "darkgreen")
points(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep4",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(pyl1.abi1)), y = pyl1.abi1["ABI1_W300A_380µM_rep4",], col = "darkgreen")

#### axes, etc
axis(1, at = c(0, 10, 20, 30, 40, 50, 60, 70), labels = c(0, 10, 20, 30, 40, 50, 60, 70), cex.axis = 2)
mtext("Time [hours]", side = 1, line = 4.5, cex = 2.5)
axis(2, at = seq(f = 0.1, t = 0.6, by = 0.1), labels = seq(f = 0.1, t = 0.6, by = 0.1), cex.axis = 2, las = 2)
mtext("OD", side = 2, line = 4.5, cex = 2.5)

dev.off()


## GID1A-GAI ##
pdf("../../results/FigureS1/FigureS1B_GluePCA_benchmark_GID1A.pdf", width = 15, height = 8)
par(mar = c(6,8,3,0.5), mfrow = c(1,2))

### WT

#### without GA3
plot(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep1",], type = "p", pch = 16, cex = 0.6, bty = "n", 
     ylim = c(0.15,1.2), xaxt = "n", yaxt = "n", ylab = "", xlab = "", main = "WT", cex.main = 2.8, 
     xlim = c(0,70))
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep1",])
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep2",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep2",])
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep3",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep3",])
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep4",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_0nM_rep4",])

#### with GA3
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep1",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep1",], col = "darkgreen")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep2",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep2",], col = "darkgreen")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep3",],pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep3",], col = "darkgreen")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep4",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["WT_1mM_rep4",], col = "darkgreen")

#### axes, etc
axis(1, at = c(0, 10, 20, 30, 40, 50, 60, 70), labels = c(0, 10, 20, 30, 40, 50, 60, 70), cex.axis = 2)
mtext("Time [hours]", side = 1, line = 4.5, cex = 2.5)
axis(2, at = seq(f = 0.2, t = 1.2, by = 0.2), labels = seq(f = 0.2, t = 1.2, by = 0.2), cex.axis = 2, las = 2)
mtext("OD", side = 2, line = 4.5, cex = 2.5)
legend("topleft", legend = c("1 mM GA3", "no GA3"), pch = 16, 
       col = c("darkgreen", "black"), bty = "n", cex = 1.4)


### GID1A-GAI (GAI E54R)

#### without GA3
plot(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep1",], type = "p", pch = 16, cex = 0.6, bty = "n", 
     ylim = c(0.15,1.2), xaxt = "n", yaxt = "n", ylab = "", xlab = "", col = "black", 
     main = "GAI E54R", cex.main = 2.8, xlim = c(0,70))
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep1",], col = "black")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep2",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep2",], col = "black")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep3",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep3",], col = "black")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep4",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_0nM_rep4",], col = "black")

#### with GA3
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep1",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep1",], col = "darkgreen")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep2",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep2",], col = "darkgreen")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep3",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep3",], col = "darkgreen")
points(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep4",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(gid1a.gai)), y = gid1a.gai["GAI_E54R_1mM_rep4",], col = "darkgreen")

#### axes, etc
axis(1, at = c(0, 10, 20, 30, 40, 50, 60, 70), labels = c(0, 10, 20, 30, 40, 50, 60, 70), cex.axis = 2)
mtext("Time [hours]", side = 1, line = 4.5, cex = 2.5)
axis(2, at = seq(f = 0.2, t = 1.2, by = 0.2), labels = seq(f = 0.2, t = 1.2, by = 0.2), cex.axis = 2, las = 2)
mtext("OD", side = 2, line = 4.5, cex = 2.5)

dev.off()


## FKBP12-FRB ##
pdf("../../results/FigureS1/FigureS1C_GluePCA_benchmark_FRB.pdf", width = 15, height = 8)
par(mar = c(6,8,3,0.5), mfrow = c(1,2))

### WT

#### without Rapamycin
plot(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep1",], type = "p", pch = 16, cex = 0.6, bty = "n", 
     ylim = c(0.1,0.5), xaxt = "n", yaxt = "n", ylab = "", xlab = "", main = "WT", cex.main = 2.8, xlim = c(0,70))
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep1",])
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep2",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep2",])
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep3",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep3",])
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep4",], pch = 16, cex = 0.6)
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_0nM_rep4",])

#### with Rapamycin
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uM_rep1",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uM_rep1",], col = "darkgreen")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uM_rep2",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uM_rep2",], col = "darkgreen")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uM_rep3",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uM_rep3",], col = "darkgreen")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uM_rep4",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["WT_10uµM_rep4",], col = "darkgreen")

#### axes, etc
axis(1, at = c(0, 10, 20, 30, 40, 50, 60, 70), labels = c(0, 10, 20, 30, 40, 50, 60, 70), cex.axis = 2)
mtext("Time [hours]", side = 1, line = 4.5, cex = 2.5)
axis(2, at = seq(f = 0.1, t = 0.5, by = 0.1), labels = seq(f = 0.1, t = 0.5, by = 0.1), cex.axis = 2, las = 2)
mtext("OD", side = 2, line = 4.5, cex = 2.5)
legend("topleft", legend = c("10 µM Rapamycin", "no Rapamycin"), pch = 16, 
       col = c("darkgreen", "black"), bty = "n", cex = 1.4)


### FKBP12-FRB (FRB PLF + Y2088A)

#### without Rapamycin
plot(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep1",], type = "p", pch = 16, cex = 0.6, bty = "n", 
     ylim = c(0.1,0.5), xaxt = "n", yaxt = "n", ylab = "", xlab = "", col = "black", main = "FRB PLF Y2088A", cex.main = 2.8, 
     xlim = c(0,70))
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep1",], col = "black")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep2",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep2",], col = "black")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep3",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep3",], col = "black")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep4",], pch = 16, cex = 0.6, col = "black")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_0nM_rep4",], col = "black")

#### with Rapamycin
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep1",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep1",], col = "darkgreen")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep2",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep2",], col = "darkgreen")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep3",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep3",], col = "darkgreen")
points(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep4",], pch = 16, cex = 0.6, col = "darkgreen")
lines(x = as.numeric(colnames(fkbp12.frb)), y = fkbp12.frb["FRB_PLF_Y2088A_10uM_rep4",], col = "darkgreen")

#### axes, etc
axis(1, at = c(0, 10, 20, 30, 40, 50, 60, 70), labels = c(0, 10, 20, 30, 40, 50, 60, 70), cex.axis = 2)
mtext("Time [hours]", side = 1, line = 4.5, cex = 2.5)
axis(2, at = seq(f = 0.1, t = 0.5, by = 0.1), labels = seq(f = 0.1, t = 0.5, by = 0.1), cex.axis = 2, las = 2)
mtext("OD", side = 2, line = 4.5, cex = 2.5)

dev.off()


## 3. Version ##
################

