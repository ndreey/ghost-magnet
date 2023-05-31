library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggpubr)
library(pheatmap)
source("analysis_tools.R")

################################################################################
# Analyses the genome fraction
################################################################################

# STANDARD VECTORS
# HC percentages
hc <- c(0, 0.254, 0.429, 0.566, 0.664, 0.759, 0.818, 0.878, 0.923, 0.964, 0.983)
# Sample names
sample_id <- c("X00", "X01", "X02", "X03", "X04", "X05", "X06", "X07",
               "X08", "X090", "X095")
hc2 <- as.character(hc*100)


# Load in data
# GSA
gfrac <- read.table("assmble/gsa/gold/TSV/gsa_Genome_fraction.tsv",
                    sep = "\t", header = TRUE, row.names = 1,
                    na.strings = "-")
# Replace NAs
gfrac[is.na(gfrac)] <- 0

# Clean up data and add columns
gfrac <- clean_data(gfrac, "gsa")


# Index out the groups
endo <- gfrac[gfrac$group != "host",1:11]
host <- gfrac[gfrac$group == "host",1:11]
mock <- gfrac[,1:11]
rf <- gfrac[gfrac$group == "rfungi",1:11]
omf <- gfrac[gfrac$group == "OMF",1:11]
fungi <- gfrac[gfrac$group %in% c("rfungi", "OMF"), 1:11]
ba_ar <- gfrac[gfrac$group == "ba_ar",1:11]
pl_vi_unk <- gfrac[gfrac$group == "pl_vi_unk",1:11]


# Boxplot
boxplot(mock, main = "Genome Fraction")

m2 <- mock*(1/100)

range(m2)

library(viridis)


# Heatmap all genomes
pheatmap(m2, scale = "row", 
         col = viridis(256, begin = 0, end = 1, option = "G"),
         cluster_cols = F, 
         cluster_rows = T, 
         fontsize_row = 7,
         fontsize_col = 9,
         angle_col = 0, 
         cutree_rows = 2,
         treeheight_row = 0,
         labels_col = hc2,
         main = "Genome Fraction",
         xlab = "Host Contamination (%)")

?pheatmap
# Heat map, mock groups
gsa_host <- host
omf_med <- apply(omf, 2, median)
rf_med <- apply(rf, 2, median)
baar_med <- apply(ba_ar, 2, median)
pl_med <- apply(pl_vi_unk, 2, median)
endo_med <- apply(endo, 2, median)

heat_group <- rbind(host = gsa_host, OMF = omf_med, rfungi = rf_med, ba_ar = baar_med,
                    pl_vi_unk = pl_med, endophytes = endo_med)

# hatmap with endophytes
pheatmap(heat_group, 
         scale = "row",
         cluster_cols = F,
         angle_col = 315,
         col = brewer.pal(9,"PuBu"),
         treeheight_row = 0,
         cutree_rows = 2,
         main = "Genome Fraction",
         fontsize_row = 12,
         fontsize_col = 12,
         cellwidth = 30,
         cellheight = 100)

# Heatmap without endophyte
pheatmap(heat_group[1:5,], 
         scale = "row",
         cluster_cols = F,
         cluster_rows = T,
         angle_col = 315,
         col = brewer.pal(9,"PuBu"),
         treeheight_row = 0,
         cutree_rows = 2,
         gaps_row = c(1, 5),
         main = "Genome Fraction",
         fontsize_row = 12,
         fontsize_col = 12,
         cellwidth = 30,
         cellheight = 100)

################## Analysis #######################
# GSA
median <- apply(mock, 2, median, na.rm = TRUE)
shapiro.test(median) # W = 0.80552, p-value = 0.01112 Not norm..
# Correlation
cor.test(hc, median, method ="kendall", exact = FALSE) 
#p-value = 5.511e-07 r = -0.9636364 




################################# RAW ##########################################
# Load in data
# raw
raw_gfrac <- read.table("assmble/gsa/raw_assembly/summary/TSV/raw_Genome_fraction.tsv",
                    sep = "\t", header = TRUE, row.names = 1,
                    na.strings = "-")
# Replace NAs
raw_gfrac[is.na(raw_gfrac)] <- 0

# Clean up data and add columns
raw_gfrac <- clean_data(raw_gfrac, "raw")


# Index out the groups
r_endo <- raw_gfrac[raw_gfrac$group != "host",1:11]
r_host <- raw_gfrac[raw_gfrac$group == "host",1:11]
r_mock <- raw_gfrac[,1:11]
r_rf <- raw_gfrac[raw_gfrac$group == "rfungi",1:11]
r_omf <- raw_gfrac[raw_gfrac$group == "OMF",1:11]
r_ba_ar <- raw_gfrac[raw_gfrac$group == "ba_ar",1:11]
r_pl_vi_unk <- raw_gfrac[raw_gfrac$group == "pl_vi_unk",1:11]


# Boxplot
boxplot(mock, main = "Genome Fraction")



# Heat map, mock groups
r_median <- apply(r_mock, 2, median)
raw_host <- r_host
r_omf_med <- apply(r_omf, 2, median)
r_rf_med <- apply(r_rf, 2, median)
r_baar_med <- apply(r_ba_ar, 2, median)
r_pl_med <- apply(r_pl_vi_unk, 2, median)
r_endo_med <- apply(r_endo, 2, median)

r_heat_group <- rbind(host = raw_host,
                      rfungi = r_rf_med,
                      ba_ar = r_baar_med,
                      OMF = r_omf_med,
                      pl_vi_unk = r_pl_med,
                      endophytes = r_endo_med)
# Heatmap with endophytes as their own group as well
pheatmap(r_heat_group, 
         scale = "row",
         cluster_cols = F,
         cluster_rows = T,
         angle_col = 315,
         col = brewer.pal(9,"PuBu"),
         treeheight_row = 0,
         cutree_rows = 2,
         gaps_row = c(1, 5),
         main = "Genome Fraction",
         fontsize_row = 12,
         fontsize_col = 12,
         cellwidth = 30,
         cellheight = 100)

# Heatmap without endophyte
pheatmap(r_heat_group[1:5,], 
         scale = "row",
         cluster_cols = F,
         cluster_rows = T,
         angle_col = 315,
         col = brewer.pal(9,"PuBu"),
         treeheight_row = 0,
         cutree_rows = 2,
         gaps_row = c(1, 5),
         main = "Genome Fraction",
         fontsize_row = 12,
         fontsize_col = 12,
         cellwidth = 30,
         cellheight = 100)

################## Analysis #######################
# Seeing if the whole mock dataset has a -correlation
shapiro.test(r_median) # W = 0.54954, p-value = 6.062e-06 Not norm..
# Correlation
cor.test(hc, r_median, method ="kendall", exact = FALSE) 
#p-value = 0.001152 r = -0.8101627  

# Checking just endos # W = 0.53508, p-value = 4.036e-06 not norm
shapiro.test(r_endo_med)
cor.test(hc, r_endo_med, method = "kendall", exact = F)
# p-value = 0.000813  r = -0.8420754



