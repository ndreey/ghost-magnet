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

# Load in data
# raw
raw_gfract <- read.table("assmble/gsa/raw_assembly/summary/TSV/raw_Genome_fraction.tsv",
                    sep = "\t", header = TRUE, row.names = 1,
                    na.strings = "-")
# Replace NAs
raw_gfract[is.na(raw_gfract)] <- 0

# Clean up data and add columns
raw_gfract <- clean_data(raw_gfract, "raw")


# Index out the groups
raw_endo <- raw_gfract[raw_gfract$group != "host",1:11]
raw_host <- raw_gfract[raw_gfract$group == "host",1:11]
raw_mock <- raw_gfract[,1:11]
raw_rf <- raw_gfract[raw_gfract$group == "rfungi",1:11]
raw_omf <- raw_gfract[raw_gfract$group == "OMF",1:11]
raw_fungi <- raw_gfract[raw_gfract$group %in% c("rfungi", "OMF"), 1:11]
raw_ba_ar <- raw_gfract[raw_gfract$group == "ba_ar",1:11]
raw_pl_vi_unk <- raw_gfract[raw_gfract$group == "pl_vi_unk",1:11]



# Load in data
# x1000c
x1000c_gfract <- read.table("assmble/x1000c_quast/Genome_fraction.tsv",
                         sep = "", header = TRUE, row.names = 1,
                         na.strings = "-")
# Replace NAs
x1000c_gfract[is.na(x1000c_gfract)] <- 0


# Addin missing samples to filtered.
x1000c_gfract$X01 <- gfrac$X01_sensi.contigs
x1000c_gfract$X00 <- gfrac$X00_sensi.contigs
x1000c_gfract <- x1000c_gfract[,c(11,10,1:9)]


# Clean up data and add columns
x1000c_gfract <- clean_data(x1000c_gfract, "x1000c")


# Index out the groups
x1000c_endo <- x1000c_gfract[x1000c_gfract$group != "host",1:11]
x1000c_host <- x1000c_gfract[x1000c_gfract$group == "host",1:11]
x1000c_mock <- x1000c_gfract[,1:11]
x1000c_rf <- x1000c_gfract[x1000c_gfract$group == "rfungi",1:11]
x1000c_omf <- x1000c_gfract[x1000c_gfract$group == "OMF",1:11]
x1000c_fungi <- x1000c_gfract[x1000c_gfract$group %in% c("rfungi", "OMF"), 1:11]
x1000c_ba_ar <- x1000c_gfract[x1000c_gfract$group == "ba_ar",1:11]
x1000c_pl_vi_unk <- x1000c_gfract[x1000c_gfract$group == "pl_vi_unk",1:11]




# Load in data
# x700c
x700c_gfract <- read.table("assmble/x700c_quast/Genome_fraction.tsv",
                            sep = "", header = TRUE, row.names = 1,
                            na.strings = "-")
# Replace NAs
x700c_gfract[is.na(x700c_gfract)] <- 0


# Addin missing samples to filtered.
x700c_gfract$X00 <- gfrac$X00_sensi.contigs
x700c_gfract <- x700c_gfract[,c(11,1:10)]


# Clean up data and add columns
x700c_gfract <- clean_data(x700c_gfract, "x700c")


# Index out the groups
x700c_endo <- x700c_gfract[x700c_gfract$group != "host",1:11]
x700c_host <- x700c_gfract[x700c_gfract$group == "host",1:11]
x700c_mock <- x700c_gfract[,1:11]
x700c_rf <- x700c_gfract[x700c_gfract$group == "rfungi",1:11]
x700c_omf <- x700c_gfract[x700c_gfract$group == "OMF",1:11]
x700c_fungi <- x700c_gfract[x700c_gfract$group %in% c("rfungi", "OMF"), 1:11]
x700c_ba_ar <- x700c_gfract[x700c_gfract$group == "ba_ar",1:11]
x700c_pl_vi_unk <- x700c_gfract[x700c_gfract$group == "pl_vi_unk",1:11]





# Load in data
# x700c
maxgfract <- read.table("assmble/maxbin_quast/Genome_fraction.tsv",
                           sep = "", header = TRUE, row.names = 1,
                           na.strings = "-")
# Replace NAs
maxgfract[is.na(maxgfract)] <- 0


# Addin missing samples to filtered.
maxgfract$X00 <- gfrac$X00_sensi.contigs
maxgfract <- maxgfract[,c(11,1:10)]


# Clean up data and add columns
maxgfract <- clean_data(maxgfract, "x700c")


# Index out the groups
max_endo <- maxgfract[maxgfract$group != "host",1:11]
max_host <- maxgfract[maxgfract$group == "host",1:11]
max_mock <- maxgfract[,1:11]
max_rf <- maxgfract[maxgfract$group == "rfungi",1:11]
max_omf <- maxgfract[maxgfract$group == "OMF",1:11]
max_fungi <- maxgfract[maxgfract$group %in% c("rfungi", "OMF"), 1:11]
max_ba_ar <- maxgfract[maxgfract$group == "ba_ar",1:11]
max_pl_vi_unk <- maxgfract[maxgfract$group == "pl_vi_unk",1:11]



par(mfrow = c(2, 2))

# Boxplot
boxplot(mock, main = "700 Genome Fraction")



# Heatmap all genomes
pheatmap(raw_mock, 
         scale = "none", 
         col=brewer.pal(9,"PuBu"), 
         cluster_cols = FALSE, 
         cluster_rows = T, 
         fontsize_row = 7,
         fontsize_col = 9,
         angle_col = 0, 
         cutree_rows = 2,
         treeheight_row = 0,
         main = "Genome Fraction")



# Heat map, mock groups
gsa_host <- raw_host
omf_med <- apply(raw_omf, 2, median)
rf_med <- apply(raw_rf, 2, median)
baar_med <- apply(raw_ba_ar, 2, median)
pl_med <- apply(raw_pl_vi_unk, 2, median)
endo_med <- apply(raw_endo, 2, median)

# Heat map, mock groups
gsa_host <- x1000c_host
omf_med <- apply(x1000c_omf, 2, median)
rf_med <- apply(x1000c_rf, 2, median)
baar_med <- apply(x1000c_ba_ar, 2, median)
pl_med <- apply(x1000c_pl_vi_unk, 2, median)
endo_med <- apply(x1000c_endo, 2, median)

# Heat map, mock groups
gsa_host <- x700c_host
omf_med <- apply(x700c_omf, 2, median)
rf_med <- apply(x700c_rf, 2, median)
baar_med <- apply(x700c_ba_ar, 2, median)
pl_med <- apply(x700c_pl_vi_unk, 2, median)
endo_med <- apply(x700c_endo, 2, median)

# Heat map, mock groups
gsa_host <- max_host
omf_med <- apply(max_omf, 2, median)
rf_med <- apply(max_rf, 2, median)
baar_med <- apply(max_ba_ar, 2, median)
pl_med <- apply(max_pl_vi_unk, 2, median)
endo_med <- apply(max_endo, 2, median)

heat_group <- rbind(host = gsa_host, OMF = omf_med, rfungi = rf_med, ba_ar = baar_med,
                    pl_vi_unk = pl_med, endophytes = endo_med)




# hatmap with endophytes
pheatmap(heat_group, 
         scale = "none",
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
         scale = "none",
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
raw_raw_gfract <- read.table("assmble/gsa/raw_assembly/summary/TSV/raw_Genome_fraction.tsv",
                        sep = "\t", header = TRUE, row.names = 1,
                        na.strings = "-")
# Replace NAs
raw_raw_gfract[is.na(raw_raw_gfract)] <- 0

# Clean up data and add columns
raw_raw_gfract <- clean_data(raw_raw_gfract, "raw")


# Index out the groups
r_endo <- raw_raw_gfract[raw_raw_gfract$group != "host",1:11]
r_host <- raw_raw_gfract[raw_raw_gfract$group == "host",1:11]
r_mock <- raw_raw_gfract[,1:11]
r_rf <- raw_raw_gfract[raw_raw_gfract$group == "rfungi",1:11]
r_omf <- raw_raw_gfract[raw_raw_gfract$group == "OMF",1:11]
r_ba_ar <- raw_raw_gfract[raw_raw_gfract$group == "ba_ar",1:11]
r_pl_vi_unk <- raw_raw_gfract[raw_raw_gfract$group == "pl_vi_unk",1:11]


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



