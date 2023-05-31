library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggpubr)
source("analysis_tools.R")



# Load in data
# GSA
dupe <- read.table("assmble/gsa/gold/TSV/gsa_Duplication_ratio.tsv",
                        sep = "\t", header = TRUE, row.names = 1,
                        na.strings = "-")

# raw assembly
raw_dupe <- read.table("assmble/gsa/raw_assembly/summary/TSV/raw_Duplication_ratio.tsv",
                         sep = "\t", header = TRUE, row.names = 1,
                         na.strings = "-")

# 1000c assembly
x1000c_dupe <- read.table("assmble/gsa/1000cc_filtered/TSV/Duplication_ratio.tsv",
                            sep = "\t", header = TRUE, row.names = 1,
                            na.strings = "-")
# 700c assembly
x700c_dupe <- read.table("assmble/gsa/700cc_filtered/TSV/Duplication_ratio.tsv",
                           sep = "\t", header = TRUE, row.names = 1,
                           na.strings = "-")

# maxbin assembly
maxbin_dupe <- read.table("assmble/gsa/maxbin/TSV/Duplication_ratio.tsv",
                            sep = "\t", header = TRUE, row.names = 1,
                            na.strings = "-")

# Quick check
dim(dupe) 
dim(raw_dupe)
dim(x1000c_dupe)
dim(x700c_dupe)
dim(maxbin_dupe)

# Addin missing samples to filtered.
x1000c_dupe$X01 <- raw_dupe$X01_sensi.contigs
x1000c_dupe$X00 <- raw_dupe$X00_sensi.contigs
x1000c_dupe <- x1000c_dupe[,c(11,10,1:9)]

x700c_dupe$X00 <- raw_dupe$X00_sensi.contigs
x700c_dupe <- x700c_dupe[,c(11,1:10)]

maxbin_dupe$X00 <- raw_dupe$X00_sensi.contigs
maxbin_dupe <- maxbin_dupe[,c(11,1:10)]


# Clean up data and add columns
dupe <- clean_data(dupe, "gsa")
raw_dupe <- clean_data(raw_dupe, "raw")
x1000c_dupe <- clean_data(x1000c_dupe, "x1000c")
x700c_dupe <- clean_data(x700c_dupe, "x700c")
maxbin_dupe <- clean_data(maxbin_dupe, "maxbin")


# STANDARD VECTORS
# HC percentages
hc <- c(0, 0.254, 0.429, 0.566, 0.664, 0.759, 0.818, 0.878, 0.923, 0.964, 0.983)
# Sample names
sample_id <- c("X00", "X01", "X02", "X03", "X04", "X05", "X06", "X07",
               "X08", "X090", "X095")


# Index out the groups
endo <- dupe[dupe$group != "host",1:11]
host <- dupe[dupe$group == "host",1:11]
mock <- dupe[,1:11]
rf <- dupe[dupe$group == "rfungi",1:11]
omf <- dupe[dupe$group == "OMF",1:11]
fungi <- dupe[dupe$group %in% c("rfungi", "OMF"), 1:11]
ba_ar <- dupe[dupe$group == "ba_ar",1:11]
pl_vi_unk <- dupe[dupe$group == "pl_vi_unk",1:11]

raw_endo <- raw_dupe[raw_dupe$group != "host",1:11]
raw_host <- raw_dupe[raw_dupe$group == "host",1:11]
raw_mock <- raw_dupe[,1:11]
raw_rf <- raw_dupe[raw_dupe$group == "rfungi",1:11]
raw_omf <- raw_dupe[raw_dupe$group == "OMF",1:11]
raw_fungi <- raw_dupe[raw_dupe$group %in% c("rfungi", "OMF"), 1:11]
raw_ba_ar <- raw_dupe[raw_dupe$group == "ba_ar",1:11]
raw_pl_vi_unk <- raw_dupe[raw_dupe$group == "pl_vi_unk",1:11]

x1000c_endo <- x1000c_dupe[x1000c_dupe$group != "host",1:11]
x1000c_host <- x1000c_dupe[x1000c_dupe$group == "host",1:11]
x1000c_mock <- x1000c_dupe[,1:11]
x1000c_rf <- x1000c_dupe[x1000c_dupe$group == "rfungi",1:11]
x1000c_omf <- x1000c_dupe[x1000c_dupe$group == "OMF",1:11]
x1000c_fungi <- x1000c_dupe[x1000c_dupe$group %in% c("rfungi", "OMF"), 1:11]
x1000c_ba_ar <- x1000c_dupe[x1000c_dupe$group == "ba_ar",1:11]
x1000c_pl_vi_unk <- x1000c_dupe[x1000c_dupe$group == "pl_vi_unk",1:11]

x700c_endo <- x700c_dupe[x700c_dupe$group != "host",1:11]
x700c_host <- x700c_dupe[x700c_dupe$group == "host",1:11]
x700c_mock <- x700c_dupe[,1:11]
x700c_rf <- x700c_dupe[x700c_dupe$group == "rfungi",1:11]
x700c_omf <- x700c_dupe[x700c_dupe$group == "OMF",1:11]
x700c_fungi <- x700c_dupe[x700c_dupe$group %in% c("rfungi", "OMF"), 1:11]
x700c_ba_ar <- x700c_dupe[x700c_dupe$group == "ba_ar",1:11]
x700c_pl_vi_unk <- x700c_dupe[x700c_dupe$group == "pl_vi_unk",1:11]

maxbin_endo <- maxbin_dupe[maxbin_dupe$group != "host",1:11]
maxbin_host <- maxbin_dupe[maxbin_dupe$group == "host",1:11]
maxbin_mock <- maxbin_dupe[,1:11]
maxbin_rf <- maxbin_dupe[maxbin_dupe$group == "rfungi",1:11]
maxbin_omf <- maxbin_dupe[maxbin_dupe$group == "OMF",1:11]
maxbin_fungi <- maxbin_dupe[maxbin_dupe$group %in% c("rfungi", "OMF"), 1:11]
maxbin_ba_ar <- maxbin_dupe[maxbin_dupe$group == "ba_ar",1:11]
maxbin_pl_vi_unk <- maxbin_dupe[maxbin_dupe$group == "pl_vi_unk",1:11]


par(mfrow = c(1, 3))

# Endophytes
boxplot(endo, main = "Endophytes") 
boxplot(raw_endo, main = "raw Endophytes")

# All
boxplot(mock, main = "mock community")
boxplot(dupe[dupe$group !="pl_vi_unk",1:11], main = "pl_vi_unk excluded")

boxplot(raw_mock, main = "raw mock data")
boxplot(x1000c_mock,  main = "x1000c")
boxplot(x700c_mock,  main = "x700c")
boxplot(maxbin_mock,  main = "maxbin")

# rfungi
boxplot(rf, main = "rfungi")
boxplot(raw_rf, main = "raw rfungi")

# OMF
boxplot(omf, main = "OMF")
boxplot(raw_omf, main = "raw OMF")

# Fungi
boxplot(fungi, main = "fungi")
boxplot(raw_fungi, main = "raw Fungi")

# ba_ar
boxplot(ba_ar, main = "ba_ar")
boxplot(raw_ba_ar, main = "raw ba_ar")

# pl_vi_unk
boxplot(pl_vi_unk, main = "pl_vi_unk")
boxplot(raw_pl_vi_unk, main = "raw pl_vi_unk")

par(mfrow = c(1, 1))


################## Analysis #######################
# GSA
median <- apply(mock, 2, median, na.rm = TRUE)
sd <- apply(mock, 2, sd, na.rm = TRUE)


mock_no_pl <- dupe[dupe$group != "pl_vi_unk",1:11]

median2 <- apply(mock_no_pl, 2, median, na.rm = TRUE)
sd2 <- apply(mock_no_pl, 2, sd, na.rm = TRUE)

# shapiro.test(median) # Error in shapiro.test(median) : all 'x' values are identical
# X00  X01  X02  X03  X04  X05  X06  X07  X08 X090 X095 
#  1    1    1    1    1    1    1    1    1    1    1 

# RAW
raw_median <- apply(raw_mock, 2, median, na.rm = TRUE)
shapiro.test(raw_median)  # W = 0.46028, p-value = 5.065e-07 Not norm

# Correlation
cor.test(hc, raw_median, method ="kendall", exact = FALSE) 
#p-value = 0.02616 r = 0.5877538 

# x1000c
x1000c_median <- apply(x1000c_mock, 2, median, na.rm = TRUE)
shapiro.test(x1000c_median)  # W = 0.48861, p-value = 1.106e-06 Not norm

# Correlation
cor.test(hc, x1000c_median, method ="kendall", exact = FALSE) 
#p-value = 0.02616 r = 0.5877538 

# x700c
x700c_median <- apply(x700c_mock, 2, median, na.rm = TRUE)
shapiro.test(x700c_median)  # W = 0.34499, p-value = 2.243e-08 Not norm

# Correlation
cor.test(hc, x700c_median, method ="kendall", exact = FALSE) 
#p-value = 0.2059 r = 0.3411211

# Maxbin
maxbin_median <- apply(maxbin_mock, 2, median, na.rm = TRUE)
# shapiro.test(maxbin_median)  # All values are the same
# X00  X01  X02  X03  X04  X05  X06  X07  X08 X090 X095 
#  1    1    1    1    1    1    1    1    1    1    1


plot_maxbin_stats(mock, "gsa")
plot_maxbin_stats(raw_mock, "raw")
plot_maxbin_stats(x1000c_mock, "x1000c")
plot_maxbin_stats(x700c_mock, "x700c")
plot_maxbin_stats(maxbin_mock, "maxbin")




# Convert to long-format
gsa_long <- pivot_longer(data = dupe, cols = colnames(dupe)[1:11],
                         names_to = "sample", values_to="value")

raw_long <- pivot_longer(data = raw_dupe, cols = colnames(raw_dupe)[1:11],
                         names_to = "sample", values_to="value")

x1000c_long <- pivot_longer(data = x1000c_dupe, cols = colnames(x1000c_dupe)[1:11],
                            names_to = "sample", values_to="value")

x700c_long <- pivot_longer(data = x700c_dupe, cols = colnames(x700c_dupe)[1:11],
                           names_to = "sample", values_to="value")

maxbin_long <- pivot_longer(data = maxbin_dupe, cols = colnames(maxbin_dupe)[1:11],
                            names_to = "sample", values_to="value")

# View the first few rows
head(raw_long)
head(gsa_long)
head(x1000c_long)
head(x700c_long)

merged_long <- rbind(gsa_long, raw_long, x1000c_long, x700c_long, maxbin_long)
merged_long$repli <- factor(merged_long$repli, levels = c("gsa", "raw", "x1000c",
                                                          "x700c", "maxbin"))
merged_long$hc <- hc[match(merged_long$sample, sample_id)]


# Boxplot
dupe_p1 <- ggplot(merged_long, aes(x=sample, y = value, fill = repli)) +
  geom_boxplot(alpha=0.4, position = "dodge") +
  theme_bw() +
  theme(axis.text.x=element_blank(),legend.position = "top") +
  labs(y = "Ratio", x = "")

# Boxplot
dupe_p2 <- ggplot(merged_long, aes(x=sample, y = value, fill = repli)) +
  geom_boxplot(alpha=0.4, position = "dodge") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(subtitle = "ylim: (0.995,1.02)", y = "Ratio") +
  ylim(0.995,1.02)

ggarrange(dupe_p1, dupe_p2, labels = c("A", "B"), ncol = 1, nrow = 2)

