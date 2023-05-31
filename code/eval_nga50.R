library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggpubr)
source("analysis_tools.R")

# Load in data
# GSA
nga50 <- read.table("assmble/gsa/gold/TSV/gsa_NGA50.tsv",
                   sep = "\t", header = TRUE, row.names = 1,
                   na.strings = "-")

# raw assembly
raw_nga50 <- read.table("assmble/gsa/raw_assembly/summary/TSV/raw_NGA50.tsv",
                       sep = "\t", header = TRUE, row.names = 1,
                       na.strings = "-")

# 1000c assembly
x1000c_nga50 <- read.table("assmble/gsa/1000cc_filtered/TSV/NGA50.tsv",
                          sep = "\t", header = TRUE, row.names = 1,
                          na.strings = "-")
# 700c assembly
x700c_nga50 <- read.table("assmble/gsa/700cc_filtered/TSV/NGA50.tsv",
                         sep = "\t", header = TRUE, row.names = 1,
                         na.strings = "-")

# maxbin assembly
maxbin_nga50 <- read.table("assmble/gsa/maxbin/TSV/NGA50.tsv",
                          sep = "\t", header = TRUE, row.names = 1,
                          na.strings = "-")
# Quick check
dim(nga50) 
dim(raw_nga50)
dim(x1000c_nga50)
dim(x700c_nga50)
dim(maxbin_nga50)

# Addin missing samples to filtered.
x1000c_nga50$X01 <- raw_nga50$X01_sensi.contigs
x1000c_nga50$X00 <- raw_nga50$X00_sensi.contigs
x1000c_nga50 <- x1000c_nga50[,c(11,10,1:9)]

x700c_nga50$X00 <- raw_nga50$X00_sensi.contigs
x700c_nga50 <- x700c_nga50[,c(11,1:10)]

maxbin_nga50$X00 <- raw_nga50$X00_sensi.contigs
maxbin_nga50 <- maxbin_nga50[,c(11,1:10)]


# Clean up data and add columns
nga50 <- clean_data(nga50, "gsa")
raw_nga50 <- clean_data(raw_nga50, "raw")
x1000c_nga50 <- clean_data(x1000c_nga50, "x1000c")
x700c_nga50 <- clean_data(x700c_nga50, "x700c")
maxbin_nga50 <- clean_data(maxbin_nga50, "maxbin")


# STANDARD VECTORS
# HC percentages
hc <- c(0, 0.254, 0.429, 0.566, 0.664, 0.759, 0.818, 0.878, 0.923, 0.964, 0.983)
# Sample names
sample_id <- c("X00", "X01", "X02", "X03", "X04", "X05", "X06", "X07",
               "X08", "X090", "X095")


# Index out the groups
endo <- nga50[nga50$group != "host",1:11]
host <- nga50[nga50$group == "host",1:11]
mock <- nga50[,1:11]
rf <- nga50[nga50$group == "rfungi",1:11]
omf <- nga50[nga50$group == "OMF",1:11]
fungi <- nga50[nga50$group %in% c("rfungi", "OMF"), 1:11]
ba_ar <- nga50[nga50$group == "ba_ar",1:11]
pl_vi_unk <- nga50[nga50$group == "pl_vi_unk",1:11]

raw_endo <- raw_nga50[raw_nga50$group != "host",1:11]
raw_host <- raw_nga50[raw_nga50$group == "host",1:11]
raw_mock <- raw_nga50[,1:11]
raw_rf <- raw_nga50[raw_nga50$group == "rfungi",1:11]
raw_omf <- raw_nga50[raw_nga50$group == "OMF",1:11]
raw_fungi <- raw_nga50[raw_nga50$group %in% c("rfungi", "OMF"), 1:11]
raw_ba_ar <- raw_nga50[raw_nga50$group == "ba_ar",1:11]
raw_pl_vi_unk <- raw_nga50[raw_nga50$group == "pl_vi_unk",1:11]

x1000c_endo <- x1000c_nga50[x1000c_nga50$group != "host",1:11]
x1000c_host <- x1000c_nga50[x1000c_nga50$group == "host",1:11]
x1000c_mock <- x1000c_nga50[,1:11]
x1000c_rf <- x1000c_nga50[x1000c_nga50$group == "rfungi",1:11]
x1000c_omf <- x1000c_nga50[x1000c_nga50$group == "OMF",1:11]
x1000c_fungi <- x1000c_nga50[x1000c_nga50$group %in% c("rfungi", "OMF"), 1:11]
x1000c_ba_ar <- x1000c_nga50[x1000c_nga50$group == "ba_ar",1:11]
x1000c_pl_vi_unk <- x1000c_nga50[x1000c_nga50$group == "pl_vi_unk",1:11]

x700c_endo <- x700c_nga50[x700c_nga50$group != "host",1:11]
x700c_host <- x700c_nga50[x700c_nga50$group == "host",1:11]
x700c_mock <- x700c_nga50[,1:11]
x700c_rf <- x700c_nga50[x700c_nga50$group == "rfungi",1:11]
x700c_omf <- x700c_nga50[x700c_nga50$group == "OMF",1:11]
x700c_fungi <- x700c_nga50[x700c_nga50$group %in% c("rfungi", "OMF"), 1:11]
x700c_ba_ar <- x700c_nga50[x700c_nga50$group == "ba_ar",1:11]
x700c_pl_vi_unk <- x700c_nga50[x700c_nga50$group == "pl_vi_unk",1:11]

maxbin_endo <- maxbin_nga50[maxbin_nga50$group != "host",1:11]
maxbin_host <- maxbin_nga50[maxbin_nga50$group == "host",1:11]
maxbin_mock <- maxbin_nga50[,1:11]
maxbin_rf <- maxbin_nga50[maxbin_nga50$group == "rfungi",1:11]
maxbin_omf <- maxbin_nga50[maxbin_nga50$group == "OMF",1:11]
maxbin_fungi <- maxbin_nga50[maxbin_nga50$group %in% c("rfungi", "OMF"), 1:11]
maxbin_ba_ar <- maxbin_nga50[maxbin_nga50$group == "ba_ar",1:11]
maxbin_pl_vi_unk <- maxbin_nga50[maxbin_nga50$group == "pl_vi_unk",1:11]


par(mfrow = c(1, 2))

# Endophytes
boxplot(endo, main = "Endophytes") 
boxplot(raw_endo, main = "raw Endophytes")

# All
boxplot(mock, main = "mock data")
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
boxplot(fungi, main = "Fungi")
boxplot(raw_fungi, main = "raw Fungi")

# ba_ar
boxplot(ba_ar, main = "ba_ar")
boxplot(raw_ba_ar, main = "raw ba_ar")

# pl_vi_unk
boxplot(pl_vi_unk, main = "pl_vi_unk")
boxplot(raw_pl_vi_unk, main = "raw pl_vi_unk")

par(mfrow = c(1, 1))



# Analysis
# GSA
median <- apply(mock, 2, median, na.rm = TRUE)
shapiro.test(median) # W = 0.94097, p-value = 0.6206 norm dist
boxplot(median)
hist(median)

# Sig. neg. cor.
cor.test(hc, median, method ="pearson") 
# p-value = 0.0001428   r = -0.9611134


# RAW
raw_median <- apply(raw_mock, 2, median, na.rm = TRUE)
shapiro.test(raw_median) 
boxplot(raw_median)
hist(raw_median)

# 
cor.test(hc, raw_median, method ="kendall") 
#



# Stats
median <- apply(mock, 2, median, na.rm = TRUE)
variance <- apply(mock, 2, var, na.rm = TRUE)
sd <- apply(mock, 2, sd, na.rm = TRUE)

Q1 <- vector()
Q3 <- vector()
IQR <- vector()


for (sample in 1:11) {
  print(sample)
  q1_temp <- quantile(mock[,sample], 0.25, na.rm = TRUE)
  print(q1_temp)
  q2_temp <- quantile(mock[,sample], 0.75, na.rm = TRUE)
  iqr_temp <- q2_temp - q1_temp
  
  Q1 <- c(Q1,q1_temp)
  Q3 <- c(Q3,q2_temp)
  IQR <- c(IQR,iqr_temp)
  
}

stats <- data.frame(median = median, variance = variance, sd = sd,
                    sample = sample_id, Q1 = Q1, Q3 = Q3, iqr = IQR)

plot(hc, median, main = "mock df\nMedian", xlab = "Host Contamination (%)")
plot(hc, variance, main = "mock df\nVariance")


# Boxplot
boxplot(mock[1:8], names = sample_id[1:8], xlab = "Sample", ylim=c(0,12000),
        ylab = "NGA50 (bp)", main = "NGA50 Disribution\nof\nGSA Samples")


# There were only NA after X07
# Convert to long-format
raw_long <- pivot_longer(data = raw_nga50[c(1:8, 12, 13)], cols = colnames(raw_nga50)[1:8],
                         names_to = "sample", values_to="value")

gsa_long <- pivot_longer(data = nga50[c(1:8, 12, 13)], cols = colnames(nga50)[1:8],
                         names_to = "sample", values_to="value")

# View the first few rows
head(raw_long)
head(gsa_long)

merged_long <- rbind(raw_long, gsa_long)

# Boxplot
ggplot(merged_long, aes(x=sample, y = value, fill = repli)) +
  geom_boxplot(alpha=0.4) +
  theme_light() +
  labs(title = "NGA50", y = "NGA50 (bp)")





