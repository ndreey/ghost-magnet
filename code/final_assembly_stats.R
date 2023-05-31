library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
source("analysis_tools.R")

# Load in data
# GSA
#nga50 <- read.table("assmble/gsa/raw_assembly",
#                    sep = "\t", header = TRUE, row.names = 1,
#                    na.strings = "-")
#

# raw assembly
raw_nga50 <- read.table("assmble/raw/summary/TSV/num_contigs.tsv",
                        sep = "\t", header = TRUE, row.names = 1,
                        na.strings = "-")

# 1000c assembly
x1000c_nga50 <- read.table("assmble/x1000c_quast/num_contigs.tsv",
                           sep = "", header = TRUE, row.names = 1,
                           na.strings = "-")
# 700c assembly
x700c_nga50 <- read.table("assmble/x700c_quast/num_contigs.tsv",
                          sep = "", header = TRUE, row.names = 1,
                          na.strings = "-")

# maxbin assembly
maxbin_nga50 <- read.table("assmble/maxbin_quast/num_contigs.tsv",
                           sep = "", header = TRUE, row.names = 1,
                           na.strings = "-")

not_aln <- data.frame(raw = as.vector(raw_nga50[101,]), 
                 x1000c = as.vector(x1000c_nga50[101,]), 
                 x700c = as.vector(x700c_nga50[101,]),
                 maxbin = as.vector(maxbin_nga50[101,]))

raw_nga50 <- raw_nga50[-101,]
x1000c_nga50 <- x1000c_nga50[-101,]
x700c_nga50 <- x700c_nga50[-101,]
maxbin_nga50 <- maxbin_nga50[-101,]

# Quick check
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
#nga50 <- clean_data(nga50, "gsa")
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






library(reshape2)

# Melt the dataframe to convert it into long format
melted_data <- melt(raw_nga50, id.vars = c("group", "repli"), variable.name = "sample", value.name = "values")
me2 <- melt(x1000c_nga50, id.vars = c("group", "repli"), variable.name = "sample", value.name = "values")
me3 <- melt(x700c_nga50, id.vars = c("group", "repli"), variable.name = "sample", value.name = "values")
me4<- melt(maxbin_nga50, id.vars = c("group", "repli"), variable.name = "sample", value.name = "values")

counts_long <- rbind(melted_data, me2, me3, me4)




ggplot(counts_long, aes(x=sample, y = values, fill = repli)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_point(aes(color = group))+
  theme_bw() +
  labs(title = "Num of contigs", y = "Counts") +
  ylim(0,6000)




par(mfrow = c(2, 2))

boxplot(raw_endo,
        main = "Endophyte contigs\nraw",
        ylab = "Counts",
        xlab = "Sample",
        xaxt = "n",
        cex = 1)

boxplot(raw_mock,
        main = "Mock Community contigs\nraw",
        ylab = "Counts",
        xlab = "S",
        xaxt = "n",
        cex = 1)



boxplot(x1000c_endo,
        main = "Endophyte contigs\nx1000c",
        ylab = "",
        xlab = "",
        xaxt = "n",
        yaxt = "n",
        cex = 1)

boxplot(x1000c_mock,
        main = "Mock Community contigs\nx1000c",
        ylab = "",
        xlab = "",
        xaxt = "n",
        yaxt = "n",
        cex = 1)





boxplot(x700c_endo,
        main = "Endophyte contigs\nx700c",
        ylab = "Counts",
        xlab = "Sample",
        cex = 1)

boxplot(x700c_mock,
        main = "Mock Community contigs\nx700c",
        ylab = "Counts",
        xlab = "Sample",
        cex = 1)





boxplot(maxbin_endo,
        main = "Endophyte contigs\nmaxbin",
        ylab = "",
        xlab = "Sample",
      yxaxt = "n",
        cex = 1)

boxplot(maxbin_mock,
        main = "Mock Community contigs\nmaxbin",
        ylab = "",
        xlab = "Sample",
        yaxt = "n",
        cex = 1)


par(mfrow = c(4, 2))

# Endophytes
boxplot(raw_endo, main = "Endophyte contigs\nraw", ylab = "Counts")
boxplot(x1000c_endo, main = "raw Endophytes")
boxplot(x700c_endo, main = "raw Endophytes")
boxplot(maxbin_endo, main = "Endophytes") 

# All
#boxplot(mock, main = "mock data")
boxplot(raw_mock, main = "raw")
boxplot(x1000c_mock,  main = "x1000c")
boxplot(x700c_mock,  main = "x700c")
boxplot(maxbin_mock,  main = "maxbin")


# Endophytes
boxplot(raw_host, main = "raw Endophytes")
boxplot(x1000c_host, main = "raw Endophytes")
boxplot(x700c_host, main = "raw Endophytes")
boxplot(maxbin_host, main = "Endophytes") 




# rfungi
boxplot(raw_rf, main = "raw rfungi")
boxplot(x1000c_rf, main = "x1000c rfungi")
boxplot(x700c_rf, main = "x700c rfungi")
boxplot(maxbin_rf, main = "maxbin rfungi")

# OMF
boxplot(raw_omf, main = "raw OMF")
boxplot(x1000c_omf, main = "x100c OMF")
boxplot(x700c_omf, main = "x700c OMF")
boxplot(maxbin_omf, main = "maxbin OMF")

# Fungi
boxplot(raw_fungi, main = "raw Fungi")
boxplot(x1000c_fungi, main = "x100c Fungi")
boxplot(x700c_fungi, main = "x700c Fungi")
boxplot(maxbin_fungi, main = "maxbin Fungi")

# ba_ar
boxplot(ba_ar, main = "ba_ar")
boxplot(raw_ba_ar, main = "raw ba_ar")

# pl_vi_unk
boxplot(pl_vi_unk, main = "pl_vi_unk")
boxplot(raw_pl_vi_unk, main = "raw pl_vi_unk")

par(mfrow = c(1, 1))

# RAW
raw_median <- apply(raw_mock, 2, median, na.rm = TRUE)
shapiro.test(raw_median) 
boxplot(raw_median)
hist(raw_median)


maxbin_median <- apply(maxbin_mock, 2, median, na.rm = TRUE)
shapiro.test(maxbin_median) 
boxplot(maxbin_median)
hist(maxbin_median)
