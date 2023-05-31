library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
# Long table cols
num <- vector()
prcnt <- vector()
repli <- vector()
c_ln <- vector()

####################### 00raw ##############################
raw_00 <- read.table("contigs/raw/00_sensi.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)

x00_all <- raw_00$V2


# Mean and median
mean(x00_all)        # 691.9308
median(x00_all)      # 450

# Range of lengths
range(x00_all)       # 200 16053

par(mfrow = c(1, 2))

# Distribution
hist(log(x00_all), xlab = "Length of contig (bp)\nLogscaled",
     main = "x00\nContig length distribution")

boxplot(x00_all, ylim = c(0,900), main = "Boxplot x00\nylim: (0, 900)")



# Number
for (c_len in c(500, 700, 1000, 1500)) {
  
  # Number and percentage of contigs equal or less than c_len
  x1 = length(x00_all[x00_all >= c_len])
  x2 = length(x00_all[x00_all >= c_len]) / length(x00_all)
  print(paste("Num contigs >=", c_len, ":    ", x1, sep = ""))
  print(paste("Percentage:           ", round(x2*100, 2), "%", sep = ""))
  print("")
  num <- c(num, x1)
  prcnt <- c(prcnt, x2)
  repli <- c(repli,"x00")
  c_ln <- c(c_ln, c_len)
  
}




####################### raw ##############################

raw_01 <- read.table("contigs/raw/01_sensi.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
raw_02 <- read.table("contigs/raw/02_sensi.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
raw_03 <- read.table("contigs/raw/03_sensi.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
raw_04 <- read.table("contigs/raw/04_sensi.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
raw_05 <- read.table("contigs/raw/05_sensi.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
raw_06 <- read.table("contigs/raw/06_sensi.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
raw_07 <- read.table("contigs/raw/07_sensi.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
raw_08 <- read.table("contigs/raw/08_sensi.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
raw_090 <- read.table("contigs/raw/090_sensi.contigs.fa_stats.tsv", sep = "\t",
                      header = FALSE)
raw_095 <- read.table("contigs/raw/095_sensi.contigs.fa_stats.tsv", sep = "\t",
                      header = FALSE)

raw_all <- c(raw_01$V2, raw_02$V2, raw_03$V2, raw_04$V2, raw_05$V2,
             raw_06$V2, raw_07$V2, raw_08$V2, raw_090$V2, raw_095$V2)


# Mean and median
mean(raw_all)        # 497.3246
median(raw_all)      # 412

# Range of lengths
range(raw_all)       # 200 16703

par(mfrow = c(1, 2))

# Distribution
hist(log(raw_all), xlab = "Length of contig (bp)\nLogscaled",
     main = "raw\nContig length distribution")

boxplot(raw_all, ylim = c(0,900), main = "Boxplot raw\nylim: (0, 900)")



# Number
for (c_len in c(500, 700, 1000, 1500)) {
  
  # Number and percentage of contigs equal or less than c_len
  x1 = length(raw_all[raw_all >= c_len])
  x2 = length(raw_all[raw_all >= c_len]) / length(raw_all)
  print(paste("Num contigs >=", c_len, ":    ", x1, sep = ""))
  print(paste("Percentage:           ", round(x2*100, 2), "%", sep = ""))
  print("")
  
  num <- c(num, x1)
  prcnt <- c(prcnt, x2)
  repli <- c(repli,"raw")
  c_ln <- c(c_ln, c_len)
}




####################### x1000c ##############################
x1000c_02 <- read.table("contigs/x1000c/02_1000c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x1000c_03 <- read.table("contigs/x1000c/03_1000c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x1000c_04 <- read.table("contigs/x1000c/04_1000c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x1000c_05 <- read.table("contigs/x1000c/05_1000c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x1000c_06 <- read.table("contigs/x1000c/06_1000c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x1000c_07 <- read.table("contigs/x1000c/07_1000c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x1000c_08 <- read.table("contigs/x1000c/08_1000c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x1000c_090 <- read.table("contigs/x1000c/090_1000c.contigs.fa_stats.tsv", sep = "\t",
                         header = FALSE)
x1000c_095 <- read.table("contigs/x1000c/095_1000c.contigs.fa_stats.tsv", sep = "\t",
                         header = FALSE)

x1000c_all <- c(x1000c_02$V2, x1000c_03$V2, x1000c_04$V2, x1000c_05$V2,
                x1000c_06$V2, x1000c_07$V2, x1000c_08$V2, x1000c_090$V2, x1000c_095$V2)


# Mean and median
mean(x1000c_all)        # 480.9691
median(x1000c_all)      # 399

# Range of lengths
range(x1000c_all)       # 200 16703

par(mfrow = c(1, 2))

# Distribution
hist(log(x1000c_all), xlab = "Length of contig (bp)\nLogscaled",
     main = "x1000c\nContig length distribution")

boxplot(x1000c_all, ylim = c(0,900), main = "Boxplot x1000c\nylim: (0, 900)")



# Number
for (c_len in c(500, 700, 1000, 1500)) {
  
  # Number and percentage of contigs equal or less than c_len
  x1 = length(x1000c_all[x1000c_all >= c_len])
  x2 = length(x1000c_all[x1000c_all >= c_len]) / length(x1000c_all)
  print(paste("Num contigs >=", c_len, ":    ", x1, sep = ""))
  print(paste("Percentage:           ", round(x2*100, 2), "%", sep = ""))
  print("")
  
  num <- c(num, x1)
  prcnt <- c(prcnt, x2)
  repli <- c(repli,"x1000c")
  c_ln <- c(c_ln, c_len)
}




####################### x700c ##############################
x700c_01 <- read.table("contigs/x700c/01_700c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x700c_02 <- read.table("contigs/x700c/02_700c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x700c_03 <- read.table("contigs/x700c/03_700c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x700c_04 <- read.table("contigs/x700c/04_700c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x700c_05 <- read.table("contigs/x700c/05_700c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x700c_06 <- read.table("contigs/x700c/06_700c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x700c_07 <- read.table("contigs/x700c/07_700c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x700c_08 <- read.table("contigs/x700c/08_700c.contigs.fa_stats.tsv", sep = "\t",
                        header = FALSE)
x700c_090 <- read.table("contigs/x700c/090_700c.contigs.fa_stats.tsv", sep = "\t",
                         header = FALSE)
x700c_095 <- read.table("contigs/x700c/095_700c.contigs.fa_stats.tsv", sep = "\t",
                         header = FALSE)

x700c_all <- c(x700c_01$V2, x700c_02$V2, x700c_03$V2, x700c_04$V2, x700c_05$V2,
                x700c_06$V2, x700c_07$V2, x700c_08$V2, x700c_090$V2, x700c_095$V2)


# Mean and median
mean(x700c_all)        # 462.6866
median(x700c_all)      # 396

# Range of lengths
range(x700c_all)       # 200 16703

par(mfrow = c(1, 2))

# Distribution
hist(log(x700c_all), xlab = "Length of contig (bp)\nLogscaled",
     main = "x700c\nContig length distribution")

boxplot(x700c_all, ylim = c(0,900), main = "Boxplot x700c\nylim: (0, 900)")


# Number
for (c_len in c(500, 700, 1000, 1500)) {
  
  # Number and percentage of contigs equal or less than c_len
  x1 = length(x700c_all[x700c_all >= c_len])
  x2 = length(x700c_all[x700c_all >= c_len]) / length(x700c_all)
  print(paste("Num contigs >=", c_len, ":    ", x1, sep = ""))
  print(paste("Percentage:           ", round(x2*100, 2), "%", sep = ""))
  print("")
  
  num <- c(num, x1)
  prcnt <- c(prcnt, x2)
  repli <- c(repli,"x700c")
  c_ln <- c(c_ln, c_len)
}




####################### Maxbin ##############################
max_01 <- read.table("contigs/max_bin/01_max_bin.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
max_02 <- read.table("contigs/max_bin/02_max_bin.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
max_03 <- read.table("contigs/max_bin/03_max_bin.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
max_04 <- read.table("contigs/max_bin/04_max_bin.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
max_05 <- read.table("contigs/max_bin/05_max_bin.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
max_06 <- read.table("contigs/max_bin/06_max_bin.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
max_07 <- read.table("contigs/max_bin/07_max_bin.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
max_08 <- read.table("contigs/max_bin/08_max_bin.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
max_090 <- read.table("contigs/max_bin/090_max_bin.contigs.fa_stats.tsv", sep = "\t",
                      header = FALSE)
max_095 <- read.table("contigs/max_bin/095_max_bin.contigs.fa_stats.tsv", sep = "\t",
                      header = FALSE)

max_all <- c(max_01$V2, max_02$V2, max_03$V2, max_04$V2, max_05$V2,
             max_06$V2, max_07$V2, max_08$V2, max_090$V2, max_095$V2)


# Mean and median
mean(max_all)        # 468.600
median(max_all)      # 402

# Range of lengths
range(max_all)       # 200 16703

par(mfrow = c(1, 2))

# Distribution
hist(log(max_all), xlab = "Length of contig (bp)\nLogscaled",
     main = "max_bin\nContig length distribution")

boxplot(max_all, ylim = c(0,900), main = "Boxplot max_bin\nylim: (0, 900)")


# Number
for (c_len in c(200,500, 700, 1000, 1500)) {
  
  # Number and percentage of contigs equal or less than c_len
  x1 = length(max_095$V2[max_095$V2 >= c_len])
  x2 = length(max_095$V2[max_095$V2 >= c_len]) / length(max_all)
  print(paste("Num contigs >=", c_len, ":    ", x1, sep = ""))
  print(paste("Percentage:           ", round(x2*100, 2), "%", sep = ""))
  print("")
  
  #num <- c(num, x1)
  #prcnt <- c(prcnt, x2)
  #repli <- c(repli,"max")
  #c_ln <- c(c_ln, c_len)
}

#c(500, 700, 1000, 1500)
#c("1500", "1000", "700", "500")

###### PLOT #######

repli <- factor(repli, levels = c("x00", "raw", "x1000c", "x700c", "max"))
c_ln <- as.character(c_ln)
c_ln <- factor(c_ln, levels = c("500", "700", "1000", "1500"))

df1_long <- data.frame(prcnt, num, repli, c_ln)

# Percentage
ggplot(df1_long, aes(x = c_ln, y = prcnt, fill = repli)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", alpha = 0.6) +
  theme_light() +
  scale_fill_brewer(palette = "Pastel1") +
  labs(title = "Percentage of all reads that are equal or less than x", 
       y = "Fraction", x = "Contig lengths (bp)")




# Counts
ggplot(df1_long, aes(x = c_ln, y = num, fill = repli)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  theme_light() +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Number of reads that are equal or less than x", 
       y = "Count", x = "Contig lengths (bp)")

################ PLAYGROUND ######################

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggpubr)
source("analysis_tools.R")

# Load in data
# GSA
gfrac <- read.table("assmble/gsa/gold/TSV/gsa_Genome_fraction.tsv",
                    sep = "\t", header = TRUE, row.names = 1,
                    na.strings = "-")

# raw assembly
raw_gfrac <- read.table("assmble/gsa/raw_assembly/summary/TSV/raw_Genome_fraction.tsv",
                        sep = "\t", header = TRUE, row.names = 1,
                        na.strings = "-")

# 1000c assembly
x1000c_gfrac <- read.table("assmble/gsa/1000cc_filtered/TSV/Genome_fraction.tsv",
                           sep = "\t", header = TRUE, row.names = 1,
                           na.strings = "-")
# 700c assembly
x700c_gfrac <- read.table("assmble/gsa/700cc_filtered/TSV/Genome_fraction.tsv",
                          sep = "\t", header = TRUE, row.names = 1,
                          na.strings = "-")

# maxbin assembly
maxbin_gfrac <- read.table("assmble/gsa/maxbin/TSV/Genome_fraction.tsv",
                           sep = "\t", header = TRUE, row.names = 1,
                           na.strings = "-")
# Quick check
dim(gfrac) 
dim(raw_gfrac)
dim(x1000c_gfrac)
dim(x700c_gfrac)
dim(maxbin_gfrac)

# Addin missing samples to filtered.
x1000c_gfrac$X01 <- raw_gfrac$X01_sensi.contigs
x1000c_gfrac$X00 <- raw_gfrac$X00_sensi.contigs
x1000c_gfrac <- x1000c_gfrac[,c(11,10,1:9)]

x700c_gfrac$X00 <- raw_gfrac$X00_sensi.contigs
x700c_gfrac <- x700c_gfrac[,c(11,1:10)]

maxbin_gfrac$X00 <- raw_gfrac$X00_sensi.contigs
maxbin_gfrac <- maxbin_gfrac[,c(11,1:10)]


# Clean up data and add columns
gfrac <- clean_data(gfrac, "gsa")
raw_gfrac <- clean_data(raw_gfrac, "raw")
x1000c_gfrac <- clean_data(x1000c_gfrac, "x1000c")
x700c_gfrac <- clean_data(x700c_gfrac, "x700c")
maxbin_gfrac <- clean_data(maxbin_gfrac, "maxbin")


# STANDARD VECTORS
# HC percentages
hc <- c(0, 0.254, 0.429, 0.566, 0.664, 0.759, 0.818, 0.878, 0.923, 0.964, 0.983)
# Sample names
sample_id <- c("X00", "X01", "X02", "X03", "X04", "X05", "X06", "X07",
               "X08", "X090", "X095")


gfrac[is.na(gfrac)] <- 0
raw_gfrac[is.na(raw_gfrac)] <- 0
x1000c_gfrac[is.na(x1000c_gfrac)] <- 0
x700c_gfrac[is.na(x700c_gfrac)] <- 0 
maxbin_gfrac[is.na(maxbin_gfrac)] <- 0


# Index out the groups
endo <- gfrac[gfrac$group != "host",1:11]
host <- gfrac[gfrac$group == "host",1:11]
mock <- gfrac[,1:11]
rf <- gfrac[gfrac$group == "rfungi",1:11]
omf <- gfrac[gfrac$group == "OMF",1:11]
fungi <- gfrac[gfrac$group %in% c("rfungi", "OMF"), 1:11]
ba_ar <- gfrac[gfrac$group == "ba_ar",1:11]
pl_vi_unk <- gfrac[gfrac$group == "pl_vi_unk",1:11]

raw_endo <- raw_gfrac[raw_gfrac$group != "host",1:11]
raw_host <- raw_gfrac[raw_gfrac$group == "host",1:11]
raw_mock <- raw_gfrac[,1:11]
raw_rf <- raw_gfrac[raw_gfrac$group == "rfungi",1:11]
raw_omf <- raw_gfrac[raw_gfrac$group == "OMF",1:11]
raw_fungi <- raw_gfrac[raw_gfrac$group %in% c("rfungi", "OMF"), 1:11]
raw_ba_ar <- raw_gfrac[raw_gfrac$group == "ba_ar",1:11]
raw_pl_vi_unk <- raw_gfrac[raw_gfrac$group == "pl_vi_unk",1:11]

x1000c_endo <- x1000c_gfrac[x1000c_gfrac$group != "host",1:11]
x1000c_host <- x1000c_gfrac[x1000c_gfrac$group == "host",1:11]
x1000c_mock <- x1000c_gfrac[,1:11]
x1000c_rf <- x1000c_gfrac[x1000c_gfrac$group == "rfungi",1:11]
x1000c_omf <- x1000c_gfrac[x1000c_gfrac$group == "OMF",1:11]
x1000c_fungi <- x1000c_gfrac[x1000c_gfrac$group %in% c("rfungi", "OMF"), 1:11]
x1000c_ba_ar <- x1000c_gfrac[x1000c_gfrac$group == "ba_ar",1:11]
x1000c_pl_vi_unk <- x1000c_gfrac[x1000c_gfrac$group == "pl_vi_unk",1:11]

x700c_endo <- x700c_gfrac[x700c_gfrac$group != "host",1:11]
x700c_host <- x700c_gfrac[x700c_gfrac$group == "host",1:11]
x700c_mock <- x700c_gfrac[,1:11]
x700c_rf <- x700c_gfrac[x700c_gfrac$group == "rfungi",1:11]
x700c_omf <- x700c_gfrac[x700c_gfrac$group == "OMF",1:11]
x700c_fungi <- x700c_gfrac[x700c_gfrac$group %in% c("rfungi", "OMF"), 1:11]
x700c_ba_ar <- x700c_gfrac[x700c_gfrac$group == "ba_ar",1:11]
x700c_pl_vi_unk <- x700c_gfrac[x700c_gfrac$group == "pl_vi_unk",1:11]

maxbin_endo <- maxbin_gfrac[maxbin_gfrac$group != "host",1:11]
maxbin_host <- maxbin_gfrac[maxbin_gfrac$group == "host",1:11]
maxbin_mock <- maxbin_gfrac[,1:11]
maxbin_rf <- maxbin_gfrac[maxbin_gfrac$group == "rfungi",1:11]
maxbin_omf <- maxbin_gfrac[maxbin_gfrac$group == "OMF",1:11]
maxbin_fungi <- maxbin_gfrac[maxbin_gfrac$group %in% c("rfungi", "OMF"), 1:11]
maxbin_ba_ar <- maxbin_gfrac[maxbin_gfrac$group == "ba_ar",1:11]
maxbin_pl_vi_unk <- maxbin_gfrac[maxbin_gfrac$group == "pl_vi_unk",1:11]


################## Analysis #######################
# GSA
median <- apply(mock, 2, median, na.rm = TRUE)
shapiro.test(median) # W = 0.73521, p-value = 0.001334
# Correlation
cor.test(hc, median, method ="kendall", exact = TRUE) 
#p-value = 5.01e-08    r = -1

# RAW
raw_median <- apply(raw_mock, 2, median, na.rm = TRUE)
shapiro.test(raw_median)  # W = 0.54954, p-value = 6.062e-06 Not norm
# Correlation
cor.test(hc, raw_median, method ="kendall", exact = FALSE) 
# p-value = 0.001152   r = -0.8101627 

# x1000c
x1000c_median <- apply(x1000c_mock, 2, median, na.rm = TRUE)
shapiro.test(x1000c_median)  # W = 0.53515, p-value = 4.044e-06 Not norm
# Correlation
cor.test(hc, x1000c_median, method ="kendall", exact = FALSE) 
# p-value = 0.001152   r = -0.8101627  

# x700c
x700c_median <- apply(x700c_mock, 2, median, na.rm = TRUE)
shapiro.test(x700c_median)  # W = 0.54025, p-value = 4.668e-06 Not norm
# Correlation
cor.test(hc, x700c_median, method ="kendall", exact = FALSE) 
# p-value = 0.0006225  r = -0.8528029  

# Maxbin
maxbin_median <- apply(maxbin_mock, 2, median, na.rm = TRUE)
shapiro.test(maxbin_median)  # W = 0.53609, p-value = 4.152e-06 not norm
# Correlation
cor.test(hc, maxbin_median, method ="kendall", exact = FALSE) 
# p-value = 0.001152   r = -0.8101627 



################## BOXPLOT ############
# Convert to long-format
gsa_long <- pivot_longer(data = gfrac, cols = colnames(gfrac)[1:11],
                         names_to = "sample", values_to="value")

raw_long <- pivot_longer(data = raw_gfrac, cols = colnames(raw_gfrac)[1:11],
                         names_to = "sample", values_to="value")

x1000c_long <- pivot_longer(data = x1000c_gfrac, cols = colnames(x1000c_gfrac)[1:11],
                            names_to = "sample", values_to="value")

x700c_long <- pivot_longer(data = x700c_gfrac, cols = colnames(x700c_gfrac)[1:11],
                           names_to = "sample", values_to="value")

maxbin_long <- pivot_longer(data = maxbin_gfrac, cols = colnames(maxbin_gfrac)[1:11],
                            names_to = "sample", values_to="value")

# View the first few rows
head(raw_long)
head(gsa_long)
head(x1000c_long)
head(x700c_long)
head(maxbin_long)

merged_long <- rbind(gsa_long, raw_long, x1000c_long, x700c_long, maxbin_long)
merged_long$repli <- factor(merged_long$repli, levels = c("gsa", "raw", "x1000c",
                                                          "x700c", "maxbin"))
merged_long$hc <- hc[match(merged_long$sample, sample_id)]


# Boxplot
gfrac_p1 <- ggplot(merged_long, aes(x=sample, y = value, fill = repli)) +
  geom_boxplot(alpha=0.4, position = "dodge") +
  theme_bw() +
  theme(axis.text.x=element_blank(),legend.position = "top") +
  labs(y = "Fraction (%)", x = "")


# Boxplot
gfrac_p2 <- ggplot(merged_long[merged_long$repli != "gsa",], aes(x=sample, y = value, fill = repli)) +
  geom_boxplot(alpha=0.4, position = "dodge") +
  ylim(0,20) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "Fraction (%)", subtitle = "GSA excluded", y = "")

ggarrange(gfrac_p1, gfrac_p2, labels = c("A", "B"), ncol = 1, nrow = 2)




######################### MEDIAN PLOTTING ####################################

# Building the data.frame
hc_long <- rep(round(hc,2), 5)
sample_id_long <- rep(sample_id, 5)
medians <- c(median, raw_median, x1000c_median, x700c_median, maxbin_median)
repli <- vector()
for (name in unique(merged_long$repli)) {
  repli <- c(repli, rep(name, 11))
}
median_long <- data.frame(repli, hc = hc_long, sample = sample_id_long, median = medians)


# Spaghetti plot
median_long$repli <- factor(median_long$repli, levels = c("gsa", "raw", "x1000c",
                                                          "x700c", "maxbin"))
ggplot(median_long, aes(x=hc, y = median, color = repli, linetype = repli)) +
  geom_line(lwd=0.9, alpha = 0.6) +
  theme_bw() +
  scale_color_manual(values = c("darkgreen", "red", "orange", "purple", "blue"))



# PLOTTING SEPARATE GRAPHS TO MIND THE SPAGHETTI
# MAXBIN
median_long$repli <- factor(median_long$repli, levels = c("gsa", "raw",
                                                          "x1000c", "x700c", "maxbin"))
maxbin_p <- ggplot(median_long[median_long$repli != "gsa",],
                   aes(x=hc, y = median, group = repli)) +
  geom_line(aes(linetype = repli, color = repli, size = repli)) +
  scale_linetype_manual(values=c(1,1,1,1)) +
  scale_color_manual(values=c("grey", "grey", "grey", "darkgreen")) +
  scale_size_manual(values=c(1, 1, 1, 0.8)) +
  scale_x_continuous(breaks=hc_long) + 
  theme_bw() + 
  labs(subtitletitle = "maxbin", y = "") +
  theme(axis.text.x=element_text(angle=-45, hjust = 0, vjust = 0, size = 8), 
        axis.ticks.y=element_blank(),axis.text.y=element_blank(),
        legend.position = "none") 


# RAW
median_long$repli <- factor(median_long$repli, levels = c("gsa", "x1000c",
                                                          "x700c", "maxbin", "raw"))
raw_p <- ggplot(median_long[median_long$repli != "gsa",],
                aes(x=hc, y = median, group = repli)) +
  geom_line(aes(linetype = repli, color = repli, size = repli)) +
  scale_linetype_manual(values=c(1,1,1,1)) +
  scale_color_manual(values=c("grey", "grey", "grey", "darkgreen")) +
  scale_size_manual(values=c(1, 1, 1, 0.8)) +
  scale_x_continuous(breaks=hc_long) + 
  theme_bw() +
  labs(subtitletitle = "raw", x = "") +
  theme(axis.text.x=element_blank(), legend.position = "none") 

# x1000c
median_long$repli <- factor(median_long$repli, levels = c("gsa", "raw", "x700c",
                                                          "maxbin", "x1000c"))
x1000c_p <- ggplot(median_long[median_long$repli != "gsa",],
                   aes(x=hc, y = median, group = repli)) +
  geom_line(aes(linetype = repli, color = repli, size = repli)) +
  scale_linetype_manual(values=c(1,1,1,1)) +
  scale_color_manual(values=c("grey", "grey", "grey", "darkgreen")) +
  scale_size_manual(values=c(1, 1, 1, 0.8)) +
  scale_x_continuous(breaks=hc_long) + 
  theme_bw() +
  labs(subtitletitle = "x1000c", x = "", y = "") +
  theme(axis.text.x=element_blank(), axis.ticks.y=element_blank(),
        axis.text.y=element_blank(), legend.position = "none") 


# x700c
median_long$repli <- factor(median_long$repli, levels = c("gsa", "raw",
                                                          "x1000c", "maxbin","x700c"))
x700c_p <- ggplot(median_long[median_long$repli != "gsa",],
                  aes(x=hc, y = median, group = repli)) +
  geom_line(aes(linetype = repli, color = repli, size = repli)) +
  scale_linetype_manual(values=c(1,1,1,1)) +
  scale_color_manual(values=c("grey", "grey", "grey", "darkgreen")) +
  scale_size_manual(values=c(1, 1, 1, 0.8)) +
  scale_x_continuous(breaks=hc_long) + 
  theme_bw() + 
  labs(subtitletitle = "x700c") +
  theme(axis.text.x=element_text(angle=-45, hjust = 0, vjust = 0, size = 8), 
        legend.position = "none") 

ggarrange(raw_p, x1000c_p, x700c_p, maxbin_p, labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

