library(ggridges)
library(ggplot2)
library(RColorBrewer)

################################################################################
# R module to analyse the distribution of contigs length for GSA and RAW
################################################################################


############################### GSA ############################################

# Read in the data and build table in long format
gsa_00 <- read.table("contigs/gsa/00_gsa.fasta_stats.tsv", sep = "\t",
                     header = FALSE)
gsa_01 <- read.table("contigs/gsa/01_gsa.fasta_stats.tsv", sep = "\t",
                     header = FALSE)
gsa_02 <- read.table("contigs/gsa/02_gsa.fasta_stats.tsv", sep = "\t",
                     header = FALSE)
gsa_03 <- read.table("contigs/gsa/03_gsa.fasta_stats.tsv", sep = "\t",
                     header = FALSE)
gsa_04 <- read.table("contigs/gsa/04_gsa.fasta_stats.tsv", sep = "\t",
                     header = FALSE)
gsa_05 <- read.table("contigs/gsa/05_gsa.fasta_stats.tsv", sep = "\t",
                     header = FALSE)
gsa_06 <- read.table("contigs/gsa/06_gsa.fasta_stats.tsv", sep = "\t",
                     header = FALSE)
gsa_07 <- read.table("contigs/gsa/07_gsa.fasta_stats.tsv", sep = "\t",
                     header = FALSE)
gsa_08 <- read.table("contigs/gsa/08_gsa.fasta_stats.tsv", sep = "\t",
                     header = FALSE)
gsa_090 <- read.table("contigs/gsa/090_gsa.fasta_stats.tsv", sep = "\t",
                      header = FALSE)
gsa_095 <- read.table("contigs/gsa/095_gsa.fasta_stats.tsv", sep = "\t",
                      header = FALSE)

gsa_00$sample <- rep("X00", nrow(gsa_00))
gsa_01$sample <- rep("X01", nrow(gsa_01))
gsa_02$sample <- rep("X02", nrow(gsa_02))
gsa_03$sample <- rep("X03", nrow(gsa_03))
gsa_04$sample <- rep("X04", nrow(gsa_04))
gsa_05$sample <- rep("X05", nrow(gsa_05))
gsa_06$sample <- rep("X06", nrow(gsa_06))
gsa_07$sample <- rep("X07", nrow(gsa_07))
gsa_08$sample <- rep("X08", nrow(gsa_08))
gsa_090$sample <- rep("X90", nrow(gsa_090))
gsa_095$sample <- rep("X95", nrow(gsa_095))

# Generate long table
gsa_long <- data.frame(length = c(gsa_00$V2, gsa_01$V2, gsa_02$V2, gsa_03$V2, 
                                  gsa_04$V2, gsa_05$V2, gsa_06$V2, gsa_07$V2, 
                                  gsa_08$V2, gsa_090$V2, gsa_095$V2),
                       sample = c(gsa_00$sample, gsa_01$sample, gsa_02$sample, 
                                  gsa_03$sample, gsa_04$sample, gsa_05$sample,
                                  gsa_06$sample, gsa_07$sample, gsa_08$sample,
                                  gsa_090$sample, gsa_095$sample))

gsa_long$repli <- rep("gsa", nrow(gsa_long))



# Analysing common statistics for distribution
# Checking subset of samples
median(gsa_00$V2) # 274
IQR(gsa_00$V2) # 131


median(gsa_01$V2) # 263
IQR(gsa_01$V2) # 57

median(gsa_02$V2) # 264
IQR(gsa_02$V2) # 58


median(gsa_090$V2) # 270
IQR(gsa_090$V2) # 74

# Checking compiled version of all samples
median(gsa_long$length)  # 268
IQR(gsa_long$length)   # 62
range(gsa_long$length) # [148 155810]
length(gsa_long$length) # 45130914

t <- gsa_long[gsa_long$length <=150,]
t2 <- gsa_long[gsa_long$length <=300,]


# The amount of data that is removed if we use min length of 1000.
1 - (nrow(gsa_02[gsa_02$V2 > 1000,]) / nrow(gsa_02) )   # 98% is discarded
1 - (nrow(gsa_02[gsa_02$V2 > 700,]) / nrow(gsa_02) )   # 96.96% is discarded
1 - (nrow(gsa_02[gsa_02$V2 > 500,]) / nrow(gsa_02) )   # 96.96% is discarded



nrow(gsa_02[gsa_02$V2 <= 150,])
nrow(gsa_02)

# nrow(gsa_02[gsa_02$V2 > 700,])
# [1] 105841
# > nrow(gsa_02[gsa_02$V2 < 700,])
# [1] 3387000
# > nrow(gsa_02[gsa_02$V2 > 1000,])
# [1] 67734
# > nrow(gsa_02[gsa_02$V2 < 1000,])
# [1] 3425226

# Vector holding the number of contigs for each sample.
gsa_length <- c(length(gsa_00$V2), length(gsa_01$V2), length(gsa_02$V2),
                    length(gsa_03$V2), length(gsa_04$V2), length(gsa_05$V2),
                    length(gsa_06$V2), length(gsa_07$V2), length(gsa_08$V2),
                    length(gsa_090$V2), length(gsa_095$V2))

median(gsa_length) # 4624252

# Ridge plot showing that the distribution is very similar with two spikes at
# 150 and 250-300.
#ggplot(gsa_long, aes(x = length, y = sample, fill = sample)) +
  geom_density_ridges() +
  scale_fill_brewer(palette = "Spectral") +
  theme_ridges() + 
  theme(legend.position = "none") +
  labs(title = "Distribution of GSA contigs") +
  scale_x_continuous(breaks = seq(0,750, by = 50), limits = c(100,600)) 
  



############################### RAW ############################################


raw_00 <- read.table("contigs/raw/00_sensi.contigs.fa_stats.tsv", sep = "\t",
                     header = FALSE)
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


raw_00$sample <- rep("X00", nrow(raw_00))
raw_01$sample <- rep("X01", nrow(raw_01))
raw_02$sample <- rep("X02", nrow(raw_02))
raw_03$sample <- rep("X03", nrow(raw_03))
raw_04$sample <- rep("X04", nrow(raw_04))
raw_05$sample <- rep("X05", nrow(raw_05))
raw_06$sample <- rep("X06", nrow(raw_06))
raw_07$sample <- rep("X07", nrow(raw_07))
raw_08$sample <- rep("X08", nrow(raw_08))
raw_090$sample <- rep("X90", nrow(raw_090))
raw_095$sample <- rep("X95", nrow(raw_095))


raw_long <- data.frame(length = c(raw_00$V2, raw_01$V2, raw_02$V2, raw_03$V2,
                                  raw_04$V2, raw_05$V2, raw_06$V2, raw_07$V2,
                                  raw_08$V2, raw_090$V2, raw_095$V2),
                       sample = c(raw_00$sample, raw_01$sample, raw_02$sample,
                                  raw_03$sample, raw_04$sample, raw_05$sample,
                                  raw_06$sample, raw_07$sample, raw_08$sample,
                                  raw_090$sample, raw_095$sample))
raw_long$repli <- rep("raw", nrow(raw_long))


median(raw_long$length)  # 418
IQR(raw_long$length)   # 218
range(raw_long$length) # 200 16703
length(raw_long$length) # 481587

1 - (nrow(raw_long[raw_long$length > 500,]) / length(raw_long$length)) 
# 84.8% discarded


raw_lengths <- c(length(raw_00$V2), length(raw_01$V2), length(raw_02$V2),
                    length(raw_03$V2), length(raw_04$V2), length(raw_05$V2),
                    length(raw_06$V2), length(raw_07$V2), length(raw_08$V2),
                    length(raw_090$V2), length(raw_095$V2))


median(raw_lengths) # 29789



# Ridge plot showing that the distribution is very similar with two spikes at
# 150 and 250-300.
ggplot(raw_long, aes(x = length, y = sample, fill = sample)) +
  geom_density_ridges() +
  scale_fill_brewer(palette = "Spectral") +
  theme_ridges() + 
  theme(legend.position = "none") +
  labs(title = "Distribution of raw contigs") +
  scale_x_continuous(breaks = seq(0,17000, by = 200), limits = c(150,1500)) 


# Raw VS GSA

rvg_long <- rbind(raw_long, gsa_long)

rvg_long$log.len <- log(rvg_long$length)

# Violing plot, with the boxplot it kinda looks meh.
ggplot(rvg_long[rvg_long$sample == "X05",], aes(x=repli, y=length, fill = repli)) + 
  geom_violin() +
  #scale_fill_brewer(palette = "Paired") +
  theme_bw() + 
  theme(legend.position = "none") +
  labs(title = "Distribution of contigs", x = "replicate") 
  # scale_y_continuous(breaks = seq(0,2000, by = 200), limits = c(100,1200)) 
  # geom_boxplot(width=0.015, fill = "lightgrey")


# Much better representation i think.
ggplot(rvg_long[rvg_long$sample == "X05",], aes(x = length, y = repli, fill = repli)) +
  geom_density_ridges() +
  scale_fill_brewer(palette = "Paired") +
  theme_ridges() + 
  theme(legend.position = "none") +
  labs(title = "Contig length distribution", y = "", x = "length (bp)") +
  scale_x_continuous(breaks = seq(0,1200, by = 100), limits = c(100,1200)) 



########### host MAG stats #########



x700 <- c(262531, 3500344, 11378057, 20412080, 24598053,
          31683965, 49710402, 57299033, 60263965, 63913566)
shapiro.test(x700) # norm
mean(x700)  # 32302200 32.3M

x1000c <-c(0, 99268, 354516, 808275, 1148938, 4917286,
           6112569, 7392535, 7785221, 8519125)

shapiro.test(x1000c) # not norm
mean(x1000c)
median(x1000c) # 3033112



######### playground ######
vector <- c(282042138, 279746860, 269473090, 268284664, 
            246193872, 228527597, 234157245, 221050550, 
            217457960, 224569803, 231454572)

boxplot(vector)
sd(vector)
mean(vector)
shapiro.test(vector)

