# ANALYSIS TO DETERMINE CORRELATION

# Packages
library(ggplot2)
library(dplyr)


# Missmatches per 100kbp
missamble <- read.table("assmble/gsa/TSV/gsa_num_misassemblies.tsv",
                        sep = "\t", header = TRUE, row.names = 1,
                        na.strings = "-")

# Replace NA with 0
#missamble[is.na(missamble)] <- 0


dim(missamble)   # 100 11

# Remove ".contigs" suffix from row names 
rownames(missamble) <- gsub(".contigs", "", rownames(missamble))
# Remove numbers and dot prefix from row names
rownames(missamble) <- sub("^\\d+\\.", "", rownames(missamble))

# Add the grouping from mock_data.
# Order rownames
missamble <- missamble[order(rownames(missamble)),]

# Mock dataframe 
mock_df <- read.table("mock_genomes.txt", sep = "\t", header = TRUE, 
                      row.names = 1)
# Order row names
mock_df <- mock_df[order(rownames(mock_df)),]

# Add group and order gsa
missamble$group <- mock_df$group

missamble <- missamble[order(missamble$group),]

# Looks good
head(missamble)

# HC percentages
hc <- c(0, 0.254, 0.429, 0.566, 0.664, 0.759, 0.818, 0.878, 0.923, 0.964, 0.983)


# Analysis

shapiro.test(missamble$X00_gsa)    # Not normally dist, tested more samples

range(missamble$X095_gsa, na.rm = TRUE)

endo <- missamble[missamble$group != "host",1:11]
mock <- missamble[,1:11]
rf <- missamble[missamble$group == "rfungi",1:11]
omf <- missamble[missamble$group == "OMF",1:11]
fungi <- missamble[missamble$group %in% c("rfungi", "OMF"), 1:11]
ba_ar <- missamble[missamble$group == "ba_ar",1:11]
pl_vi_unk <- missamble[missamble$group == "pl_vi_unk",1:11]
sample_id <- c("00", "01", "02", "03", "04", "05", "06", "07",
               "08", "090", "095")

# Endophytes
boxplot(endo, main = "Endophytes", names = sample_id) # we see that it goes down to 1.

# par(mfrow = c(2, 3))

# All
boxplot(mock, main = "mock data",names = sample_id) # we see that it goes down to 1.

# rfungi
boxplot(rf, main = "rfungi", names = sample_id)

# OMF
boxplot(omf, main = "OMF",names = sample_id)

# Fungi
boxplot(fungi, main = "Fungi",names = sample_id)

# ba_ar
boxplot(ba_ar, main = "ba_ar",names = sample_id)

# pl_vi_unk
boxplot(pl_vi_unk, main = "pl_vi_unk",names = sample_id)

par(mfrow = c(1, 1))

boxplot(missamble[missamble$group != "pl_vi_unk",1:11], names = sample_id,
        main = "Duplication ratio\npl_vi_unk excluded")

# mock
median <- apply(mock, 2, median, na.rm = TRUE)
variance <- apply(mock, 2, var, na.rm = TRUE)
sd <- apply(mock, 2, sd, na.rm = TRUE)
IQR <- apply(mock, 2, function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  return(IQR)
})

stats <- data.frame(median = median, variance = variance, sd = sd,
                    sample = sample_id, iqr = IQR)

plot(hc, median, main = "mock df\nMedian", xlab = "Host Contamination (%)")
plot(hc, variance, main = "mock df\nVariance")

plot(hc, stats$median, pch = 16, cex = 1.1, ylim = c(0, 1), ylab = "Median",
     main = "Median for GSA\nIQR Whiskers", xlab = "Host Contamination (%)")
arrows(x0=hc, y0=stats$median-stats$iqr, x1=hc, y1=stats$median+stats$iqr,
       code = 3, angle = 90, length = 0.1)


ggplot(stats) +
  geom_bar( aes(x=sample, y = median), stat = "identity", fill = "skyblue",
            alpha = 0.7) +
  geom_errorbar (aes(x=sample, ymin=median-sd, ymax=median+sd), width=0.4, 
                 colour="orange", alpha=0.9, size=1)
  


par(mfrow = c(1,2))

#Host
plot(hc, missamble[missamble$group == "host",1:11],
     main = "Num of Misassemblies\nof the\nHost genome", ylab = "Misassemblies", 
     xlab = "Proportion of Host Contamination")

plot(hc, stats$median, pch = 16, cex = 1.1, ylim = c(0, 1), ylab = "Median",
     main = "Num of Misassemblies\nof\nMock Data",
     xlab = "Proportion of Host Contamination")
arrows(x0=hc, y0=stats$median-stats$iqr, x1=hc, y1=stats$median+stats$iqr,
       code = 3, angle = 90, length = 0.1)

# Kendall's tau
cor.test(hc, median, method ="kendall", exact = T)

sd(median)  # 0
# No correlation




