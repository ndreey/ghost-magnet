library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(tidyr)
source("analysis_tools.R")

samples <- c("X00", "X01", "X02", "X03", "X04", "X05", "X06", "X07",
             "X08", "X090", "X095")
hc <- c(0, 0.25, 0.43, 0.57, 0.66, 0.76, 0.81, 0.88, 0.92, 0.96, 0.98)

# Store the bin vectors holding the most abundant organism group.
# CAREFUL
samp_list <- list()
repli1 <- vector()
sample1 <- vector()
counts1 <- vector()
group1 <- vector()

for (sample in samples) {
  
  filename <- paste0("binning/maxbin_host/", sample, "_blast.tsv")
  cat("Starting: ", filename, "\n")
  
  maxbin_bins <- read.table(filename, sep = "\t", header = FALSE)
  
  colnames(maxbin_bins) <- c("qseqid", "sseqid", "pident", "evalue", "bitscore", 
                             "bin")
  # RAW X00 
  res <- vector()
  for (bin in unique(maxbin_bins$bin)) {
    # Gets so each unique qseqid matches to highest score qseq.
    abundance <- get_highest_score(maxbin_bins[maxbin_bins$bin == bin,])
    
    # Returns the most abundant group in the bin
    bin_group <- get_abundant(abundance)
    
    # If its even steven, take the one with highest score
    if (length(bin_group) >= 2) {
      tmp <- abundance[abundance$bitscore == max(abundance$bitscore),2]
      bin_group <- regex_split(tmp$sseqid)
    }
    
    res <- c(res, bin_group)
  }
  
  if (length(unique(maxbin_bins$bin)) < length(res)) {
    cat("WARNING    ", sample)
  }
  
  
  
  samp_list[[sample]] <- res
  cat("Complete: ", filename, "\n")
}


# Change me                <<<<<<<<
replicate <- "maxbin"

for (i in 1:11) {
  
  for (grp in c("host", "rfungi", "OMF", "ba_ar", "pl_vi_unk")) {
    counts1 <- c(counts1, sum(samp_list[[i]] == grp))
    group1 <- c(group1, grp)
    sample1 <- c(sample1, samples[i])
    repli1 <- c(repli1, replicate)
    
  }
  
}








##################### LONG TABLE ARENA ##########

# For proportion of each group
df_long <- data.frame(repli = repli1, sample = sample1, counts = counts1,
                      group = group1)

df_long <- read.table("binning/df_long_mockdb.tsv", sep = "\t", header = TRUE)

df_long$group <- factor(df_long$group, levels = c("rfungi", "host", "OMF", "ba_ar",
                                                  "pl_vi_unk"))
# Write data frame to TSV file
#write.table(df_long, "df_long_mockdb.tsv", sep = "\t", quote = FALSE, 
#            row.names = FALSE)




# For count line plot
host_long <- df_long[df_long$group == "host",]

# Subset of the data frame
subset_df <- df_long[df_long$group != "host",]

# Sum counts by sample
sum_counts <- summarise(group_by(subset_df, repli, sample), total_counts = sum(counts))

colnames(sum_counts) <- c("repli", "sample", "counts")


count_long <- data.frame(repli = c(sum_counts$repli, host_long$repli), 
                         sample = c(sum_counts$sample, host_long$sample),
                         counts = c(sum_counts$counts, host_long$counts),
                         group = c(rep("endophyte",nrow(sum_counts)),
                                   as.character(host_long$group)),
                         hc = c(hc, hc))

count_long$repli <- factor(count_long$repli,
                           levels = c("raw", "x1000c", "x700c", "maxbin"))

count_long$group <- factor(count_long$group, levels = c("endophyte", "host"))



# For proportion of endophytes
prop_df <- data.frame(repli = character(),
                      prop = numeric(),
                      n_bin = numeric(),
                      n_endo = numeric(),
                      sample = character(),
                      stringsAsFactors = FALSE)


for (r in c("raw", "x1000c", "x700c", "maxbin")) {
  cat("Repli:", r, "\n")
  
  for (s in samples) {
    cat("    Sample:", s, "\n")
    
    n_bin <- sum(count_long[count_long$repli == r & count_long$sample == s, "counts"])
    n_endo <- sum(count_long[count_long$repli == r & count_long$sample == s & count_long$group == "endophyte", "counts"])
    prop <- n_endo / n_bin
    
    prop_df <- rbind(prop_df, data.frame(repli = r,
                                         prop = prop,
                                         n_bin = n_bin,
                                         n_endo = n_endo,
                                         sample = s))
  }
}

prop_df$repli <- factor(prop_df$repli, levels = c("raw", "x1000c", "x700c", "maxbin"))
prop_df$hc <- rep(hc, 4)

############### PLOT ARENA #################


# Plot of percentage stack plot
pp_raw <- ggplot(df_long[df_long$repli == "raw",], aes(fill=group, y=counts, x=sample)) + 
  geom_bar(position="fill", stat="identity", alpha = 1) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  labs(subtitle = "raw", y = "Fraction", x = "") +
  theme(axis.text.x=element_blank(), legend.position = "none", legend.title=element_blank()) 

pp_x1000c <- ggplot(df_long[df_long$repli == "x1000c",], aes(fill=group, y=counts, x=sample)) + 
  geom_bar(position="fill", stat="identity", alpha = 1) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  labs(subtitle = "x1000c", y = "", x = "") + 
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), 
        axis.ticks.x=element_blank(), axis.text.x=element_blank(),
        legend.position = "none", legend.title=element_blank())

pp_x700c <- ggplot(df_long[df_long$repli == "x700c",], aes(fill=group, y=counts, x=sample)) + 
  geom_bar(position="fill", stat="identity", alpha = 1) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  labs(subtitle = "x700c", y = "Fraction", x = "sample") +
  theme(legend.position = "none", legend.title=element_blank())

pp_max <- ggplot(df_long[df_long$repli == "maxbin",], aes(fill=group, y=counts, x=sample)) + 
  geom_bar(position="fill", stat="identity", alpha = 1) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  labs(subtitle = "maxbin", y = "", x = "sample") +
  theme(legend.position = "none", axis.ticks.y=element_blank(), 
        axis.text.y=element_blank(), legend.title=element_blank())

ggarrange(pp_raw, pp_x1000c, pp_x700c, pp_max, labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")



# Plot line plot
lp_raw <- ggplot(count_long[count_long$group == "endophyte",], 
       aes(x = hc, y = counts, color = repli)) + 
  geom_line(aes(linetype = repli, color = repli, size = repli)) +
  scale_x_continuous(breaks=hc) +
  scale_linetype_manual(values=c(1, 1, 1, 1)) +
  scale_color_manual(values=c("darkgreen", "grey", "grey", "grey")) +
  scale_size_manual(values=c(1, 1, 1,1)) +
  theme_bw() + 
  labs(subtitle = "raw", x = "") +
  theme(axis.text.x=element_blank(), legend.position = "none", legend.title=element_blank())

lp_x1000c <- ggplot(count_long[count_long$group == "endophyte",], 
                 aes(x = hc, y = counts, color = repli)) + 
  geom_line(aes(linetype = repli, color = repli, size = repli)) +
  scale_x_continuous(breaks=hc) +
  scale_linetype_manual(values=c(1, 1, 1, 1)) +
  scale_color_manual(values=c("grey", "darkgreen", "grey", "grey")) +
  scale_size_manual(values=c(1, 1, 1,1)) +
  theme_bw() + 
  labs(subtitle = "x1000c", x="", y="") +
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), 
        axis.ticks.x=element_blank(), axis.text.x=element_blank(),
        legend.position = "none", legend.title=element_blank())

lp_x700c <- ggplot(count_long[count_long$group == "endophyte",], 
                 aes(x = hc, y = counts, color = repli)) + 
  geom_line(aes(linetype = repli, color = repli, size = repli)) +
  scale_x_continuous(breaks=hc) +
  scale_linetype_manual(values=c(1, 1, 1, 1)) +
  scale_color_manual(values=c("grey", "grey", "darkgreen", "grey")) +
  scale_size_manual(values=c(1, 1, 1,1)) +
  theme_bw() + 
  labs(subtitle = "x700c") +
  theme(legend.position = "none", legend.title=element_blank(),
        axis.text.x=element_text(angle=-45, hjust = 0, vjust = 0, size = 8.5))

lp_max <- ggplot(count_long[count_long$group == "endophyte",], 
                 aes(x = hc, y = counts, color = repli)) + 
  geom_line(aes(linetype = repli, color = repli, size = repli)) +
  scale_x_continuous(breaks=hc) +
  scale_linetype_manual(values=c(1, 1, 1, 1)) +
  scale_color_manual(values=c("grey", "grey", "grey", "darkgreen")) +
  scale_size_manual(values=c(1, 1, 1,1)) +
  theme_bw() + 
  labs(subtitle = "maxbin", y = "") +
  theme(legend.position = "none", axis.ticks.y=element_blank(), 
        axis.text.y=element_blank(),legend.title=element_blank(),
        axis.text.x=element_text(angle=-45, hjust = 0, vjust = 0, size = 8.5))


ggarrange(lp_raw, lp_x1000c, lp_x700c, lp_max, labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, common.legend = FALSE)





bp_raw <- ggplot(count_long[count_long$repli == "raw",], aes(fill=group, y=counts, x=sample)) + 
  geom_bar(position="dodge", stat="identity", alpha = 0.9) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  labs(subtitle = "raw", y = "Counts", x = "") +
  theme(axis.text.x=element_blank(), legend.position = "none", legend.title=element_blank())

bp_x1000c <- ggplot(count_long[count_long$repli == "x1000c",], aes(fill=group, y=counts, x=sample)) + 
  geom_bar(position="dodge", stat="identity", alpha = 1) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  labs(subtitle = "x1000c", y = "", x = "") + 
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), 
        axis.ticks.x=element_blank(), axis.text.x=element_blank(),
        legend.position = "none", legend.title=element_blank()) 

bp_x700c <- ggplot(count_long[count_long$repli == "x700c",], aes(fill=group, y=counts, x=sample)) + 
  geom_bar(position="dodge", stat="identity", alpha = 1) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  labs(subtitle = "x700c", y = "Counts", x = "sample") +
  theme(legend.position = "none", legend.title=element_blank())

bp_max <- ggplot(count_long[count_long$repli == "maxbin",], aes(fill=group, y=counts, x=sample)) + 
  geom_bar(position="dodge", stat="identity", alpha = 1) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  labs(subtitle = "maxbin", y = "", x = "sample") +
  theme(legend.position = "none", axis.ticks.y=element_blank(), 
        axis.text.y=element_blank(), legend.title=element_blank())

ggarrange(bp_raw, bp_x1000c, bp_x700c, bp_max, labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")
 

# Proportion barplot
ggplot(prop_df[prop_df$repli != "maxbin",], aes(x=sample, y = prop, fill = repli)) +
  geom_bar(position="dodge", stat = "identity") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  labs(subtitle = "Proportion of bins that belong to endophyte", y = "Proportion", x = "sample") +
  theme(legend.position = "bottom", legend.title=element_blank())

# Proportion lineplot
plp_raw <- ggplot(prop_df, aes(x=hc, y = prop, color = repli)) +
  geom_line(aes(linetype = repli, color = repli, size = repli)) +
  scale_x_continuous(breaks=hc) +
  scale_linetype_manual(values=c(1, 1, 1, 1)) +
  scale_color_manual(values=c("darkgreen", "grey", "grey", "grey")) +
  scale_size_manual(values=c(1, 1, 1,1)) +
  theme_bw() + 
  labs(subtitle = "raw", x = "") +
  theme(axis.text.x=element_blank(), legend.position = "none", legend.title=element_blank())

plp_x1000c <- ggplot(prop_df, aes(x=hc, y = prop, color = repli)) +
  geom_line(aes(linetype = repli, color = repli, size = repli)) +
  scale_x_continuous(breaks=hc) +
  scale_linetype_manual(values=c(1, 1, 1, 1)) +
  scale_color_manual(values=c("grey", "darkgreen", "grey", "grey")) +
  scale_size_manual(values=c(1, 1, 1,1)) +
  theme_bw() + 
  labs(subtitle = "x1000c", x="", y="") +
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), 
        axis.ticks.x=element_blank(), axis.text.x=element_blank(),
        legend.position = "none", legend.title=element_blank())


plp_x700c <- ggplot(prop_df, aes(x=hc, y = prop, color = repli)) +
  geom_line(aes(linetype = repli, color = repli, size = repli)) +
  scale_x_continuous(breaks=hc) +
  scale_linetype_manual(values=c(1, 1, 1, 1)) +
  scale_color_manual(values=c("grey", "grey", "darkgreen", "grey")) +
  scale_size_manual(values=c(1, 1, 1,1)) +
  theme_bw() + 
  labs(subtitle = "x700c") +
  theme(legend.position = "none", legend.title=element_blank(),
        axis.text.x=element_text(angle=-45, hjust = 0, vjust = 0, size = 8.5))


plp_max <- ggplot(prop_df, aes(x = hc, y = prop, color = repli)) + 
  geom_line(aes(linetype = repli, color = repli, size = repli)) +
  scale_x_continuous(breaks=hc) +
  scale_linetype_manual(values=c(1, 1, 1, 1)) +
  scale_color_manual(values=c("grey", "grey", "grey", "darkgreen")) +
  scale_size_manual(values=c(1, 1, 1,1)) +
  theme_bw() + 
  labs(subtitle = "maxbin", y = "") +
  theme(legend.position = "none", axis.ticks.y=element_blank(), 
        axis.text.y=element_blank(),legend.title=element_blank(),
        axis.text.x=element_text(angle=-45, hjust = 0, vjust = 0, size = 8.5))

ggarrange(plp_raw, plp_x1000c, plp_x700c, plp_max, labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, common.legend = FALSE)


################# STATISTICS ##########################


########## ENDOPHYTE AND HOST BIN COUNT CORRELATION SIGNIFICANCE ##############
raw_host <- count_long[count_long$repli == "raw" & count_long$group == "host",]
shapiro.test(raw_host$counts) # W = 0.85633, p-value = 0.0516 norm dist.
# Correlation
cor.test(hc, raw_host$counts, method ="pearson") 
# p-value = 4.881e-06   r = 0.9547685

raw_endo <- count_long[count_long$repli == "raw" & count_long$group == "endophyte",]
shapiro.test(raw_endo$counts) # W = 0.74899, p-value = 0.002016 not norm
# Correlation
cor.test(hc, raw_endo$counts, method ="kendall", exact = FALSE) 
# p-value = 0.0002058   r = -0.8874245


x1000c_host <- count_long[count_long$repli == "x1000c" & count_long$group == "host",]
shapiro.test(x1000c_host$counts) # W = 0.79299, p-value = 0.007608 not norm dist.
# Correlation
cor.test(hc, x1000c_host$counts, method ="kendall", exact = FALSE) 
# p-value = 0.183   r = 0.3148688 

x1000c_endo <- count_long[count_long$repli == "x1000c" & count_long$group == "endophyte",]
shapiro.test(x1000c_endo$counts) # W = 0.74691, p-value = 0.001894 not norm
# Correlation
cor.test(hc, x1000c_endo$counts, method ="kendall", exact = FALSE) 
# p-value = 0.0001089   r = -0.9251873


x700c_host <- count_long[count_long$repli == "x700c" & count_long$group == "host",]
shapiro.test(x700c_host$counts) # W = 0.91075, p-value = 0.2489 norm dist.
# Correlation
cor.test(hc, x700c_host$counts, method ="pearson") 
# p-value = 0.5157   r = 0.2200149  

x700c_endo <- count_long[count_long$repli == "x700c" & count_long$group == "endophyte",]
shapiro.test(x700c_endo$counts) # W = 0.77317, p-value = 0.004178 not norm
# Correlation
cor.test(hc, x700c_endo$counts, method ="kendall", exact = FALSE) 
# p-value = 0.0001089   r = -0.9251873 


max_host <- count_long[count_long$repli == "maxbin" & count_long$group == "host",]
shapiro.test(max_host$counts) # Error in shapiro.test(max_host$counts) : all 'x' values are identical.
# Correlation NONE

max_endo <- count_long[count_long$repli == "maxbin" & count_long$group == "endophyte",]
shapiro.test(max_endo$counts) # W = 0.87771, p-value = 0.09731 norm
# Correlation
cor.test(hc, max_endo$counts, method ="pearson") 
# p-value = 0.002421   r = -0.811581 


host_cor <- c(0.9547685, 0.3148688, 0.2200149, 0)
host_p <- c(4.881e-06, 0.183, 0.5157, 0)
host_repli <- c("raw", "x1000c", "x700c", "maxbin")
host_df <- data.frame(host_repli, host_cor, host_p)

endo_cor <- c(-0.8874245, -0.9251873, -0.9251873, -0.811581)
endo_p <- c(0.0002058, 0.0001089, 0.0001089, 0.002421)
endo_repli <- c("raw", "x1000c", "x700c", "maxbin")
endo_df <- data.frame(endo_repli, endo_cor, endo_p)


endo_df$endo_repli <- factor(endo_df$endo_repli, levels = c("raw", "x1000c", "x700c", "maxbin"))


############ PROPORTION CORRELATION ##########
prop_df

raw_prop <- prop_df[prop_df$repli == "raw",]
shapiro.test(raw_prop$prop) # W = 0.76515, p-value = 0.003279 not norm dist.
# Correlation
cor.test(hc, raw_prop$prop, method ="kendall", exact = FALSE)

x1000c_prop <- prop_df[prop_df$repli == "x1000c",]
shapiro.test(x1000c_prop$prop) # W = 0.76515, p-value = 0.003279 not norm dist.
# Correlation
cor.test(hc, x1000c_prop$prop, method ="kendall", exact = FALSE)

x700c_prop <- prop_df[prop_df$repli == "x700c",]
shapiro.test(x700c_prop$prop) # W = 0.76515, p-value = 0.003279 not norm dist.
# Correlation
cor.test(hc, x700c_prop$prop, method ="pearson", exact = T)


for (r in c("raw", "x1000c", "x700c", "maxbin")){
  cat("Repli:", r, "\n")
  
  for (s in samples) {
    cat("    Sample:    ", s, "\n")
    
    n_bin <- sum(count_long[count_long$repli == r & count_long$sample == s,]$counts)
    n_endo <- count_long[count_long$repli == r & count_long$sample == s & count_long$group == "endophyte",]$counts
    prop <- n_endo / n_bin
    cat("    Bins:    ", n_bin, "\n")
    cat("    % endo:    ", prop, "\n")
    
  }
}




######################################## PLAYGROUND #################################

playgrnd <- count_long


# Perform the Kruskal-Wallis test for host bins
host_test <- kruskal.test(counts ~ repli, data = subset(count_long, group == "host"))
print(host_test)

# Perform the Kruskal-Wallis test for host bins
host_test <- kruskal.test(counts ~ repli, data = subset(count_long, group == "endophyte"))
print(host_test)

testme <- data.frame(maxbin = count_long[count_long$group == "endophyte" & count_long$repli == "maxbin",3],
            raw = count_long[count_long$group == "endophyte" & count_long$repli == "raw",3],
            x1000c = count_long[count_long$group == "endophyte" & count_long$repli == "x1000c",3],
            x700c = count_long[count_long$group == "endophyte" & count_long$repli == "x700c",3])

wilcox.test(testme$maxbin, testme$raw, exact = FALSE)





df_long[df_long$repli == "raw" & df_long$sample == "X00",]



















# First, we define the data frame
df <- data.frame(
  HC = c(0.00, 25.39, 42.91, 56.64, 66.40, 75.93, 81.78, 87.84, 92.33, 96.42, 98.28),
  OMF = c(56.13, 42.66, 33.09, 24.88, 19.68, 13.45, 10.59, 6.89, 4.44, 2.08, 0.99),
  rfungi = c(16.04, 11.67, 8.77, 6.75, 5.08, 3.87, 2.78, 1.92, 1.17, 0.54, 0.26),
  ba_ar = c(27.62, 20.12, 15.11, 11.63, 8.76, 6.68, 4.80, 3.31, 2.02, 0.94, 0.45),
  pl_vi_unk = c(0.18, 0.13, 0.10, 0.07, 0.06, 0.04, 0.03, 0.02, 0.01, 0.00, 0.00)
)

# Divide all values by 100
df <- df / 100

# For each row, scale by 1 divided by the corresponding OMF value
df <- df / df$OMF

# Print the data frame
print(df)

# Calculate standard deviations
sd_HC <- sd(df$HC)
sd_rfungi <- sd(df$rfungi)
sd_ba_ar <- sd(df$ba_ar)
sd_pl_vi_unk <- sd(df$pl_vi_unk)

print(paste("SD for HC: ", round(sd_HC, 4)))
print(paste("SD for rfungi: ", round(sd_rfungi, 4)))
print(paste("SD for ba_ar: ", round(sd_ba_ar, 4)))
print(paste("SD for pl_vi_unk: ", round(sd_pl_vi_unk, 4)))


boxplot(df[,-c(1,2)])





