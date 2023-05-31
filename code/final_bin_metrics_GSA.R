library(ggplot2)
library(RColorBrewer)
library(ggpubr)
source("analysis_tools.R")

################################################################################
# Analyzing the bin metrics from the binning of the GSA with 700 contig 
# threshold. Metrics are derived from AMBER
################################################################################
samples <- c("X00", "X01", "X02", "X03", "X04", "X05", "X06", "X07",
             "X08", "X090", "X095")
hc <- c(0, 0.25, 0.43, 0.57, 0.66, 0.76, 0.81, 0.88, 0.92, 0.96, 0.98)


# Create an empty list to store results for each sample
prop_long_list <- list()

# Iterate over samples
for (sample in samples) {
  filename <- paste0("gsa_amber/x700c/", sample, "_bin_metrics.tsv")
  cat("Starting: ", filename, "\n")
  
  df <- read.table(filename, sep = "\t", header = TRUE)
  mockdf <- read.table("mock_genomes.txt", sep = "\t", header = TRUE)
  
  df <- df[df$Tool != "Gold standard", c("BINID", "genome_id", "precision_bp", "recall_bp")]
  df <- retain_binid(df)
  df <- mock_groupify(df, mockdf)
  
  # Calculate proportions
  df_prop <- data.frame(table(df$group))
  colnames(df_prop) <- c("group", "freq")
  total_freq <- sum(df_prop$freq)
  
  # Add row for "endophyte" group
  endophyte_freq <- sum(df_prop[df_prop$group != "host", "freq"])
  tmp_endo <- data.frame(group = "endophyte", freq = endophyte_freq)
  df_prop <- rbind(df_prop, tmp_endo)
  
  # Calculate the proportion
  df_prop$proportion <- round(df_prop$freq / total_freq, 4)
  
  # Add sample column
  df_prop$sample <- rep(sample, nrow(df_prop))
  

  # Append the results for the current sample to the list
  prop_long_list[[sample]] <- df_prop
}

# Combine results for all samples into a single dataframe
prop_long <- do.call(rbind, prop_long_list)

# Reset row names
rownames(prop_long) <- NULL

prop_long$group <- factor(prop_long$group,
                          levels = c("endophyte", "host", "OMF", "rfungi", 
                                     "ba_ar", "pl_vi_unk"))

prop_endo_host <- prop_long[prop_long$group == "endophyte" | prop_long$group == "host",]

################################## PLOTS #######################################
group_colors <- brewer.pal(6,"Set2")


pp_group <- ggplot(prop_long[prop_long$group != "endophyte",],
                          aes(fill=group, y=freq, x=sample)) + 
  geom_bar(position="fill", stat="identity", alpha = 1) +
  theme_bw() +
  scale_fill_manual(values=group_colors[2:6]) + 
  labs(subtitle = "Mock Groups", y = "Fraction", x = "sample") +
  theme(legend.position = "none", legend.title=element_blank())

pp_endo_vs_host <- ggplot(prop_endo_host, aes(fill=group, y=freq, x=sample)) + 
  geom_bar(position="fill", stat="identity", alpha = 1) +
  theme_bw() +
  scale_fill_manual(values=group_colors[1:2]) + 
  labs(subtitle = "Endophytes vs Host", y = "", x = "sample") +
  theme(legend.position = "none", axis.ticks.y=element_blank(), 
        axis.text.y=element_blank(), legend.title=element_blank())


ggarrange(pp_group, pp_endo_vs_host, labels = c("A", "B"),
          ncol = 2, nrow = 1, common.legend = TRUE, legend="bottom")

pp_bar_eh <- ggplot(prop_endo_host, aes(fill=group, y=proportion, x=sample)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_manual(values=group_colors[1:2]) +
  labs(subtitle = "Proportion of Host and Endophytic bins ", y = "", 
       x = "sample") +
  theme(legend.position = "none", axis.ticks.y=element_blank(), 
        axis.text.y=element_blank(), legend.title=element_blank())

ggarrange(pp_group,pp_bar_eh, labels = c("A", "B"),
          ncol = 1, nrow = 2, common.legend = TRUE, legend="bottom")



################################## Analysis ####################################

shapiro.test(prop_endo_host[prop_endo_host$group == "endophyte",3])
# W = 0.88176, p-value = 0.1095 norm dist

cor.test(hc, prop_endo_host[prop_endo_host$group == "endophyte",3], method = "pearson")
# p-value = 3.813e-09 r = -0.9908819


