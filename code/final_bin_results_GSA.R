library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(reshape2)
library(cowplot)
source("analysis_tools.R")

################################################################################
# Analyzing the bin metrics from the binning of the GSA with 700 contig 
# threshold. Metrics are derived from AMBER
################################################################################
samples <- c("X00", "X01", "X02", "X03", "X04", "X05", "X06", "X07",
             "X08", "X090", "X095")
hc <- c(0, 0.25, 0.43, 0.57, 0.66, 0.76, 0.81, 0.88, 0.92, 0.96, 0.98)


# Create an empty list to store results for each sample
res_long_list <- list()

# Iterate over samples
for (sample in samples) {
  filename <- paste0("gsa_amber/x700c/", sample, "_results.tsv")
  cat("Starting: ", filename, "\n")
  
  df <- read.table(filename, sep = "\t", header = TRUE)
  # mockdf <- read.table("mock_genomes.txt", sep = "\t", header = TRUE)
  
  df <- df[2, c("Average.completeness..bp.", "Average.purity..bp.",
               "Percentage.of.binned.bp", "Adjusted.Rand.index..bp.")]
  colnames(df) <- c("avg.comp", "avg.purity", "ARI", "prcnt.binned.bp")
  

  # Add sample column
  df$sample <- rep(sample, nrow(df))
  
  
  # Append the results for the current sample to the list
  res_long_list[[sample]] <- df
}

# Combine results for all samples into a single dataframe
res_long <- do.call(rbind, res_long_list)

# Reset row names
rownames(res_long) <- NULL
res_long$hc <- hc

freq <- c(0, 19, 61, 79, 138, 161, 185, 205, 222, 237, 235)
proportion <- c(0, 0.1462, 0.3720, 0.4540, 0.6866, 0.7931, 0.8150, 0.9071, 0.9737, 0.9958, 1.0000)

res_long$freq.bins <- freq
res_long$host.bin.prop <- proportion


####################### ANALYSIS #####################################

shapiro.test(res_long$avg.comp)
x1 <- cor.test(hc, res_long$avg.comp, method = "pearson") # -0.980211 

shapiro.test(res_long$avg.purity)
x2 <- cor.test(hc, res_long$avg.purity, method = "pearson") # 0.9695893

shapiro.test(res_long$ARI)
x3 <- cor.test(hc, res_long$ARI, method = "kendall") # -0.8909091  

shapiro.test(res_long$prcnt.binned.bp)
x4 <- cor.test(hc, res_long$prcnt.binned.bp, method = "pearson") # -0.906889 

# key data
shapiro.test(kdata[kdata$group == "precision", "value"]) # Norm
x5 <- cor.test(hc, kdata[kdata$group == "precision", "value"], method = "pearson")


# Create a dataframe with correlation estimates and p-values
df_stats <- data.frame(
  cor = c(x1$estimate, x2$estimate, x3$estimate, x4$estimate, x5$estimate),
  p.value = c(x1$p.value, x2$p.value, x3$p.value, x4$p.value, x5$p.value)
)

# Set row names
rownames(df_stats) <- c(colnames(res_long[1:4]),"precision")

# Print the dataframe
df_stats



################################# PLOT ########################################

plot(res_long$avg.comp, res_long$avg.purity)



########################### GSA vs BIN #########################################

# Create an empty list to store results for each sample
acc_list <- list()

# Iterate over samples
for (sample in samples) {
  filename <- paste0("gsa_amber/x700c/", sample, "_bin_metrics.tsv")
  cat("Starting: ", filename, "\n")
  
  df <- read.table(filename, sep = "\t", header = TRUE)
  mockdf <- read.table("mock_genomes.txt", sep = "\t", header = TRUE)
  
  df <- df[, c("genome_id", "Tool")]
  
  df$sample <- rep(sample, nrow(df))
  
  # Append the results for the current sample to the list
  acc_list[[sample]] <- df
}

# Combine results for all samples into a single dataframe
acc_long <- do.call(rbind, acc_list)

# Reset row names
rownames(acc_long) <- NULL

# generate a df holding accuracy
acc_df <- data.frame()

for (sample in samples) {
  
  
  bin_count <- length(unique(acc_long[acc_long$Tool != "Gold standard" & acc_long$sample == sample,]$genome_id))

  gsa_count <- length(unique(acc_long[acc_long$Tool == "Gold standard" & acc_long$sample == sample,]$genome_id))

  accuracy <- bin_count/gsa_count
  
  acc_df <- rbind(acc_df, data.frame(sample = sample, 
                                           gsa.count = gsa_count, 
                                           bin.count = bin_count, 
                                           accuracy = accuracy))
  
}
acc_df

##################### PLOT ######################################
group_colors <- brewer.pal(6,"Set2")
#long_acc_df$group <- factor(long_acc_df$group, levels = c(""))


ggplot(long_acc_df[long_acc_df$group != "accuracy",], 
       aes(fill=group, y=count, x=sample)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_manual(values=group_colors[1:2]) +
  labs(subtitle = "Proportion of Host and Endophytic bins ", y = "Fraction", 
       x = "sample") +
  theme(legend.position = "none", legend.title=element_blank())

# THE DATA FINALLY CAVED, THIS IS BIIIIIIIG
# The proportion of bins goes down yes, but the purity increases.
# Indicating, that removal of bins without a classifier might likely be host.
plot(hc, long_acc_df[long_acc_df$group == "accuracy",]$count,
     pch=20, cex = 2)
plot(hc, res_long$avg.purity, pch=20, cex = 2)

# Melt the dataframe to long format and customize the variable names
long_acc_df <- melt(acc_df, id.vars = "sample", variable.name = "group", value.name = "count")



########################### The Preciscion vs Purity plot ######################
group_colors <- brewer.pal(8,"Set2")

# Building a fake plot to steal the legend from.
fake <- data.frame(group = c("endophyte", "host", "OMF", "rfungi", 
                             "ba_ar", "pl_vi_unk", "precision", "purity"),
                   val = seq(1,8,by=1))
fake$group <- factor(fake$group, levels = fake$group)
fake_p <- ggplot(fake, aes(x=val, y=group, fill = group)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=group_colors)


# The KEY DATA
kdata <- data.frame(value = c(long_acc_df[long_acc_df$group == "accuracy",]$count,
                              res_long$avg.purity),
                    hc = hc,
                    group = c(rep("precision", 11), rep("purity", 11)))

kdata$group <- factor(kdata$group, levels = c("precision", "purity"))


the_lp <- ggplot(kdata, aes(x=hc, y = value, color = group)) +
  geom_line(size = 2) +
  theme_bw() +
  scale_color_manual(values=group_colors[7:8]) +
  labs(subtitle = "Purity and Precision vs HC", x = "hc", y = "Fraction")
       

legend <- as_ggplot(get_legend(fake_p))


ggarrange(pp_group,pp_bar_eh, the_lp, legend, labels = c("A", "B", "C"),
          ncol = 2, nrow = 2, common.legend = T, legend="none")


######################## PLAYGRND# #################

acc_df
res_long
####



res_longer <- data.frame(value = c(res_long$avg.comp, res_long$avg.purity, res_long$ARI,
                                    precision = kdata[kdata$group == "precision",]$value),
                         metric = c(rep("avg.comp",11), rep("avg.purity", 11),
                                     rep("ARI",11), rep("precision",11)))

ggplot(res_longer, aes(x=metric, y = value)) +
  geom_violin(aes(fill = metric)) +
  geom_boxplot(width=0.1, color="black", fill = "white", alpha=0.7) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  labs(subtitle = "Distribution of the accumulated results from each sample", y = "Fraction", x = "") +
  theme(legend.position = "none", axis.text.x = element_text(size = 13))
  

res_long$avg.purity
?wilcox.test

# Calculate the differences
differences <- rep(1,11) - res_long$avg.purity

res <- SignTest(differences)


# Count the signs
positive_signs <- sum(differences > 0)
negative_signs <- sum(differences < 0)

# Perform the Sign test (using binomial test)
p_value <- binom.test(min(positive_signs, negative_signs), positive_signs + negative_signs)$p.value

# Print the results
cat("Sign test results:\n")
cat("Positive signs:", positive_signs, "\n")
cat("Negative signs:", negative_signs, "\n")
cat("p-value:", p_value, "\n")
