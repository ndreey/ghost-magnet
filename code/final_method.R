library(ggplot2)
library(RColorBrewer)
library(ggpubr)
source("analysis_tools.R")


######################## THE ANALYSIS OF METHOD ################################

samples <- c("X00", "X01", "X02", "X03", "X04", "X05", "X06", "X07",
             "X08", "X090", "X095")
hc <- c(0, 0.25, 0.43, 0.57, 0.66, 0.76, 0.81, 0.88, 0.92, 0.96, 0.98)

# Data files were generated with eval_maxbin_bin.R
# Here we aquire the long formated table of the output
df_long <- read.table("binning/df_long_mockdb.tsv", sep = "\t", header = TRUE)

# Factor so colors match with other plots
df_long$group <- factor(df_long$group, levels = c("host", "OMF", "rfungi", 
                                                   "ba_ar", "pl_vi_unk"))


# Keeps the counts for host bins across replis
host_long <- df_long[df_long$group == "host",]


## Generate the data frame holding proportions
# Subset of the data frame
subset_df <- df_long[df_long$group != "host",]

# Sum counts by sample
sum_counts <- summarise(group_by(subset_df, repli, sample), 
                        total_counts = sum(counts))
colnames(sum_counts) <- c("repli", "sample", "counts")

# Build to counts to long format
count_long <- data.frame(repli = c(sum_counts$repli, host_long$repli), 
                         sample = c(sum_counts$sample, host_long$sample),
                         counts = c(sum_counts$counts, host_long$counts),
                         group = c(rep("endophyte",nrow(sum_counts)),
                                   as.character(host_long$group)),
                         hc = c(hc, hc))
# Factor accordingly
count_long$repli <- factor(count_long$repli,
                           levels = c("raw", "x1000c", "x700c", "maxbin"))
count_long$group <- factor(count_long$group, levels = c("endophyte", "host"))


# Empty df to store proportions
prop_df <- data.frame(repli = character(),
                      prop = numeric(),
                      n_bin = numeric(),
                      n_endo = numeric(),
                      sample = character(),
                      h.prop = numeric(),
                      stringsAsFactors = FALSE)

# Loop through to generate proportions
for (r in c("raw", "x1000c", "x700c", "maxbin")) {
  cat("Repli:", r, "\n")
  
  for (s in samples) {
    cat("    Sample:", s, "\n")
    
    n_bin <- sum(count_long[count_long$repli == r & 
                              count_long$sample == s, "counts"])
    
    n_endo <- sum(count_long[count_long$repli == r &
                               count_long$sample == s &
                               count_long$group == "endophyte", "counts"])
    prop <- n_endo / n_bin
    
    h_prop <- 1 - prop
    
    # Append df to existing
    prop_df <- rbind(prop_df, data.frame(repli = r,
                                         prop = prop,
                                         n_bin = n_bin,
                                         n_endo = n_endo,
                                         sample = s,
                                         h.prop = h_prop))
  }
}

# Factor and adding the hc levels
prop_df$repli <- factor(prop_df$repli, levels = c("raw", "x1000c", "x700c",
                                                  "maxbin"))
prop_df$hc <- rep(hc, 4)

df_long$group

########################## PLOT ARENA ##########################################
# Color vector to subset from
group_colors <- brewer.pal(8,"Set2")
#display.brewer.all()
df_long$group


df_long$hc <- NA
df_long[df_long$sample == "X00","hc"] <- hc[1]*100


# Generate stacked bar plots for the proportion of mock groups
# Plot of percentage stack plot
# raw
raw_df <- df_long[df_long$repli == "raw",]
pp_raw <- ggplot(raw_df, aes(fill = group, y = counts, x = sample)) + 
  geom_bar(position = "fill", stat = "identity",  color = "black") +
  scale_fill_manual(values = group_colors[2:6]) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  labs(title = "raw", y = "Fraction", x = "") 

# x1000c
x1000c_df <- df_long[df_long$repli == "x1000c",]  
pp_x1000c <- ggplot(x1000c_df, aes(fill = group, y = counts, x = sample)) + 
  geom_bar(position = "fill", stat = "identity",  color = "black") +
  scale_fill_manual(values = group_colors[2:6]) +
  theme_bw() +
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), 
        axis.ticks.x=element_blank(), axis.text.x=element_blank(),
        legend.position = "none", legend.title=element_blank()) +
  labs(subtitle = "x1000c", y = "", x = "") 

# x700c    
x700c_df <- df_long[df_long$repli == "x700c",]
pp_x700c <- ggplot(x700c_df, aes(fill = group, y = counts, x = sample)) + 
  geom_bar(position = "fill", stat = "identity", color = "black") +
  scale_fill_manual(values = group_colors[2:6]) +
  theme_bw() +
  theme(legend.position = "none", legend.title=element_blank(),
        axis.ticks.y=element_blank(), axis.text.y=element_blank(), 
        axis.ticks.x=element_blank(), axis.text.x=element_blank()) +
  labs(title = "x700c", y = "", x = "") 
  
# maxbin
maxbin_df <- df_long[df_long$repli == "maxbin",]
pp_max <- ggplot(maxbin_df, aes(fill = group, y = counts, x = sample)) + 
  geom_bar(position = "fill", stat = "identity", color = "black") +
  scale_fill_manual(values=group_colors[2:6]) +
  theme_bw() +
  theme(legend.position = "none", legend.title=element_blank(),
        axis.ticks.y=element_blank(), axis.text.y=element_blank(), 
        axis.ticks.x=element_blank(), axis.text.x=element_blank()) + 
  labs(title = "maxbin", y = "", x = "") 

# Arranging them accordingly
ggarrange(pp_raw, pp_x1000c, pp_x700c, pp_max, labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")





#Generate stacked bar plots for the proportion of mock groups
# Plot of percentage stack plot
# raw
raw_df <- df_long[df_long$repli == "raw",]
pp_raw2 <- ggplot(raw_df, aes(fill = group, y = counts, x = as.character(hc))) + 
  geom_bar(position = "dodge", stat = "identity",  color = "black") +
  scale_fill_manual(values = group_colors[2:6]) +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "none") +
  labs(subtitle = "", y = "Counts", x = "Host Contamination (%)") +
  ylim(0,55)

# x1000c
x1000c_df <- df_long[df_long$repli == "x1000c",]  
pp_x1000c2 <- ggplot(x1000c_df, aes(fill = group, y = counts, x = sample)) + 
  geom_bar(position = "dodge", stat = "identity",  color = "black") +
  scale_fill_manual(values = group_colors[2:6]) +
  theme_bw() +
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(),
        legend.position = "none", legend.title=element_blank(),
        axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
  labs(subtitle = "x1000c", y = "", x = "") +
  ylim(0,55)

# x700c    
x700c_df <- df_long[df_long$repli == "x700c",]
pp_x700c2 <- ggplot(x700c_df, aes(fill = group, y = counts, x = as.character(hc))) + 
  geom_bar(position = "dodge", stat = "identity", color = "black") +
  scale_fill_manual(values = group_colors[2:6]) +
  theme_bw() +
  theme(legend.position = "none", legend.title=element_blank(),
        axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
  labs(subtitle = "", y = "", x = "Host Contamination (%)") +
  ylim(0,55)

# maxbin
maxbin_df2 <- df_long[df_long$repli == "maxbin",]
pp_max2 <- ggplot(maxbin_df, aes(fill = group, y = counts, x = as.character(hc))) + 
  geom_bar(position = "dodge", stat = "identity", color = "black") +
  scale_fill_manual(values=group_colors[2:6]) +
  theme_bw() +
  theme(legend.position = "none", legend.title=element_blank(),
        axis.ticks.y=element_blank(), axis.text.y=element_blank())+ 
  labs(subtitle = "", y = "", x = "Host Contamination (%)") +
  ylim(0,55)

# Arranging them accordingly
ggarrange(pp_raw2, pp_x1000c2, pp_x700c2, pp_max2, labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")

ggarrange(pp_raw, pp_x700c, pp_max, pp_raw2, pp_x700c2, pp_max2,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom")

################## RADAR #####################


# Library
library(fmsb)

raw_radar <- aggregate(counts ~ group, data = df_long[df_long$repli == "raw", ], FUN = sum)
# Reshape the summarized dataframe from long to wide format
wide_df1 <- pivot_wider(raw_radar, names_from = group, values_from = counts)
wide_df1 <- rbind(rep(400,5), rep(0,5),wide_df)
radarchart(wide_df1)

# Summarize count values for each group
x1000c_radar <- aggregate(counts ~ group, data = df_long[df_long$repli == "x1000c", ], FUN = sum)
wide_df2 <- pivot_wider(x1000c_radar, names_from = group, values_from = counts)
wide_df2 <- rbind(rep(300,5), rep(0,5),wide_df2)
radarchart(wide_df2)

# Summarize count values for each group
x700c_radar <- aggregate(counts ~ group, data = df_long[df_long$repli == "x700c", ], FUN = sum)
wide_df3 <- pivot_wider(x700c_radar, names_from = group, values_from = counts)
wide_df3 <- rbind(rep(200,5), rep(0,5),wide_df3)
radarchart(wide_df3)
# Summarize count values for each group
max_radar <- aggregate(counts ~ group, data = df_long[df_long$repli == "maxbin", ], FUN = sum)
wide_df4 <- pivot_wider(max_radar, names_from = group, values_from = counts)
wide_df4 <- rbind(rep(150,5), rep(0,5),wide_df4)
radarchart(wide_df4)


# Generate barplots showing counts how each repli differed in counts
# Raw
bp_raw <- ggplot(count_long[count_long$repli == "raw",], 
                 aes(fill=group, y=counts, x=sample)) + 
  geom_bar(position="dodge", stat="identity",  color = "black") +
  scale_fill_manual(values = group_colors[1:2]) +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.position = "none", legend.title=element_blank()) +
  labs(subtitle = "raw", y = "Counts", x = "") 

# x1000c
bp_x1000c <- ggplot(count_long[count_long$repli == "x1000c",],
                    aes(fill=group, y=counts, x=sample)) + 
  geom_bar(position="dodge", stat="identity",  color = "black") +
  scale_fill_manual(values = group_colors[1:2]) +
  theme_bw() +
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(), 
        axis.ticks.x=element_blank(), axis.text.x=element_blank(),
        legend.position = "none", legend.title=element_blank()) +
  labs(subtitle = "x1000c", y = "", x = "") 

# x700c
bp_x700c <- ggplot(count_long[count_long$repli == "x700c",], 
                   aes(fill=group, y=counts, x=sample)) + 
  geom_bar(position="dodge", stat="identity", color = "black") +
  scale_fill_manual(values = group_colors[1:2]) +
  theme_bw() +
  theme(legend.position = "none", legend.title=element_blank()) +
  labs(subtitle = "x700c", y = "Counts", x = "sample") 
  
# maxbin
bp_max <- ggplot(count_long[count_long$repli == "maxbin",], aes(fill=group, y=counts, x=sample)) + 
  geom_bar(position="dodge", stat="identity", color = "black") +
  scale_fill_manual(values = group_colors[1:2]) +
  theme_bw() +
  theme(legend.position = "none", legend.title=element_blank(),
        axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
  labs(subtitle = "maxbin", y = "", x = "sample") 

# Arrangement
ggarrange(bp_raw, bp_x1000c, bp_x700c, bp_max, labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")


# Generate barplots showing how the proportion of endophytes differed between
# raw, x1000c and x700c
# Proportion barplot


pp2 <- ggplot(prop_df[prop_df$repli != "maxbin",],
       aes(x=sample, y = prop, fill = repli)) +
  geom_bar(position = "fill", stat = "identity", color = "black") +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  labs(subtitle = "raw", y = "Fraction", x = "") 

ggarrange(pp2, pp1, labels = c("A", "B"),
            ncol = 1, nrow = 2, common.legend = TRUE, legend="bottom")  


# Violin
# Violin plot
violino <- ggplot(prop_df, aes(x=repli, y=prop, fill = repli)) + 
  geom_violin() +
  geom_boxplot(width=0.1, color="black", fill = "lightgrey", alpha=1) +
  scale_fill_manual(values=group_colors[1:4]) +
  theme_bw() + 
  theme(legend.position = "none", legend.title = element_blank()) +
  labs(title = "", x = "replicate", y = "Proportion of Endophytic Bins") 
   



pp1 <- ggplot(prop_df[prop_df$repli != "maxbin",], 
              aes(x=sample, y = prop, fill = repli)) +
  geom_bar(position="dodge", stat = "identity", color = "black") +
  theme_bw() +
  scale_fill_manual(values=group_colors[1:4]) +
  labs(subtitle = "", 
       y = "Proportion of Endophytic Bins", x = "sample") +
  theme(legend.position = "bottom", legend.title=element_blank())


ggarrange(violino, pp1,labels = c("A", "B"),
          ncol = 1, nrow = 2, 
          common.legend = TRUE, 
          legend="top")

  
# Mnjaaa
ggplot(prop_df, aes(x=sample, y=prop, shape=repli,color=repli)) + 
  geom_point(size=3) +
  theme_bw()


################### Analysis ######################


shapiro.test(prop_df[prop_df$repli == "raw",2]) # p-value = 0.003279 Not normal
cor.test(hc, prop_df[prop_df$repli == "raw",2], method = "kendall", exact = F)
# p-value = 0.0001574 r = -0.8975491 
median(prop_df[prop_df$repli == "raw",2])
sd(prop_df[prop_df$repli == "raw",2])


shapiro.test(prop_df[prop_df$repli == "x1000c",2]) # p-value = 0.007389 Not normal
cor.test(hc, prop_df[prop_df$repli == "x1000c",2], method = "kendall", exact = F)
# p-value = 8.269e-05 r = -0.934947  
median(prop_df[prop_df$repli == "x1000c",2])
sd(prop_df[prop_df$repli == "x1000c",2])

shapiro.test(prop_df[prop_df$repli == "x700c",2]) # p-value  = 0.06005 normal
cor.test(hc, prop_df[prop_df$repli == "x700c",2], method = "kendall", exact = F)
# p-value = 4.241e-05 r = -0.9723449 
median(prop_df[prop_df$repli == "x700c",2])
sd(prop_df[prop_df$repli == "x700c",2])
# none as it si 1.
shapiro.test(prop_df[prop_df$repli == "maxbin",2]) # 
cor.test(hc, prop_df[prop_df$repli == "maxbin",2], method = "kendall", exact = F)
median(prop_df[prop_df$repli == "maxbin",2])
sd(prop_df[prop_df$repli == "maxbin",2])

kruskal.test(prop~repli,prop_df) 
# Kruskal-Wallis chi-squared = 21.09, df = 3, p-value = 0.0001008

wilcox.test(prop_df[prop_df$repli == "raw",2], 
            prop_df[prop_df$repli == "x700c",2],
            exact = F)

wilcox.test(prop_df[prop_df$repli == "x1000c",2],
            prop_df[prop_df$repli == "x700c",2],
            exact = F)


prop_df



pairwise.wilcox.test(prop_df$prop, prop_df$repli, exact = F,
                     p.adjust = "BH")


# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  prop_df$prop and prop_df$repli 
# 
#          raw     x1000c  x700c  
#   x1000c 0.86823 -       -      
#   x700c  0.82568 0.82871 -      
#   maxbin 0.00018 0.00018 0.00018
# 
# P value adjustment method: BH 


write.csv(prop_df, "prop_df.csv", row.names = FALSE)




round(prop_df[,-c(1,5)],3)


