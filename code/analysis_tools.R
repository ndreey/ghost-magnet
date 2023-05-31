library(dplyr)



samples <- c("X00", "X01", "X02", "X03", "X04", "X05", "X06", "X07",
             "X08", "X090", "X095")
hc <- c(0, 0.25, 0.43, 0.57, 0.66, 0.76, 0.81, 0.88, 0.92, 0.96, 0.98)

# As the TSV files follow same structure i will create "quality of life functions"
# This module is to be sourced to each metric analysis

######## ASSEMBLY EVAL ##################
clean_data <- function(df, replicate) {
  
  # Remove ".contigs" suffix from row names
  rownames(df) <- gsub(".contigs", "", rownames(df))
  # Remove numbers and dot prefix from row names
  rownames(df) <- sub("^\\d+\\.", "", rownames(df))
  
  # Add better 
  sample_id <- c("X00", "X01", "X02", "X03", "X04", "X05", "X06", "X07",
                 "X08", "X090", "X095")
  
  colnames(df) <- sample_id
  
  # Order rownames
  df <- df[order(rownames(df)),]
  
  # Add group and order gsa
  mock_df <- read.table("mock_genomes.txt", sep = "\t", header = TRUE, 
                        row.names = 1)
  # Order row names
  mock_df <- mock_df[order(rownames(mock_df)),]
  
  df$group <- mock_df$group
  
  df <- df[order(df$group),]
  
  df$repli <- rep(replicate, nrow(df))
  
  return(df)
}


plot_maxbin_stats <- function(df, title) {
  median <- apply(df, 2, median, na.rm = TRUE)
  variance <- apply(df, 2, var, na.rm = TRUE)
  sd <- apply(df, 2, sd, na.rm = TRUE)
  Q1 <- vector()
  Q3 <- vector()
  IQR <- vector()
  
  for (sample in 1:11) {
    q1_temp <- quantile(df[, sample], 0.25, na.rm = TRUE)
    q2_temp <- quantile(df[, sample], 0.75, na.rm = TRUE)
    iqr_temp <- q2_temp - q1_temp
    
    Q1 <- c(Q1, q1_temp)
    Q3 <- c(Q3, q2_temp)
    IQR <- c(IQR, iqr_temp)
  }
  
  stats <- data.frame(median = median, variance = variance, sd = sd,
                      sample = sample_id, Q1 = Q1, Q3 = Q3, iqr = IQR)
  
  plot(hc, median, main = "mock df\nMedian", xlab = "Host Contamination (%)")
  plot(hc, variance, main = "mock df\nVariance")
  
  plot(hc, stats$median, pch = 16, cex = 1.1, ylab = "Median",
       main = paste("Duplication Ratio of", title, "\nIQR Whiskers"), xlab = "Proportion of Host Contamination")
  arrows(x0 = hc, y0 = stats$Q1, x1 = hc, y1 = stats$Q3,
         code = 3, angle = 90, length = 0.1)
}



######### BINNING EVAL ###########
mock_groupify <- function(df, mockdf) {
  # Empty column
  df$group <- NA
  
  # Loop through each row in df
  for (i in 1:nrow(df)) {
    # Map to the corresponding mockdf genomeid
    matching_row <- mockdf$genome_id == df$genome_id[i]
    df$group[i] <- mockdf$group[matching_row]
  }
  
  return(df)
}




retain_binid <- function(dataframe) {
  dataframe$BINID <- basename(dataframe$BINID)
  dataframe$BINID <- gsub(".fa", "", dataframe$BINID)
  return(dataframe)
}


regex_split <- function(contig) {
  # Decides what group the contig belongs to
  contig[grepl("RNODE", contig)] <- "pl_vi_unk" 
  contig[grepl("PRJ", contig)] <- "ba_ar"
  contig[grepl("tig", contig)] <- "rfungi"
  contig[grepl("scaffold", contig)] <- "OMF"
  contig[grepl("Rhizoctonia", contig)] <- "OMF"
  contig[grepl("Chr", contig)] <- "host"
  return(contig)
}



get_highest_score <- function(data) {
  # Gets the qseqid with highest score.
  # Group the data by qseqid.
  qseqids <- group_by(data, qseqid)
  
  # For each group, keep the row with the highest bitscore
  max_bitscore_rows <- slice(qseqids, which.max(bitscore))
  
  # Remove the grouping information
  filtr_data <- ungroup(max_bitscore_rows)
  
  return(filtr_data)
}

get_abundant <- function(data) {
  # Aquires the most abundant group.
  groups <- regex_split(data$sseqid)
  table <- data.frame(table(groups))
  most_abundant <- table[table$Freq == max(table$Freq),1]
  return(as.character(most_abundant))
}

