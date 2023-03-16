source("code/set_abundance.R")

# set the seed for reproducibility
set.seed(13371337)

# Load in mock_df
mock_df <- read.table("submission/mock_genomes.txt", 
                        header = TRUE, sep = "\t")

# proportion of endophytes
proportion <- list(rfungi = 0.226, OMF = 0.452, ba_ar = 0.297, pl_vi_unk = 0.025)

# proportion of host-contamination
host_contamination <- c(0, 0.5, 0.8, 0.95, 0.98)

for (host_abu in host_contamination){
  
  # The abundance for host_abu level of contamination.
  relative_abundance <- c(host_abu)
  for (group in names(proportion)){
    
    # Gets the subset of mock_df
    group_set <- mock_df[mock_df$group == group,1]
    
    # Gets abundances
    group_abu <- set_abundance(host_abu, length(group_set), proportion[[group]])
    
    # Stores it in vector that holds all abundances for the community.
    relative_abundance <- c(relative_abundance, group_abu)
    
  }
  
  # Store genome_id and relative_abundance in one data frame
  abu_df <- data.frame(mock_df[,1], relative_abundance )
  
  # Generates unique filenames
  filename <- paste0(host_abu,"_hc_abundance.tsv")
  
  # Write to .tsv
  write.table(abu_df, file = filename, sep = "\t", row.names = FALSE,
              col.names = FALSE)
  
  cat("The sum of abundance is: ", sum(abu_df$relative_abundance), "\n")
}



