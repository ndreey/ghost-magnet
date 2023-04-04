# Second version of code to generate abundance values for CAMISIM de novo 
# community design. It combines set_abundance.R and write_abundance_tsv.R + it
# adjusts for each genomes genome size.


# Load in mock_df
mock_df <- read.table("submission/mock_genomes.txt", 
                      header = TRUE, sep = "\t")

# The amount of bp CAMISIM is set to generate
seq_target<- 5e9

# Host genome soze
host_size <- 4186550321

# proportion of endophytes
proportion <- list(rfungi = 0.226, OMF = 0.452, ba_ar = 0.297, pl_vi_unk = 0.025)

# proportion of host-contamination
host_contamination <- c(0, 0.5, 0.8, 0.95, 0.98)


# Write abundance.tsv files for each HC level.
for (host_abu in host_contamination){
  
  # Adjusts the proportion based on host size
  adj_host <- (seq_target*host_abu) / host_size
  
  # Vector holding relative abundance for all genomes
  relative_abundance <- c(adj_host)
  
  # Calculate abundance values for endophyte group by accounting for hc and each
  # genomes size.
  for (group in names(proportion)){
    
    # Adjustment actor accounting for the abundance of hc
    adj_factor <- 1 - hc
    
    # Endophyte sequence size
    endo_seq <- seq_size * adj_factor
    
    # Gets the subset of mock_df representing group
    group_set <- mock_df[mock_df$group == group,1]
    
    # Expected proportion of bp for group 
    group_seq <- endo_seq * proportion[[group]]
    
    # Get random abundance values for all genomes in group
    rndm_abun <- runif(length(group_set), 0, 1)
    
    # Normalize the abundances to fit the expected proportion of
    # the group
    nrm_abu <- proportion[[group]] * rndm_abun / sum(rndm_abun)
    
    # Vector of the groups genome sizes
    genome_sizes <- mock_df[mock_df$group == group,5]
    
    # Adjust each genome for its genome size and add it to vector
    counter <- 1
    adj_nrm_abu <- vector()
    for (size in genome_sizes) {
      adj_abu <- group_seq * nrm_abu[counter] / size
      adj_nrm_abu <- c(adj_nrm_abu, adj_abu)
      counter <- counter + 1
    }
    
    # Add the adjusted and normalized abundance values to the relative abundance
    # vector
    relative_abundance <- c(relative_abundance, adj_nrm_abu)
    
  }
  
  # Store genome_id and relative_abundance in one data frame
  abu_df <- data.frame(mock_df[,1], relative_abundance )
  
  # Generates unique filenames
  filename <- paste0(host_abu,"_hc_abundance.tsv")
  
  # Write to .tsv
  write.table(abu_df, file = filename, sep = "\t", row.names = FALSE,
              col.names = FALSE)

}