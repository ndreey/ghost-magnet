# VERSION 2 which ACCOUNTS FOR GENOME SIZE
# 


# The amount of basepairs CAMISIM is set to generate
seq_target<- 5e9

# Host genome soze
host_size <- 4186550321

# Host contamination level
hc <- 0.5

# HC is 2.5Gb and we therefore need little more than half of host.
adj_host <- (seq_target*hc)/(host_size)


# Calculate for rfungi

adj_factor <- 1 - hc

# Endophytic seq size
endo_seq <- seq_size * adj_factor

# Proportion for group
grp_proportion <- 0.226

# Expected proportion of reads for group 
group_seq <- endo_seq * grp_proportion

# Get random values for the group size
group_size <- length(mock_df[mock_df$group == "rfungi",1])
rndm_abun <- runif(group_size, 0, 1)

# Normalize the relative abundances to fit the adjusted expected proportion of
# the community 
nrm_abu <- grp_proportion * rndm_abun / sum(rndm_abun)

# Thus the first genome should cover 0.00308 of the group_seq. We therefore 
# need to adjust each genome by genome size to cover their proportion
genome_sizes <- mock_df[mock_df$group == "rfungi",5]

counter <- 1
adj_nrm_abu <- vector()
for (size in genome_sizes) {
  
  adj_abu <- group_seq * nrm_abu[counter] / size
  adj_nrm_abu <- c(adj_nrm_abu, adj_abu)
  counter <- counter + 1
}









# set the seed for reproducibility
set.seed(13371337)

set_abundance2 <- function(host_ab, community_size, proportion, genome_size) {
  # Generates relative abundance adjusted for host-contamination
  # and genome size
  # 
  # The expected proportion of the community is adjusted to account for the 
  # proportion of host-contamination and the genome size. Then, a set of 
  # abundance values is generated based on the community size and the adjusted 
  # proportion. These abundance values are such that the sum of abundances for 
  # each community member will add up to the adjusted proportion of the 
  # community.
  # 
  # Args:
  #   host_ab: The proportion of host-contamination in decimal form (float).
  #   community_size: The size of the community (integer).
  #   proportion: The expected proportion of the community (float).
  #   genome_size: Numeric vector with genome sizes in base pairs.
  # 
  # Returns:
  #   Vector with float values representing the adjusted relative abundance
  
  # The adjust factor
  adj_factor <- 1 - host_ab
  
  # Adjust the proportion of the community to fit the host contamination
  # and genome size
  proportion <- proportion * adj_factor * genome_size / sum(proportion * genome_size)
  
  # Generate a vector of n random numbers between 0 and 1 that represent the 
  # relative abundance in the community
  rndm_abun <- runif(community_size, 0, 1)
  
  # Normalize the relative abundances to fit the adjusted expected proportion of
  # the community 
  nrm_abu <- proportion * rndm_abun / sum(rndm_abun)
  
  # Returns the randomly chosen adjusted relative abundances of each genome in
  # the community
  return(nrm_abu)
}



# Load in mock_df
mock_df <- read.table("submission/mock_genomes.txt", 
                      header = TRUE, sep = "\t")

# proportion of endophytes
proportion <- list(rfungi = 0.226, OMF = 0.452, ba_ar = 0.297, pl_vi_unk = 0.025)

# proportion of host-contamination
host_contamination <- c(0, 0.5, 0.8, 0.95, 0.98)

# vector with genome sizes
genome_size <- mock_df$size

for (host_abu in host_contamination){
  
  # The abundance for host_abu level of contamination.
  relative_abundance <- c(host_abu)
  for (group in names(proportion)){
    
    # Gets the subset of mock_df
    group_set <- mock_df[mock_df$group == group,1]
    
    # Gets abundances
    group_abu <- set_abundance(host_abu, length(group_set), proportion[[group]], genome_size[mock_df$group == group])
    
    # Stores it in vector that holds all abundances for the community.
    relative_abundance <- c(relative_abundance, group_abu)
    
  }
  
  # Store genome_id and relative_abundance in one data frame
  abu_df <- data.frame(mock_df[,1], relative_abundance )
  
  # Generates unique filenames
  filename <- paste0(host_abu,"_hc_abundance_v2.tsv")
  
  # Write to .tsv
  write.table(abu_df, file = filename, sep = "\t", row.names = FALSE,
              col.names = FALSE)
  
  cat("The sum of abundance is: ", sum(abu_df$relative_abundance), "\n")
}


