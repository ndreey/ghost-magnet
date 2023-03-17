# set the seed for reproducibility
set.seed(13371337)

set_abundance <- function(host_ab, community_size, proportion) {
  # Generates relative abundance adjusted for host-contamination
  # 
  # The expected proportion of the community is adjusted to account for the 
  # proportion of host-contamination. Then, a set of abundance values is 
  # generated based on the community size and the adjusted proportion. These 
  # abundance values are such that the sum of abundances for each community 
  # member will add up to the adjusted proportion of the community.
  # 
  # Args:
  #   host_ab: The proportion of host-contamination in decimal form (float).
  #   n_community: The size of the community (integer).
  #   proportion: The expected proportion of the community (float).
  # 
  # Returns:
  #   Vector with float values representing the adjusted relative abundance
  
  # The adjust factor
  adj_factor <- 1 - host_ab
  
  # Adjust the proportion of the community to fit the host contamination
  proportion <- proportion * adj_factor
  
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











