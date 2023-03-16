# set the seed for reproducibility
set.seed(13371337)
# Host abundance will vary
host_ab <- 0.5

# The adjust factor
adj_factor <- 1-host_ab

# Ratio of endophytes
ratio <- list(OMF = 0.452, rfungi = 0.226, bark = 0.297, plasm = 0.025)

# Subset 

# Expected amount of unique species.
i <- 30 # minimum length of vector
j <- 40 # maximum length of vector

# generate a random integer between i and j
n <- sample(i:j, size = 1) 

# Generate a random subset of n genomes
set_genomes <- mock_df[mock_df$group == "bark", 1]
n_subset <- sample(set_genomes, n)

# generate a vector of n random numbers between 0 and 1
x <- runif(n, 0, 1) 

# normalize the vector so that its sum is 0.229
nrm_x <- 0.229 * x / sum(x)

# check that the sum of the elements is approximately 0.229
sum(nrm_x)






###############################################################
datasize = 10gb

rfungi = datasize*0.297
nrow(rfungi)
sample(rfungi, rndmNUMBER(15 to nrow(rfungi)), replace = True)
return subset  min 15 max 34


# ALL entries in SUBSET must have value
free = 29.7%
# Protect first 5 iterations from grabbing all of FREE
get random entry
DIS = give random amount between 0 to free
free = free - DIS 

LAST GETS THE REST

###  # # sample numbers that sum(set) ############ think this baby
# Not sort the subset, FIRST element is due to be bigger
for nrow -1
free_relative_data = .297
free_data = datasize * free_relative_data
nrow = random(15,max_row)
random_factor = runif(1, 0, 1)     # random factor have to be betwen 0 and 1 
sample_data = free_data * random_factor
free_data = free_data - sample_data  

sample_data = free_data # ti´s twas thy last data

("Ochroconis": 0.1, "haseolina_M":0.9 ..., nrow)

29.7%


x <- c(1:10)
sample(x)  















