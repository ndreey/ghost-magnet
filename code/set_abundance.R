host
1-host

SPLIT BETWEEN ENDOIPHYTES
OMF1 1
OMF2 1
OMF3 1



rfungi 31
PRJ18505*0.5size + PRJ122*6.5size.+.+.= 29.7%



bark   35 
plasm  30

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















