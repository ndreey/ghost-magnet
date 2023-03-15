# Subset source_genomes from CAMI2 rhizosphere mock data
#

# set the seed for reproducibility
set.seed(13371337)


######################## FUNGAL GENOMES ###############################
# Grabs the subset of fungal genomes.
fungal_genomes <- read.table("camisim_setup_files/fungal_genomes.txt",
                             header = FALSE, sep = "\t")
# Lets subset the metadata
fungal_metadata <- read.table("camisim_setup_files/fungal_metadata.tsv",
                            header = FALSE, sep = "\t")
# Keeps the rows that match to fungal_genomes
# Lets also remove the third column as it holds unrelated paths
df_fungi <- fungal_metadata[fungal_metadata$V1 %in% fungal_genomes$V1,-3]

# Fix the colnames
colnames(df_fungi) <- c("genome_id", "taxid")

# Manually adding the three extra OMF genomes.
omf_genomes <- c("Tulasnella_calospora_Tulcal1", "Ceratobasidium_sp_CerAGI",
         "Rhizoctonia_solani_Rhisola1")
omf_taxid <- c(156515, 305860, 1287689)

# Create a data frame from the vectors
omf_df <- data.frame(genome_id = omf_genomes, taxid = omf_taxid)

# Append the new data to the existing data frame
df_fungi <- rbind(df_fungi, omf_df)

# Lets add "OMF" as rank to the omf_genomes and "fungi" to the remaing
df_fungi$rank <- c(rep("fungi", 31), rep("OMF",3))



######################## PLASMID GENOMES ###############################
# 399 genomes/contigs with four columns, we only need the first two cols which
# represents filename and taxid.

# Read in plasmids.tsv
plasmids <- read.table("camisim_setup_files/plasmids.tsv", 
                       header = FALSE, sep = "\t")
# Remove third column as it is equal to second.
plasmids <- plasmids[,-3]
colnames(plasmids) <- c("genome_id", "taxid", "rank")

# We dont need 399 of these genomes. Lets reduce it. 399 - 349 = 50
# Subsets 349 genomes randomly
plasmid_subset <- sample(plasmids$genome_id, 349)

# Subsets the 50 plasmids
df_plasmids <- plasmids[!(plasmids$genome_id %in% plasmid_subset),]

# Remove the ".fasta" suffix
df_plasmids$genome_id <- gsub(".fasta", "", df_plasmids$genome_id)


# Correction: reference.tsv hold both bacterial and archaea!
######################## BACTERIAL GENOMES ###############################
# I want to subset 35 genomes from the ~440 bacterial genomes 
# 
# reference.tsv holds info on 417 bacterial genomes.
# They are grouped as:"known_strain", "new_species", "new_strain", "new_genus",
# "new_order", "new_family"  
# I will take a subset from the 231 "known_strains"

# Load in the data and remove the 3rd column as it holds paths
reference <- read.table("camisim_setup_files/references.tsv", 
                       header = FALSE, sep = "\t")
colnames(reference) <- c("taxid", "ref", "path", "category")
reference <- reference[,-3]

# Randomly subset 35 of the known_strains bioproj reference
bact_subset <- sample(reference[reference$category == "known_strain",]$ref, 35)

# Subset to only keep the 35 chosen and remove third column as it is not needed.
df_bacteria <- reference[(reference$ref %in% bact_subset),-3]

# Sort out df to match the structure of the others.
# Load file to get genome ids from ref
genome_to_id <- read.table("camisim_setup_files/genome_to_id.tsv", 
                           header = FALSE, sep = "\t")
colnames(genome_to_id) = c("genome_id", "path")

# Get the row where ref matches path in genome_to_id and add it to vector.
bact_ids <- vector()

for (ref in df_bacteria$ref) {
  row_idx <- grep(ref, genome_to_id[,2])
  bact_ids <- c(bact_ids, genome_to_id$genome_id[row_idx])
}

# Structure dataframe to resemble the others.
df_bacteria$genome_id <- bact_ids
df_bacteria$rank <- rep("bacteria", 35)
df_bacteria <- df_bacteria[, c("genome_id", "taxid", "rank")]

###                                                                 ###
###   CHECK IF THIS IS TRULY 100% BACTERIA, CAN BE ARCHAEA AS WELL  ###
###                                                                 ###

# taxonomic profile file to get taxonomic group IDs
taxonomic_profile <- read.table("camisim_setup_files/taxonomic_profile_1.txt", 
                                header = FALSE, sep = "\t", skip = 5)

# Lets see if the taxid matches any taxid where the row has "Archaea"
for (taxid in df_bacteria$taxid) {
  get_match <- grep(taxid,taxonomic_profile$V1[grep("Archaea",
                                                    taxonomic_profile$V4)])
  
  # When grep find nothing it returns interger(0). Thus, any value where length
  # is not 0 means that there has been a match
  if (length(get_match) != 0) {
    cat("Taxid:", taxid, "is Archaea")
  }
}

# Taxid: 456442 is Archaea
# Thus, reference.tsv holds bacteria and archaea.
# Replace the rank value with archaea.
df_bacteria$rank <- ifelse(df_bacteria$taxid == 456442, "archaea",
                           df_bacteria$rank)


