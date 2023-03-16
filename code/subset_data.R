################################################################################
######      Subset source_genomes from CAMI2 rhizosphere mock data        ######
################################################################################

############################# FUNCTIONS ########################################
get_size <- function(genome_ids,size_report) {
  # Gets the size of the genome
  # 
  # Args:
  #   genome_ids: Vector of genome ids
  #   size_report: report_genome_size.tsv dataframe
  # 
  # Returns:
  #   Numeric vector with the sizes.
  sizes <- vector()
  for (id in genome_ids) {
    sizes <- c(sizes, size_report[grep(id, size_report$genome), 2])
  }
  return(sizes)
}


############################## UTILITIES #######################################
# set the seed for reproducibility
set.seed(13371337)

# report file holding sizes of genomes
genome_size <- read.table("camisim_setup_files/report_genome_sizes.tsv",
                          header = TRUE, sep = "\t")



########################## FUNGAL GENOMES ######################################
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

# Give each fungal genome the rank fungi
df_fungi$rank <- rep("fungi", nrow(df_fungi))

# Group the fungi to rhizosphere fungi (rfungi) and OMF
df_fungi$group <- c(rep("rfungi", 31), rep("OMF",3))


# Add the genome sizes
fungi_size <- get_size(df_fungi$genome_id, genome_size)
df_fungi$size <- fungi_size



########################### PLASMID/VIRUS/UNKNOWN ##############################
# 399 genomes/contigs with four columns, we only need the first two cols which
# represents filename and taxid.

# Read in plasmids.tsv
plasmids <- read.table("camisim_setup_files/plasmids.tsv", 
                       header = FALSE, sep = "\t")
# Remove third column as it is equal to second.
plasmids <- plasmids[,-3]
colnames(plasmids) <- c("genome_id", "taxid", "rank")

# We dont need 399 of these genomes. Lets reduce it. 399 - 369 = 30
# Subsets 369 genomes randomly
plasmid_subset <- sample(plasmids$genome_id, 369)

# Subsets the 30 plasmids
df_plasmids <- plasmids[!(plasmids$genome_id %in% plasmid_subset),]

# Remove the ".fasta" suffix
df_plasmids$genome_id <- gsub(".fasta", "", df_plasmids$genome_id)

# Group the genomes 
df_plasmids$group <- rep("pl_vi_unk", nrow(df_plasmids))

# Add the genome sizes
plasmid_size <- get_size(df_plasmids$genome_id, genome_size)
df_plasmids$size <- plasmid_size



######################## BACTERIA / ARCHAEA ####################################
# I want to subset 35 genomes from the ~400 genomes 
# 
# reference.tsv holds info on 417 bacteria/archaea genomes.
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
rownames(df_bacteria) <- NULL


# Lets check what ranks each taxid has by using the taxonomic profile.
# taxonomic profile file to get taxonomic group IDs
taxonomic_profile <- read.table("camisim_setup_files/taxonomic_profile_1.txt", 
                                header = FALSE, sep = "\t", skip = 5)

archaea_taxids <- vector()
# Lets see if the taxid matches any taxid where the row has "Archaea"
for (taxid in df_bacteria$taxid) {
  get_match <- grep(taxid,taxonomic_profile$V1[grep("Archaea",
                                                    taxonomic_profile$V4)])
  
  # When grep find nothing it returns interger(0). Thus, any value where length
  # is not 0 means that there has been a match
  if (length(get_match) != 0) {
    archaea_taxids <- c(archaea_taxids, taxid)
  }
}

# Set each taxids rank
df_bacteria$rank <- ifelse(df_bacteria$taxid %in% archaea_taxids, "archaea",
                           "bacteria")

# Lets group this community
df_bacteria$group <- rep("ba_ar", nrow(df_bacteria))

# Structure the data frame
df_bacteria <- df_bacteria[, c("ref", "taxid", "rank", "group")]
colnames(df_bacteria)[1] <- "genome_id" 

# Add the genome sizes
bact_size <- get_size(df_bacteria$genome_id, genome_size)
df_bacteria$size <- bact_size




######################## BUILD MOCK DATAFRAME ###################################

# Host
df_host <- data.frame(genome_id = "Platanthera_zijinensis_chr", taxid = 2320716,
                      rank = "orchid", group = "host", size = 4186550321)

# Bind together the subsets
mock_df <- rbind(df_host, df_fungi, df_bacteria, df_plasmids)

# Create new rownames
rownames(mock_df) <- NULL

# Final structure
mock_df <- mock_df[, c("genome_id", "rank", "taxid", "group", "size")]

# Write to txt.
#write.table(mock_df, file = "mock_genomes.txt", sep = "\t", row.names = FALSE, 
#            col.names = TRUE)



