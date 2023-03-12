# Create empty lists to store data
list_genome_id <- vector()
list_size <- vector()
list_NCBI <- vector()
list_tax_group <- vector()
list_group <- vector()

# Set directory containing fasta files
fasta_dir <- "data/references/genomes/source_genomes/"

# Read mapping file to get genome IDs
genome_to_id <- read.table("genome_to_id.txt", header = TRUE, sep = "\t")

# Read report file to get genome sizes
report_genome <- read.table("report_genome.txt", header = TRUE, sep = "\t")

# Read metadata file to get NCBI IDs
metadata <- read.table("metadata.tsv", header = TRUE, sep = "\t")

# Read taxonomic profile file to get taxonomic group IDs
taxonomic_profile <- read.table("taxonomic_profile.tsv", header = TRUE, sep = "\t")

# Loop over fasta files in the directory
for (fasta_file in dir(fasta_dir, pattern = "\\.fasta$")) {
  
  # Get genome ID
  genome_id <- genome_to_id$genome_id[genome_to_id$filename == fasta_file]
  list_genome_id <- c(list_genome_id, genome_id)
  
  # Get size
  size <- report_genome$size[report_genome$filename == fasta_file]
  list_size <- c(list_size, size)
  
  # Get NCBI ID
  NCBI_ID <- metadata$NCBI_ID[metadata$genome_id == genome_id]
  list_NCBI <- c(list_NCBI, NCBI_ID)
  
  # Get taxonomic group ID
  tax_path <- taxonomic_profile$TAXPATH[taxonomic_profile$NCBI_ID == NCBI_ID]
  tax_group <- sub("\\|.*", "", tax_path) # extract first number before "|"
  list_tax_group <- c(list_tax_group, tax_group)
  
  # Get humanized group name
  if (tax_group != "2759") {
    group <- "not_euk"
  } else {
    group <- "euk"
  }
  list_group <- c(list_group, group)
}

# Create dataframe
df <- data.frame(genome_id = list_genome_id, size = list_size, 
                 taxid = list_NCBI, tax_group = list_tax_group, 
                 group = list_group)
