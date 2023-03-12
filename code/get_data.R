# Create empty lists to store data
list_genome_id <- vector()
list_size <- vector()
list_taxid <- vector()
list_tax_group <- vector()
list_group <- vector()
fix_list <- vector()

# Path for the directory containing the genomes
fasta_dir <- "data/references/genomes/source_genomes/"

# file to get genome IDs from fasta name
genome_to_id <- read.table("camisim_setup_files/genome_to_id.tsv", 
                           header = FALSE, sep = "\t")
colnames(genome_to_id) = c("genome_id", "path")

# report file to get genome sizes through fasta
genome_size <- read.table("camisim_setup_files/report_genome_sizes.tsv",
                            header = TRUE, sep = "\t")

# metadata file to get NCBI IDs
metadata <- read.table("camisim_setup_files/metadata.tsv",
                       header = TRUE, sep = "\t")

# taxonomic profile file to get taxonomic group IDs
taxonomic_profile <- read.table("camisim_setup_files/taxonomic_profile_1.txt", 
                                header = FALSE, sep = "\t", skip = 5)
count = 0
# Loop over fasta files in the directory that matches the regex for .fasta
for (fasta_file in dir(fasta_dir, pattern = "\\.fasta$")) {
  
  # Get genome ID
  # find row in genome_to_id that matches fasta_file
  row_idx <- grep(fasta_file, genome_to_id[,2])
  
  # extract genome ID from matching row
  genome_id <- genome_to_id[row_idx, 1]
  
  list_genome_id <- c(list_genome_id, genome_id)
   
  # Get size by matching genome name
  size <- genome_size$size[genome_size$genome == fasta_file]
  list_size <- c(list_size, size)
   
  # Get NCBI ID
  taxid <- metadata$NCBI_ID[metadata$genome_ID == genome_id]
  list_taxid <- c(list_taxid, taxid)
   
  # Get taxonomic group ID, V3 = TAXPATH, V1 = TAXID
  # Find the row where V1 column matches taxid exactly
  matching_row <- grep(paste0("^", taxid, "$"), taxonomic_profile$V1)
  
  if (length(matching_row) == 0) {
    # add taxid to fix_list and move to next iteration
    fix_list <- c(fix_list, taxid)
    list_tax_group <- c(list_tax_group, "fix")
    list_group <- c(list_group, "fix")
    next
  }
  
  # get the tax_path value from the matching row
  tax_path <- taxonomic_profile$V3[matching_row]
  
  # extract the first occuring number from tax_path
  tax_group <- sub(".*?(\\d+).*", "\\1", tax_path)
  
  # add the tax_group to a list
  list_tax_group <- c(list_tax_group, tax_group)
  
  # Get humanized group name
  if (tax_group != 2759) {
    group <- "not_euk"
  } else {
    group <- "euk"
  }
  # add group to a list
  list_group <- c(list_group, group)

}


# Create dataframe
df <- data.frame(genome_id = list_genome_id, size = list_size, 
                 taxid = list_taxid, tax_group = list_tax_group, 
                 group = list_group)


# Find rows where group is "fix"
fix_rows <- which(df$group == "fix")

# Loop over rows and prompt for tax group
for (i in fix_rows) {
  taxid <- df$taxid[i]
  tax_group <- df$tax_group[i]
  new_group <- readline(paste0("Taxid ", taxid, " currently has group ", 
                               tax_group, ". What should the new group be? "))
  
  # Replace "fix" with new group
  df$tax_group[i] <- new_group
  df$group[i] <- ifelse(new_group == 2759, "euk", "not_euk")
}
