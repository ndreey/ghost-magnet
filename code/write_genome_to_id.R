# Load in mock_df
mock_df <- read.table("submission/mock_genomes.txt", 
                      header = TRUE, sep = "\t")

genomes <- read.table("tes22t.txt", 
                      header = FALSE, sep = "\t")


path <- "/mnt/c/Users/andbo/thesis_andbo/01_clean_mock_data/data/references/genomes/source_genomes/"

list_path <- vector()
for (id in mock_df$genome_id) {
  
  filename <- grep(paste0(".*", id, ".*\\.(fasta|contigs\\.fasta)(?!\\.fai)$"), genomes$V1, value = TRUE, perl = TRUE)
  
  
  path_name <- paste0(path,filename)
  list_path <- c(list_path, path_name)
}

genomes_df <- data.frame(genome_id = mock_df$genome_id, file_path = list_path)

write.table(genomes_df, file = "01genome_to_id.tsv", sep = "\t", row.names = FALSE, 
            col.names = FALSE)
