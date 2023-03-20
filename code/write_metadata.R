

# Load in mock_df
mock_df <- read.table("submission/mock_genomes.txt", 
                      header = TRUE, sep = "\t")


meta_df <- data.frame(genome_ID = mock_df$genome_id,
                      OTU = mock_df$taxid,
                      NCBI_ID = mock_df$taxid,
                      novelty_category = mock_df$group)



write.table(meta_df, file = "metadata.tsv", sep = "\t", row.names = FALSE, 
                        col.names = TRUE)

