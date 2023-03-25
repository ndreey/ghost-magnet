#!/bin/bash

# Loop over all fq.gz files in the current directory
for file in *.fq.gz
do
  # Extract the file name without the extension
  name=$(basename "$file" .fq.gz)
  
  # Concatenate the contents of the current file to the output file
  zcat "$file" >> 1_hc_raw_reads.fq
  
  echo "Processed $name"
done

# Compress the output file with gzip
echo "Compressing bunch up"
pigz -p 6 -c 1_hc_raw_reads.fq > 1_hc_raw_reads.fq.gz

echo "Done!"
