#!/bin/bash

# Loop over all output subfolders
for subfolder in $(ls | grep "output$"); do
  # Extract prefix of subfolder (two digits before "_output")
  prefix=$(basename "$subfolder" "_output")
  echo "Processing $subfolder"

  # Remove any prior existing output files
  for prior in $(ls $subfolder/*all_R*); do
    rm $prior
  done

  for file in "${subfolder}/reads/"*".fq.gz"; do
    # Extract genome name and R1/R2 information from filename
    genome=$(basename "$file" ".fq.gz")
    read_num=${genome: -1}

    # Concatenate R1 and R2 files for current genome and append to output files
    if [ "$read_num" = "1" ]; then
      zcat "${subfolder}/reads/${genome}.fq.gz" >> "${subfolder}/${prefix}_all_R1.fq"
    elif [ "$read_num" = "2" ]; then
      zcat "${subfolder}/reads/${genome}.fq.gz" >> "${subfolder}/${prefix}_all_R2.fq"
    fi
    echo "Processed $genome"
  done

  # Compress the output file with pigz and remove .fq
  for bunch in $(ls $subfolder | grep "all_R"); do
    
    echo "Compressing $bunch"
    pigz -p 6 -c "${subfolder}/${bunch}" > "${subfolder}/${bunch}.gz"
    rm "${subfolder}/${bunch}"
    echo "Compressed $bunch"
  done

  echo "Bunch up complete"
done
