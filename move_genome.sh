#!/bin/bash

# Set the new path for the genome files
new_path="/mnt/c/Users/andbo/thesis_andbo/CAMISIM/platanthera_mock/reference_genomes/"

# Read the second column of the file "01_genome_to_id.tsv" and extract the old paths
old_paths=($(cut -f 2 "/mnt/c/Users/andbo/thesis_andbo/CAMISIM/platanthera_mock/01_simulation/01genome_to_id.tsv"))

# Loop through the old paths and copy the files to the new path
for old_path in "${old_paths[@]}"; do
    # Extract the filename from the old path
    filename=$(basename "$old_path")
    # Copy the file to the new path
    cp "$old_path" "${new_path}${filename}"
done

