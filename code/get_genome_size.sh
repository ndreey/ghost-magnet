#!/bin/bash

# Set the paths for input and output.
directory_path="data/references/genomes/source_genomes"
output_directory="data/processed"

echo "Calculating genome sizes of .fasta in $directory_path..."

# Forces the deletion of prior reports.
report_file="$output_directory/report_genome_sizes.tsv"
rm -f "$report_file"

# Add header to report file
echo -e "genome\tsize" > "$report_file"

# Loops through all the genomes
for file in "$directory_path"/*.fasta; do
    # Good practice to check if file exists to avoid simple errors.
    if [[ -f $file ]]; then
        file_name=$(basename "$file")
        # All characters in each row that does not begin with > is counted
        genome_sum=$(grep -v "^>" $file | tr -d '\n' | wc -c)      

        # For each file we write a tsv line
        echo -e "$file_name\t$genome_sum" >> "$report_file"
    fi
done

echo "$report_file was created"
