#!/bin/bash

# Loop through all R1 files in reads/00_raw/
for r1_file in reads/00_raw/*_all_R1.fq.gz; do
    # Extract the sample ID from the R1 file name
    sample=$(basename "$r1_file" | cut -d'_' -f1)
    echo "Processing sample $sample..."
    echo $r1_file

    # Define the corresponding R2 file name
    r2_file="reads/00_raw/${sample}_all_R2.fq.gz"

    # Run fastp with the desired parameters
    fastp --thread 6 -c -f 4 -t 2 -F 4 -T 2 -r -W 4 -M 27 -i "$r1_file" -I "$r2_file" \
    -o "reads/01_trimmed/${sample}_trim_R1.fq.gz" -O "reads/01_trimmed/${sample}_trim_R2.fq.gz" \
    -l 36 -h "analysis/fastp_reports/${sample}_fastp.html" -j "analysis/fastp_reports/${sample}_fastp.json"

    echo "Finished processing sample $sample."
done
