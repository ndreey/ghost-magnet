#!/bin/bash

directory_path="data/references/simulation_short_read"
output_directory="data/processed"

echo "Calculating abundance sums for files in $directory_path..."

report_file="$output_directory/report_abundance2.tsv"
rm -f "$report_file"

# Add header to report file
echo -e "file\tsum\tsum_abs_log" > "$report_file"

for file in "$directory_path"/abundance*.tsv; do
    if [[ -f $file ]]; then
        file_name=$(basename "$file")
        read_sum=$(awk '{sum += $2} END {print sum}' "$file")
        if (( $(echo "$col2 != 0" | bc -l) )); then
            # Add the absolute value fo the log value
            sum_logs=$(echo "$sum_logs + a(l($col2))" | bc -l)
        fi
        echo -e "$file_name\t$read_sum\t$sum_logs" >> "$report_file"
    fi
done

echo "Report generated at $report_file"