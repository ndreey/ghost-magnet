#!/bin/bash

# Set the paths for input and output.
directory_path="data/references/simulation_short_read"
output_directory="data/processed"

echo "Calculating abundance sums for files in $directory_path..."

# Forces the deletion of prior reports.
report_file="$output_directory/report_abundance.tsv"
rm -f "$report_file"

# Add header to report file
echo -e "file\tsum\tsum_abs_log" > "$report_file"

# Loops through all the abundance files
for file in "$directory_path"/abundance*.tsv; do
    # Good practice to check if file exists to avoid simple errors.
    if [[ -f $file ]]; then
        file_name=$(basename "$file")
	# For each row, sum the column 2 val and then print the sum.
	# because the awk command is inside $() the output is stored in the $read_sum variable.
        read_sum=$(awk '{sum += $2} END {print sum}' "$file")

	# For each col 2 value that is not 0, we sum the abs(log(val)).
        # Each line is split into two by the Internal Field Separator (IFS).
        sum_logs=0
        while IFS=$'\t' read -r col1 col2; do
            # Check if value in column 2 is non-zero (we use math library as 0.0 = 0)
            if (( $(echo "$col2 != 0" | bc -l) )); then
                # Add the absolute value of the logarithm of the non-zero value to the sum_logs variable
                sum_logs=$(echo "$sum_logs + a(l($col2))" | bc -l)
            fi
        # Ends the while loop.
        done < "$file"

	# For each file we write a tsv line
        echo -e "$file_name\t$read_sum\t$sum_logs" >> "$report_file"
    fi
done

echo "Report generated at $report_file"
