#!/bin/bash

# Initialize sum_logs variable
sum_logs=0

# Read in abundance1.tsv file and loop over rows
while IFS=$'\t' read -r col1 col2; do
    # Check if value in column 2 is non-zero
    if (( $(echo "$col2 != 0" | bc -l) )); then
        # Add the absolute value of the logarithm of the non-zero value to the sum_logs variable
        sum_logs=$(echo "$sum_logs + a(l($col2))" | bc -l)
    fi
done < "data/references/simulation_short_read/abundance1.tsv"

echo "Sum of absolute logarithms: $sum_logs"
