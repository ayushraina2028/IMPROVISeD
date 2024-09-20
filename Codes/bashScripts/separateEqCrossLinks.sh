#!/bin/bash

# Input file with the data
PROTEIN_NAME="1dfj"
input_file="../crosslinks/1dfjCLs/eq_dists_prots_model1.csv"

# Output files
output_file_A="../crosslinks/1dfjCLs/eq_dists_chain_A.csv"
output_file_B="../crosslinks/1dfjCLs/eq_dists_chain_B.csv"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file $input_file not found."
    exit 1
fi

# Extract the header line
header=$(head -n 1 "$input_file")

# Extract lines that contain chain_A data and prepend header
echo "$header" > "$output_file_A"
grep ",chain_A," "$input_file" >> "$output_file_A"

# Extract lines that contain chain_B data and prepend header
echo "$header" > "$output_file_B"
grep ",chain_B," "$input_file" >> "$output_file_B"

# Check if files were created successfully
if [ -f "$output_file_A" ]; then
    echo "File for chain_A data created: $output_file_A"
else
    echo "Error: chain_A data file not created."
fi

if [ -f "$output_file_B" ]; then
    echo "File for chain_B data created: $output_file_B"
else
    echo "Error: chain_B data file not created."
fi
