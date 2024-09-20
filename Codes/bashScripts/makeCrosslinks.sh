#!/bin/bash

# Hardcoded path to the folder containing the two PDB files

proteinName="1dfj"
folderLocation="../../Examples"

folder_path="$folderLocation/$proteinName"

# Find the two .pdb files in the folder
pdb_files=($(find "$folder_path" -maxdepth 1 -name "*.pdb"))

# Check if exactly two PDB files are found
if [ ${#pdb_files[@]} -ne 2 ]; then
    echo "Error: The folder must contain exactly two .pdb files."
    exit 1
fi

# Extract the base names of the PDB files (without extensions)
chain_A=$(basename "${pdb_files[0]}" .pdb)
chain_B=$(basename "${pdb_files[1]}" .pdb)

# Run the Python script to generate crosslinks
python ../GenerateCLs/makeCrosslinks.py -p "$folder_path/$chain_A.pdb,$folder_path/$chain_B.pdb" -b "$chain_A,$chain_B" -d 30

# Move the generated CSV file to the specified directory
mv "${chain_A}_crosslink30_${chain_B}_LYS_Ca.csv" "../crosslinks/${proteinName}CLs/"

# Print success message
echo "Crosslink generation and file move completed successfully."

