#!/bin/bash

# Define file paths
PROTEIN_NAME="lcn2"

FOLDER_LOCATION="../../Examples"
FOLDER_PATH="$FOLDER_LOCATION/$PROTEIN_NAME"

# Find the two .pdb files in the folder
pdb_files=($(find "$FOLDER_PATH" -maxdepth 1 -name "*.pdb"))

# Check if exactly two PDB files are found
if [ ${#pdb_files[@]} -ne 2 ]; then
    echo "Error: The folder must contain exactly two .pdb files."
    exit 1
fi

# Extract the base names of the PDB files (without extensions)
chain_A=$(basename "${pdb_files[0]}" .pdb)
chain_B=$(basename "${pdb_files[1]}" .pdb)

echo "ayush: ", $chain_A $chain_B

cp "../../Examples/${PROTEIN_NAME}/crosslinks.csv" "../crosslinks/${PROTEIN_NAME}CLs/${chain_A}_crosslink30_${chain_B}_LYS_Ca.csv"


# Output filename and python script to run
OUTPUT_CSV="${chain_A}_crosslink30_${chain_B}_LYS_Ca.csv"
PYTHON_SCRIPT="../bashPython/automate1.py"

# Run pymol in command line
echo "Running pymol with the ${PYTHON_SCRIPT}..."
pymol -cq -d "import sys; sys.argv = ['$PYTHON_SCRIPT','$FOLDER_PATH', '${chain_A}.pdb', '${chain_B}.pdb', '$OUTPUT_CSV', '$PROTEIN_NAME']; exec(open('$PYTHON_SCRIPT').read())"

CSV_SOURCE="eq_dists_prots_model1.csv"
DESTINATION_CSV="../crosslinks/${PROTEIN_NAME}CLs/"

echo "Checking if CSV file exists and moving it..."
if [ -f "$CSV_SOURCE" ]; then
    mv "$CSV_SOURCE" "$DESTINATION_CSV"
    echo "CSV file moved to ${DESTINATION_CSV}"
else
    echo "CSV file $CSV_SOURCE not found."
fi
