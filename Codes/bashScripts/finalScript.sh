#!/bin/bash

# Define file paths
PROTEIN_NAME="1dfj"
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

# Separating 2 files
input_file="../crosslinks/${PROTEIN_NAME}CLs/eq_dists_prots_model1.csv"

# Output Files
output_file_A="../crosslinks/${PROTEIN_NAME}CLs/eq_dists_chain_A.csv"
output_file_B="../crosslinks/${PROTEIN_NAME}CLs/eq_dists_chain_B.csv"

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

# Define the Python script
PYTHON_SCRIPT="../GenerateCLs/separateCols.py"

# Run the Python script
echo "Running Python script to remove columns..."
python3 "$PYTHON_SCRIPT" "$PROTEIN_NAME"

echo "Column removal completed."

# Localization
# Define the MATLAB script
MATLAB_SCRIPT="../bashMatlab"

# Run the MATLAB script using MATLAB command line
echo "Running MATLAB script..."
matlab -nodisplay -nosplash -r "addpath('$MATLAB_SCRIPT'); runLocalization('$PROTEIN_NAME'); exit;"

mv "../bashMatlab/tmp_indx1.txt" "../LocalizationBamdev/${PROTEIN_NAME}Result"
mv "../bashMatlab/tmp_indx2.txt" "../LocalizationBamdev/${PROTEIN_NAME}Result"
mv "../bashMatlab/Y_30.csv" "../LocalizationBamdev/${PROTEIN_NAME}Result"

echo "MATLAB tasks completed."

# Making X_n_index
# Run the Python script with PyMOL
path="../Registration/makeindex.py"

# write protein name into a file
FILE_NAME="protein_name.txt"
echo "$PROTEIN_NAME" > "$FILE_NAME"

pymol -c -d "run $path"

rm "$FILE_NAME"

mv "run1_Y30_chain_A_indx.txt" "../Registration/${PROTEIN_NAME}Data/"
mv "run1_Y30_chain_A_xyz.txt" "../Registration/${PROTEIN_NAME}Data/"
mv "run1_Y30_chain_B_indx.txt" "../Registration/${PROTEIN_NAME}Data/"
mv "run1_Y30_chain_B_xyz.txt" "../Registration/${PROTEIN_NAME}Data/"
mv "run1_Y30crosslink_indx.txt" "../Registration/${PROTEIN_NAME}Data/"
mv "run1_Y30crosslink_xyz.txt" "../Registration/${PROTEIN_NAME}Data/"

# Registration
# Define the path to the MATLAB script
MATLAB_SCRIPT="../bashMatlab/register.m"

FILE_NAME="protein_name.txt"
echo "$PROTEIN_NAME" > "$FILE_NAME"

# Run the MATLAB script using MATLAB in non-interactive mode
matlab -batch "run('$MATLAB_SCRIPT')"

rm "$FILE_NAME"

mv "../bashMatlab/x_n_index_run1_Y30.mat" "../Registration/${PROTEIN_NAME}Data/"
