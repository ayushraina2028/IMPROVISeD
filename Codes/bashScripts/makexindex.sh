#!/bin/bash

PROTEIN_NAME="1dfj"

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
