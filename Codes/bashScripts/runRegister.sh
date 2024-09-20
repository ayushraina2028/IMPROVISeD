#!/bin/bash

PROTEIN_NAME="1dfj"

# Define the path to the MATLAB script
MATLAB_SCRIPT="../bashMatlab/register.m"

FILE_NAME="protein_name.txt"
echo "$PROTEIN_NAME" > "$FILE_NAME"

# Run the MATLAB script using MATLAB in non-interactive mode
matlab -batch "run('$MATLAB_SCRIPT')"

rm "$FILE_NAME"

mv "../bashMatlab/x_n_index_run1_Y30.mat" "../Registration/${PROTEIN_NAME}Data/"

