#!/bin/bash

PROTEIN_NAME="1dfj"

# Define the MATLAB script
MATLAB_SCRIPT="../bashMatlab"

# Run the MATLAB script using MATLAB command line
echo "Running MATLAB script..."
matlab -nodisplay -nosplash -r "addpath('$MATLAB_SCRIPT'); addpath('../LocalizationBamdev'); runLocalization('$PROTEIN_NAME'); exit;"

mv "../bashMatlab/tmp_indx1.txt" "../LocalizationBamdev/${PROTEIN_NAME}Result"
mv "../bashMatlab/tmp_indx2.txt" "../LocalizationBamdev/${PROTEIN_NAME}Result"
mv "../bashMatlab/Y_30.csv" "../LocalizationBamdev/${PROTEIN_NAME}Result"

echo "MATLAB tasks completed."
