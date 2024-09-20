#!/bin/bash

# Define the Python script
PROTEIN_NAME="lcn2"
PYTHON_SCRIPT="../GenerateCLs/separateCols.py"

# Run the Python script
echo "Running Python script to remove columns..."
python3 "$PYTHON_SCRIPT" "$PROTEIN_NAME"

echo "Column removal completed."

