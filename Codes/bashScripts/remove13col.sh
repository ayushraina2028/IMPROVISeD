#!/bin/bash

# Define the Python script
PYTHON_SCRIPT="../GenerateCLs/separateCols.py"

# Run the Python script
echo "Running Python script to remove columns..."
python3 "$PYTHON_SCRIPT"

echo "Column removal completed."

