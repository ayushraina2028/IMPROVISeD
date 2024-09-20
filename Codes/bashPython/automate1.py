import os
import pymol
from pymol import cmd
import sys

# Get command-line arguments
pdb_folder = sys.argv[1]
pdbName1 = sys.argv[2]
pdbName2 = sys.argv[3]
output_csv = sys.argv[4]

# Define file paths
pdb_file1 = os.path.join(pdb_folder, pdbName1)
pdb_file2 = os.path.join(pdb_folder, pdbName2)

# Step 1: Start PyMOL without GUI
print("Starting PyMOL...")
pymol.finish_launching(['pymol', '-cq'])
print("PyMOL started successfully.")

# Step 2: Load the PDB files
print(f"Loading PDB files: {pdb_file1} and {pdb_file2}")
cmd.load(pdb_file1, 'chain_A')
cmd.load(pdb_file2, 'chain_B')
print("PDB files loaded successfully.")

# Step 3: Run the Python script inside PyMOL
print("Running the createCrosslinksNEq_v2.py script...")
cmd.run('../GenerateCLs/createCrosslinksNEq_v2.py')
print("Python script executed.")

# Step 4: Run the PyMOL command to create crosslinks and save the CSV
print(f"Creating crosslinks and saving to {output_csv}...")
cmd.do("createEqnCross(['chain_A', 'chain_B'], '../crosslinks/1dfjCLs/{}', 'model1')".format(output_csv))
print(f"Crosslinks created and saved to {output_csv}.")

# Step 5: Close PyMOL
print("Closing PyMOL...")
cmd.quit()
print("PyMOL closed.")

