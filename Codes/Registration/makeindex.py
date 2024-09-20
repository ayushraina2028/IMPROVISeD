# automate_pymol.py
from pymol import cmd

with open("protein_name.txt", "r") as f:
    protein_name = f.read().strip()

print(f"Received protein name: {protein_name}")

# Load PDB files
cmd.load(f"../../Examples/{protein_name}/chain_A.pdb")
cmd.load(f"../../Examples/{protein_name}/chain_B.pdb")

# Run external script
cmd.run("../Registration/makeXnIndex.py")

# Execute makeXnIndex function
cmd.do(f"makeXnIndex(['chain_A', 'chain_B'], '../LocalizationBamdev/{protein_name}Result/Y_30.csv', ['../LocalizationBamdev/{protein_name}Result/tmp_indx1.txt', '../LocalizationBamdev/{protein_name}Result/tmp_indx2.txt'], 'run1_Y30')")

