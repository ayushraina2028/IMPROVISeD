import pandas as pd
import sys
def remove_columns(file_path):
    # Load the CSV file
    df = pd.read_csv(file_path)

    # Drop the 'prot1' and 'prot2' columns
    df = df.drop(['prot1', 'prot2'], axis=1)

    # Save the modified DataFrame back to the same file
    df.to_csv(file_path, index=False)
    print(f"Processed file: {file_path}")

if __name__ == "__main__":
    protein_name = sys.argv[1]
    files = [
        f'../crosslinks/{protein_name}CLs/eq_dists_chain_A.csv',
        f'../crosslinks/{protein_name}CLs/eq_dists_chain_B.csv'
    ]

    for file in files:
        remove_columns(file)

