import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np

def read_fasta_file(fasta_file):
    """Reads the first sequence from a FASTA file."""
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            return str(record.seq)
    return None

def read_json_for_plddt(json_file):
    """Reads the 'plddt' column from a JSON file."""
    df = pd.read_json(json_file)
    return df['plddt'].tolist()

def main(fasta, plddt, output):
    sequence = read_fasta_file(fasta)
    df = pd.read_json(plddt)
    plddt_values = df['plddt'].tolist()
    pae_values = df['pae'].apply(lambda x: np.mean(x)).to_list()
    
    # Combine into a DataFrame
    df = pd.DataFrame({
        'sequence': [sequence],
        'plddt': [plddt_values],
        'pae': [pae_values]
    })
    
    # Convert DataFrame to JSON
    df.to_json(output, orient='records')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a FASTA file and a JSON file to create a dataframe.")
    parser.add_argument("--fasta", type=str, help="The path to the FASTA file containing the sequence.")
    parser.add_argument("--plddt", type=str, help="The path to the JSON file containing the 'plddt' data.")
    parser.add_argument("--output", type=str, help="The path to the OUTPUT file.")
    
    args = parser.parse_args()
    
    main(args.fasta, args.plddt, args.output)

