import numpy as np
import pandas as pd
import argparse
from Bio.SeqUtils import molecular_weight

def calculate_monoisotopic_mass(sequence):
    """Calculate the monoisotopic mass of a peptide sequence using Biopython."""
    # Calculate the mass assuming the sequence is a standard linear peptide (not circular).
    return molecular_weight(sequence, seq_type='protein', monoisotopic=True)

def compute_ppm_difference_matrix(masses):
    """Compute the ppm difference matrix for a list of monoisotopic masses."""
    n = len(masses)
    ppm_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            ppm_matrix[i][j] = abs(masses[i] - masses[j]) / masses[i] * 1e6 if masses[i] != 0 else 0
    
    return ppm_matrix

def save_matrix(matrix, output_file):
    """Save the matrix to a file."""
    np.save(output_file, matrix)

def main():
    parser = argparse.ArgumentParser(description='Compute ppm difference matrix from sequences in a CSV.')
    parser.add_argument('--csv', type=str, required=True, help='Path to the CSV file containing sequences.')
    parser.add_argument('--subgroup', type=str, required=True, help='Subgroup to filter sequences by.')
    parser.add_argument('--output', type=str, required=True, help='File path to save the output ppm matrix.')

    args = parser.parse_args()

    # Read CSV and filter sequences
    df = pd.read_csv(args.csv)
    filtered_sequences = df[df['order'] == args.subgroup]['sequence'].values

    if len(filtered_sequences) > 0:
        # Calculate monoisotopic masses
        masses = [calculate_monoisotopic_mass(seq) for seq in filtered_sequences]
        # Compute ppm matrix
        ppm_matrix = compute_ppm_difference_matrix(masses)
        # Save matrix
        save_matrix(ppm_matrix, args.output)
        print(f"PPM matrix saved to {args.output}")
    else:
        print("No sequences found for the specified subgroup.")

if __name__ == "__main__":
    main()

