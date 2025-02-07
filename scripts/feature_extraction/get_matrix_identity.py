import argparse
import numpy as np
import pandas as pd
from Bio import pairwise2

def calculate_identity(seq1, seq2, mode='aligned', include_gaps=False):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    top_alignment = max(alignments, key=lambda x: x.score)
    aligned_seq1, aligned_seq2 = top_alignment[0], top_alignment[1]
    
    if include_gaps:
        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
    else:
        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-' and b != '-')
    
    if mode == 'shortest':
        normalization_length = min(len(seq1.replace('-', '')), len(seq2.replace('-', '')))
    elif mode == 'longest':
        normalization_length = max(len(seq1.replace('-', '')), len(seq2.replace('-', '')))
    elif mode == 'aligned':
        normalization_length = len(aligned_seq1)
    else:
        raise ValueError("Invalid mode. Choose 'shortest', 'longest', or 'aligned'.")
    
    return matches / normalization_length if normalization_length > 0 else 0

def optimized_pairwise_identity_matrix(sequences):
    n = len(sequences)
    identity_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i, n):
            identity = calculate_identity(sequences[i], sequences[j])
            identity_matrix[i][j] = identity_matrix[j][i] = identity
    
    return identity_matrix

def save_matrix(matrix, output_file):
    np.save(output_file, matrix)

def main():
    parser = argparse.ArgumentParser(description='Compute pairwise sequence identity matrix from sequences in a CSV.')
    parser.add_argument('--csv', type=str, required=True, help='Path to the CSV file containing sequences.')
    parser.add_argument('--subgroup', type=str, required=True, help='Subgroup to filter sequences by.')
    parser.add_argument('--output', type=str, required=True, help='File path to save the output identity matrix.')

    args = parser.parse_args()

    # Read CSV and filter sequences
    df = pd.read_csv(args.csv)
    filtered_sequences = df[df['order'] == args.subgroup]['sequence'].values

    if len(filtered_sequences) > 0:
        identity_matrix = optimized_pairwise_identity_matrix(filtered_sequences)
        save_matrix(identity_matrix, args.output)
        print(f"Identity matrix saved to {args.output}")
    else:
        print("No sequences found for the specified subgroup.")

if __name__ == "__main__":
    main()

