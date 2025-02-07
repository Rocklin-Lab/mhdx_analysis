import argparse
import iFeatureOmegaCLI
import pandas as pd

def main():
    # Initialize argparse
    parser = argparse.ArgumentParser(description='Process input and output files for sequence analysis.')
    parser.add_argument('-i', '--fasta', type=str, required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file path')

    args = parser.parse_args()

    # Use the input file as the argument for iProtein
    protein = iFeatureOmegaCLI.iProtein(args.fasta)

    # Process descriptors and concatenate their encoding data
    descriptors = ["AAC", "GAAC", "CTDC"]
    ls = []
    for descriptor in descriptors:
        protein.get_descriptor(descriptor)
        # Apply a suffix based on the descriptor type, if necessary
        suffix_ = "_frac" if "type 1" in descriptor else "_count" if "type 2" in descriptor else ""
        protein.encodings.columns = [col + suffix_ for col in protein.encodings.columns]
        ls.append(protein.encodings)

    # Concatenate features data and reset index
    df_features = pd.concat(ls, axis=1)
    df_features = df_features.reset_index().rename(columns={'index': 'name'})

    # Output the features to a JSON file
    df_features.to_json(args.output)

if __name__ == "__main__":
    main()
