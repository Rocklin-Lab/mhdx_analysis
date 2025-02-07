import argparse
import pandas as pd
import numpy as np

def main(files, output_numeric, output_lists):
    # Initialize an empty list to store DataFrames
    dfs_numeric, dfs_lists = [], []

    # Read each JSON file and process if necessary
    for file_path in files:
        if 'json' in file_path:
            df = pd.read_json(file_path)  # Not using lines=True
        elif 'csv' in file_path:
            df = pd.read_csv(file_path)
#            df = df.drop("dssp", axis=1)

    #    print(file_path, len(df.columns))
        # Check and process 'name' column
        if 'name' in df.columns:
            df['name'] = df['name'].str.replace(r'\.pdb$', '', regex=True)
        if 'sequence' in df.columns:
            sequence = df['sequence'].values[0]
        if ('dssp' in df.columns) and ('abego' in file_path):
            dssp  = df['dssp'].values[0]
        if 'abego' in df.columns:
            abego = df['abego'].values[0]

        # Filter columns that start with 'array_'
        columns_to_drop = [col for col in df.columns if col.startswith('array_')]
        # Drop these columns from the DataFrame, if they exist
        df.drop(columns=columns_to_drop, inplace=True)

        # Rename columns: remove "array_" prefix and replace "predicted_z_scores" with "adopt"
        df.columns = df.columns.str.replace('predicted_z_scores', 'adopt', regex=False)
        if 'sasa' in file_path:
            df.columns = df.columns.str.replace('total', 'SASAtotal', regex=False)
            df.columns = df.columns.str.replace('relativeTotal', 'SASArelativeTotal', regex=False)
        numeric_cols = df.select_dtypes(include='number').columns
        lists_cols = [col for col in df.columns if any(isinstance(x, list) for x in df[col])]
        if len(numeric_cols) > 0: 
            dfs_numeric.append(df[numeric_cols])
        if len(lists_cols) > 0:
            dfs_lists.append(df[lists_cols])

#    print(len(dfs))

    # Concatenate all DataFrames, assuming similar structure/columns
#    concat_df = pd.concat(dfs, axis=1)
    
#    sequence = concat_df['sequence'].values[0][0]
#    dssp = concat_df['dssp'].values[0][0]

    ## Remove duplicated columns, keeping the first occurrence
    ##concat_df = concat_df.loc[:, ~concat_df.columns.duplicated()]
    #print(len(concat_df.columns))
    
    # Split the DataFrame into two based on data types
#    df_numeric = concat_df.select_dtypes(include=[np.number])
    df_numeric = pd.concat(dfs_numeric, axis=1)

#    print('numeric:', len(df_numeric.columns))


    # Identify and isolate columns with lists manually
 #   columns_with_lists = [col for col in concat_df.columns if any(isinstance(x, list) for x in concat_df[col])]
 #   df_lists = concat_df[columns_with_lists].copy()
    df_lists = pd.concat(dfs_lists, axis=1)

#    print('lists:', len(df_lists.columns))

    # Ensure 'sequence' column is in df_numeric if it exists
#    if 'sequence' in concat_df.columns and 'sequence' not in df_numeric.columns:
    df_numeric['sequence'] = sequence
    df_numeric['dssp'] = dssp
    df_numeric['abego'] = abego
#    if 'sequence' in concat_df.columns and 'sequence' not in df_lists.columns:
    df_lists['dssp'] = dssp
    df_lists['sequence'] = sequence
    df_lists['abego'] = abego


    # Save the DataFrames without lines=True
    df_numeric.to_json(output_numeric, orient='records')
    df_lists.to_json(output_lists, orient='records')

    print("DataFrames saved successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process JSON files and output numeric and list data.")
    parser.add_argument('files', nargs='+', help='List of file paths')
    parser.add_argument('--output_numeric', type=str, required=True, help='Output file name for numeric data')
    parser.add_argument('--output_lists', type=str, required=True, help='Output file name for list data')
    
    args = parser.parse_args()
    
    main(args.files, args.output_numeric, args.output_lists)

