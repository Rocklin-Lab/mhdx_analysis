import pandas as pd
import numpy as np
import os
import re
import argparse
from itertools import combinations


def find_structure_segments(ss, N=2):
    segments = [(match.start(), match.end()) for match in re.finditer(r'[HE]+', ss)]
    new_segments = []
    for start, end in segments:
        if end - start < N:
            continue
        e_to_h_transition = ss.find('HE', start, end)
        h_to_e_transition = ss.find('EH', start, end)
        if e_to_h_transition != -1:
            new_segments.append((start, e_to_h_transition+1))
            new_segments.append((e_to_h_transition+1, end))
        elif h_to_e_transition != -1:
            new_segments.append((start, h_to_e_transition+1))
            new_segments.append((h_to_e_transition+1, end))
        else:
            new_segments.append((start, end))
    return sorted(new_segments)

def get_largest_deltas(tuples, N):

    deltas = [(t[1] - t[0], t) for t in tuples]

    # Calculate deltas and sort tuples by their deltas
    sorted_tuples = sorted(deltas, key=lambda x: x[0], reverse=True)
    
    # If N is larger than or equal to the number of tuples, return all tuples
    if N >= len(sorted_tuples):
        return [sorted([t[1] for t in sorted_tuples])]
    
    # Determine if there's a tie at the Nth position
    Nth_delta = sorted_tuples[N-1][0]
    tied_tuples_starting_from_N = [t for t in sorted_tuples if t[0] == Nth_delta]
    
    # Count how many tuples have deltas larger than the Nth_delta
    larger_than_N_count = sum(1 for t in sorted_tuples if t[0] > Nth_delta)
    
    # If there's no tie or if N is 1, simply return the N largest deltas
    if len(tied_tuples_starting_from_N) == 1 or N == 1:
        return [sorted([t[1] for t in sorted_tuples[:N]])]
    
    # Handle ties
    # Number of ties to select from tied tuples to make up to N
    ties_to_select = N - larger_than_N_count
    tie_combinations = list(combinations(tied_tuples_starting_from_N, ties_to_select))
    
    # Generate all possible combinations, preserving the order within each list
    results = []
    for combo in tie_combinations:
        selected_tuples = sorted_tuples[:larger_than_N_count] + list(combo)
        results.append(sorted([t[1] for t in selected_tuples]))
    
    results.sort(key=lambda x: sum(t[1] - t[0] for t in x), reverse=True)
    
    return results


def get_largest_deltas(tuples, N, min_delta=2):
    deltas = [(t[1] - t[0], t) for t in tuples]

    # Filter tuples based on max_delta
    filtered_tuples = [d for d in deltas if d[0] >= min_delta]

    # Calculate deltas and sort tuples by their deltas
    sorted_tuples = sorted(filtered_tuples, key=lambda x: x[0], reverse=True)
    
    # Generate all possible combinations
    results = []
    for combo in combinations(sorted_tuples, N):
        result = sorted([t[1] for t in combo])
        results.append(result)
        
    return results


def find_topn_segments(dssp, topology):   


    n = len(topology)
    
    segments = find_structure_segments(dssp)
    
    all_possible_segments_lists = get_largest_deltas(segments, n)

    for segments_list in all_possible_segments_lists:
        if len(segments_list) != len(topology):
            continue
        match = True
        for index, (i, j) in enumerate(segments_list):
            if dssp[i] != topology[index]:
                match = False
                continue
        if match:
            return segments_list
    
    return []


def check_dssp_pf_consistency(dssp, topology):

    if topology is np.nan or topology == 'nan':
        print("Nan for topology found. Returning False")
        return False
    elif dssp is np.nan or dssp == 'nan':
        print("Nan for dssp found. Returning False")
        return False
    
    segments_list = find_topn_segments(dssp, topology)
    
    if len(segments_list) != len(topology):
        return False
    
    match = True
    for index, (i, j) in enumerate(segments_list):
        if dssp[i] != topology[index]:
            match = False
            continue
    if match:
        return True
    





def get_prop(prop, segment):
    i, j = segment
    return prop[i:j]

def assign_values(row, columns, values):
    # Assign values to the specified columns for the given row
    row[columns] = pd.Series(values)
    return row


import pandas as pd
import numpy as np

# Assuming other required functions (find_structure_segments, find_topn_segments, check_dssp_pf_consistency, etc.) are correctly implemented above

def get_properties(row, prop, prop_name=None, calculate_mean=True, calculate_max=True, calculate_min=True, calculate_iqr=True, calculate_top3_mean=True, calculate_bottom3_mean=True):
    topology = row["topology"]
    dssp = row["dssp"]
    prop_values = row[prop] if prop in row and isinstance(row[prop], list) else []

    if prop_name is None:
        prop_name = prop

    # Initialize containers for column names and their corresponding values
    cols = []
    values = []

    # Check for secondary structure consistency
    if check_dssp_pf_consistency(dssp, topology):
        segments = find_topn_segments(dssp, topology)

        # Segment-specific features
        for i, char in enumerate(topology):
            # Combine values for all segments
            segment_prop_values = get_prop(row[prop], segments[i])
            

            # Calculating statistics for each segment
            segment_stats = {
                f"{char}{i+1}_{prop_name}_mean": np.mean(segment_prop_values) if segment_prop_values else np.nan,
                f"{char}{i+1}_{prop_name}_max": np.max(segment_prop_values) if segment_prop_values else np.nan,
                f"{char}{i+1}_{prop_name}_min": np.min(segment_prop_values) if segment_prop_values else np.nan,
                f"{char}{i+1}_{prop_name}_iqr": np.percentile(segment_prop_values, 75) - np.percentile(segment_prop_values, 25) if segment_prop_values else np.nan,
                f"{char}{i+1}_{prop_name}_top3_mean": np.mean(sorted(segment_prop_values)[-3:]) if len(segment_prop_values) >= 3 else np.nan,
                f"{char}{i+1}_{prop_name}_bottom3_mean": np.mean(sorted(segment_prop_values)[:3]) if len(segment_prop_values) >= 3 else np.nan,
            }
            cols.extend(segment_stats.keys())
            values.extend(segment_stats.values())


        global_stats = {
            f"{prop_name}_global_mean": np.mean(prop_values) if len(prop_values) > 0 else np.nan,
            f"{prop_name}_global_max": np.max(prop_values) if len(prop_values) > 0 else np.nan,
            f"{prop_name}_global_min": np.min(prop_values) if len(prop_values) > 0 else np.nan,
            f"{prop_name}_global_iqr": np.percentile(prop_values, 75) - np.percentile(prop_values, 25) if len(prop_values) > 0 else np.nan,
            f"{prop_name}_global_top3_mean": np.mean(sorted(prop_values)[-3:]) if len(prop_values) >= 3 else np.nan,
            f"{prop_name}_global_bottom3_mean": np.mean(sorted(prop_values)[:3]) if len(prop_values) >= 3 else np.nan,
        }
        cols.extend(global_stats.keys())
        values.extend(global_stats.values())

        # Combining values across all segments for global secondary structure (ss) stats
        all_segment_values = np.concatenate([prop_values[s[0]:s[1]] for s in segments]) if segments else []

        ss_stats = {
            f"{prop_name}_ss_mean": np.mean(all_segment_values) if all_segment_values.size > 0 else np.nan,
            f"{prop_name}_ss_max": np.max(all_segment_values) if all_segment_values.size > 0 else np.nan,
            f"{prop_name}_ss_min": np.min(all_segment_values) if all_segment_values.size > 0 else np.nan,
            f"{prop_name}_ss_iqr": np.percentile(all_segment_values, 75) - np.percentile(all_segment_values, 25) if all_segment_values.size > 0 else np.nan,
            f"{prop_name}_ss_top3_mean": np.mean(sorted(all_segment_values)[-3:]) if all_segment_values.size >= 3 else np.nan,
            f"{prop_name}_ss_bottom3_mean": np.mean(sorted(all_segment_values)[:3]) if all_segment_values.size >= 3 else np.nan,
        }
        cols.extend(ss_stats.keys())
        values.extend(ss_stats.values())

    return pd.Series(values, index=cols)

def process_local_features(df, topology, output=None):
    """
    Process columns that are lists in df, calculate statistics, and save to a new DataFrame.
    """
    df['topology'] = topology  # Assigning topology to each row

    # Identify columns with lists
    columns_with_lists = [col for col in df.columns if df[col].apply(lambda x: isinstance(x, list)).all()]

    # Prepare an empty DataFrame to collect new properties
    new_df = pd.DataFrame()

    # Process each column with lists
    for prop in columns_with_lists:
        # Apply the get_properties function to each row for the current property
        stats_df = df.apply(get_properties, prop=prop, axis=1)
        
        # Concatenate the resulting statistics to the new_df DataFrame
        if new_df.empty:
            new_df = stats_df
        else:
            new_df = pd.concat([new_df, stats_df], axis=1)

    # Ensure no duplicate columns
    new_df = new_df.loc[:, ~new_df.columns.duplicated()]

    # Add sequence column to new_df if it's not already included
    if 'sequence' in df.columns and 'sequence' not in new_df.columns:
        new_df['sequence'] = df['sequence']

    return new_df       


def add_relative_metrics(df, topology):

     
    # Compute relativeToMax, SSmax and SSmin
    # Extract unique suffixes from column names
    suffixes = set(['_'.join(i.split("_")[1:]) for i in df.select_dtypes(include='number').columns if f"{topology[0]}1_" in i])
    
    # Initialize empty columns for min and max values
    for suffix in suffixes:
        min_col_name = f'{suffix}_SSmin'
        max_col_name = f'{suffix}_SSmax'
        df[min_col_name] = df.filter(like=suffix).min(axis=1)
        df[max_col_name] = df.filter(like=suffix).max(axis=1)
        
        # Add new columns for features relative to max
        for col in df.filter(like=suffix).columns:
            if col.endswith(suffix):
                relative_col_name = f'{col}_relativeToMax'
                df[relative_col_name] = df[col] / df[max_col_name]

    return df

def main(numeric_file_path, list_file_path, topology, output_file_path):
    # Load numeric and list DataFrames
    df_numeric = pd.read_json(numeric_file_path)
    df_lists = pd.read_json(list_file_path)

    processed_list_df = process_local_features(df_lists, topology)

    # Concatenate processed list DataFrame with numeric DataFrame
    final_df = pd.concat([df_numeric, processed_list_df], axis=1)
    # Add relative metrics
    final_df = add_relative_metrics(final_df, topology)
    # Drop eventual duplicates
    final_df = final_df.loc[:,~final_df.columns.duplicated()]
   

    # Save the final DataFrame to the specified output file
    final_df.to_json(output_file_path, orient='records')
    print(f"Final DataFrame saved to {output_file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and merge numeric and list features from JSON files.")
    parser.add_argument('--numeric_path', type=str, required=True, help="Path to the JSON file with numeric features.")
    parser.add_argument('--list_path', type=str, required=True, help="Path to the JSON file with list features.")
    parser.add_argument('--topology', type=str, required=True, help="Topology")
    parser.add_argument('--output', type=str, required=True, help="Output path for the merged DataFrame.")
    
    args = parser.parse_args()
    
    # Assuming that each row in the list DataFrame has a 'topology' column
    # If this is not the case, you'll need to adjust the logic to supply or extract topology information
    main(args.numeric_path, args.list_path, args.topology, args.output)
