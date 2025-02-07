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
    
    return False


def get_prop(prop, segment):
    i, j = segment
    return prop[i:j]

def assign_values(row, columns, values):
    # Assign values to the specified columns for the given row
    row[columns] = pd.Series(values)
    return row


def get_properties(row, prop, prop_name=None, calculate_mean=True, calculate_max=True, calculate_min=True, calculate_iqr=True, calculate_top3_mean=True, calculate_bottom3_mean=True):
#    name = row["name"]
    topology = row["topology"]
    dssp = row["dssp"]

    if prop_name is None:
        prop_name = prop

    cols, values = [], []

    cols.append(f"{prop_name}_global_mean")
    cols.append(f"{prop_name}_global_max")
    cols.append(f"{prop_name}_global_min")
    cols.append(f"{prop_name}_global_iqr")
    cols.append(f"{prop_name}_global_top3_mean")
    cols.append(f"{prop_name}_global_bottom3_mean")
        
    values.append(np.mean(row[prop]))
    values.append(np.max(row[prop]))
    values.append(np.min(row[prop]))
    values.append(np.percentile(row[prop], 75) - np.percentile(row[prop], 25))
    values.append(np.mean(sorted(row[prop])[-3:]))
    values.append(np.mean(sorted(row[prop])[:3]))
    
    
    if check_dssp_topology_consistency(dssp, topology):
        segments = find_topn_segments(dssp, topology)


        # Combine values for all segments
        all_segment_values = [get_prop(row[prop], segment) for segment in segments]
        all_segment_values = np.concatenate(all_segment_values)

        cols.append(f"{prop_name}_ss_mean")
        cols.append(f"{prop_name}_ss_max")
        cols.append(f"{prop_name}_ss_min")
        cols.append(f"{prop_name}_ss_iqr")
        cols.append(f"{prop_name}_ss_top3_mean")
        cols.append(f"{prop_name}_ss_bottom3_mean")

        values.append(np.mean(all_segment_values))
        values.append(np.max(all_segment_values))
        values.append(np.min(all_segment_values))
        values.append(np.percentile(all_segment_values, 75) - np.percentile(all_segment_values, 25))
        values.append(np.mean(sorted(all_segment_values)[-3:]))
        values.append(np.mean(sorted(all_segment_values)[:3]))


        for i, char in enumerate(topology):
            
            col_mean = f"{char}{i+1}_{prop_name}_mean"
            col_max = f"{char}{i+1}_{prop_name}_max"
            col_min = f"{char}{i+1}_{prop_name}_min"
            col_iqr = f"{char}{i+1}_{prop_name}_iqr"
            col_top3_mean = f"{char}{i+1}_{prop_name}_top3_mean"
            col_bottom3_mean = f"{char}{i+1}_{prop_name}_bottom3_mean"

            if calculate_mean:
                cols.append(col_mean)
#                print(row[prop], segments[i], get_prop(row[prop], segments[i]))
                values.append(np.mean(get_prop(row[prop], segments[i])))

            if calculate_max:
                cols.append(col_max)
                values.append(np.max(get_prop(row[prop], segments[i])))

            if calculate_min:
                cols.append(col_min)
                values.append(np.min(get_prop(row[prop], segments[i])))

            if calculate_iqr:
                cols.append(col_iqr)
                prop_values = get_prop(row[prop], segments[i])
                values.append(np.percentile(prop_values, 75) - np.percentile(prop_values, 25))

            if calculate_top3_mean:
                cols.append(col_top3_mean)
                prop_values_sorted = sorted(get_prop(row[prop], segments[i]))
                values.append(np.mean(prop_values_sorted[-3:]))

            if calculate_bottom3_mean:
                cols.append(col_bottom3_mean)
                prop_values_sorted = sorted(get_prop(row[prop], segments[i]))
                values.append(np.mean(prop_values_sorted[:3]))


        print(len(cols), len(values))
        print(cols, values)


        # Modify the DataFrame in place
        row = assign_values(row.copy(), cols, values)

    return row


def load_and_merge_dataframes(df_features_path, df_properties_path, only_designs=False):
    
    df_features = pd.read_json(df_features_path) #.drop_duplicates("sequence").reset_index(drop=True) 
    
    df_prop = pd.read_json(df_properties_path)

    print(df_features.keys())
    print(df_prop.keys())

    if ("name" in df_prop.keys()) and ("sequence" in df_prop.keys()) and ("name" in df_features.keys()):
        print("merging based on name and sequence")
        df_merged = pd.merge(df_prop, df_features, right_on=["name", "sequence"], left_on=["name", "sequence"], how="left")

    elif ("name" in df_prop.keys()) and ("name" in df_features.keys()):
        print("merging based on name only")
        df_merged = pd.merge(df_prop, df_features, right_on=["name"], left_on=["name"], how="inner")
    else:
        print("merging based on sequence only")
        df_merged = pd.merge(df_prop, df_features, right_on=["sequence"], left_on=["sequence"], how="right")


    print(df_merged[df_merged.isna().any(axis=1)])

    assert len(df_merged[df_merged.isna().any(axis=1)]) == 0

    print(df_merged.keys())
    
    print("Merged dataframe fully filled")

    if only_designs:
        print("Returning only design PF")
        df_merged = df_merged.query("topology in ['HHH', 'HEEH', 'EHEE', 'EEHEE']").drop_duplicates(["name", "topology"]).reset_index(drop=True)
        
    print("Removing sequences with inconsistent dssp")
        
    df_merged = df_merged[df_merged.apply(lambda x: check_dssp_topology_consistency(x["dssp"], x["topology"]), axis=1)].reset_index(drop=True)
    
    return df_merged

    
#def process_local_features(df, extract_from, prop_save_as, output_dir=None, calculate_mean=True, calculate_max=True, calculate_min=True, calculate_iqr=True, calculate_top3_mean=True, calculate_bottom3_mean=True):
#    l = []
#    for pf in set(df.PF):
#        for i, char in enumerate(pf):
#            l.append(f"{char}{i+1}_{prop_save_as}_mean")
#            l.append(f"{char}{i+1}_{prop_save_as}_max")
#            l.append(f"{char}{i+1}_{prop_save_as}_min")
#            l.append(f"{char}{i+1}_{prop_save_as}_iqr")
#            l.append(f"{char}{i+1}_{prop_save_as}_top3_mean")
#            l.append(f"{char}{i+1}_{prop_save_as}_bottom3_mean")
#
#    l = sorted(set(l))
#    
#    df[l] = np.nan
#    
#    df = df.apply(lambda x: get_properties(x, extract_from, prop_save_as, calculate_mean, calculate_max, calculate_min, calculate_iqr, calculate_top3_mean, calculate_bottom3_mean), axis=1)
#    
#    if not os.path.exists(output_dir):
#        os.makedirs(output_dir)
#        
#    l = []
#        
#    for pf in set(df.PF):
#    
#        df_subset = df.query(f"PF == '{pf}'").reset_index(drop=True).copy()
#        df_subset = df_subset.loc[:, ~df_subset.isna().all(axis=0)].reset_index(drop=True)
#
#        print(f"{pf} {len(df_subset)} entries processed")
#        
#        if output_dir is not None:
#
#            df_subset.to_json(f"{output_dir}/{pf}_{prop_save_as}.json")
#
#            print(f"{output_dir}/{pf}_{prop_save_as}.json saved successfully!")
#            
#        else:
#            l.append(df_subset)
#            
#    if output_dir is None:
#        return l


def process_local_features(df, extract_from_list, prop_save_as_list=None, output_dir=None, calculate_mean=True, calculate_max=True, calculate_min=True, calculate_iqr=True, calculate_top3_mean=True, calculate_bottom3_mean=True):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        

    pf_dataframes = {}

    if prop_save_as_list is None:
        prop_save_as_list = extract_from_list

    stats = ["mean", "max", "min", "iqr", "top3_mean", "bottom3_mean"]
    for extract_from, prop_save_as in zip(extract_from_list, prop_save_as_list):
        l = []
        for stat in stats:
            l.append(f"{prop_save_as}_global_{stat}")
            l.append(f"{prop_save_as}_ss_{stat}")
        for topology in set(df.topology):
            for i, char in enumerate(topology):
                l.append(f"{char}{i+1}_{prop_save_as}_mean")
                l.append(f"{char}{i+1}_{prop_save_as}_max")
                l.append(f"{char}{i+1}_{prop_save_as}_min")
                l.append(f"{char}{i+1}_{prop_save_as}_iqr")
                l.append(f"{char}{i+1}_{prop_save_as}_top3_mean")
                l.append(f"{char}{i+1}_{prop_save_as}_bottom3_mean")

        l = sorted(set(l))

        df[l] = np.nan

        df = df.apply(lambda x: get_properties(x, extract_from, prop_save_as, calculate_mean, calculate_max, calculate_min, calculate_iqr, calculate_top3_mean, calculate_bottom3_mean), axis=1)

    for pf in set(df.PF):
        df_subset = df.query(f"PF  == '{pf}'").reset_index(drop=True).copy()
        df_subset = df_subset.loc[:, ~df_subset.isna().all(axis=0)].reset_index(drop=True)

        print(f"{pf} {len(df_subset)} entries processed")

        if output_dir is not None:
            df_subset.to_json(f"{output_dir}/{pf}.json")
            print(f"{output_dir}/{pf}.json saved successfully!")

        pf_dataframes[pf] = df_subset

    if output_dir is None:
        return pf_dataframes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process local features.')
    parser.add_argument('--features_path', type=str, help='Path to the features file')
    parser.add_argument('--properties_path', type=str, help='Path to the properties file')
    parser.add_argument('--output_dir', type=str, help='Directory to save the output files')
    parser.add_argument('--extract_from', nargs='+', type=str, help='Property to extract data from (e.g., predicted_z_scores)')
    parser.add_argument('--prop_save_as', default=None, nargs='+', type=str, help='Name to save the property as (e.g., adopt)')
#    parser.add_argument('--topology', default='PF', type=str, help='Name of the column with topology information (e.g., PF with info HHH or HEEH)')
    parser.add_argument('--all', default=True, action='store_true', help='Calculate all features (mean, max, min, iqr, top3_mean, bottom3_mean)')
    parser.add_argument('--mean', action='store_true', help='Calculate mean feature')
    parser.add_argument('--max', action='store_true', help='Calculate max feature')
    parser.add_argument('--min', action='store_true', help='Calculate min feature')
    parser.add_argument('--iqr', action='store_true', help='Calculate IQR feature')
    parser.add_argument('--top3_mean', action='store_true', help='Calculate mean of top 3 highest values')
    parser.add_argument('--bottom3_mean', action='store_true', help='Calculate mean of bottom 3 lowest values')

    args = parser.parse_args()

    if args.all:
        args.mean = args.max = args.min = args.iqr = args.top3_mean = args.bottom3_mean = True

    if args.extract_from:
        extract_from_list = args.extract_from
        prop_save_as_list = args.prop_save_as
    if args.prop_save_as is not None:
        if len(prop_save_as_list) != len(extract_from_list):
            parser.error("Both --extract_from and --prop_save_as must have the same length.")

    df_features_path = args.features_path
    df_properties_path = args.properties_path
    output_dir = args.output_dir

    df_merged = load_and_merge_dataframes(df_features_path, df_properties_path)

    process_local_features(df=df_merged, extract_from_list=extract_from_list, prop_save_as_list=prop_save_as_list, output_dir=output_dir, calculate_mean=args.mean, calculate_max=args.max, calculate_min=args.min, calculate_iqr=args.iqr, calculate_top3_mean=args.top3_mean, calculate_bottom3_mean=args.bottom3_mean)

#if __name__ == "__main__":
#    parser = argparse.ArgumentParser(description='Process local features.')
#    parser.add_argument('--features_path', type=str, help='Path to the features file')
#    parser.add_argument('--properties_path', type=str, help='Path to the properties file')
#    parser.add_argument('--output_dir', type=str, help='Directory to save the output files')
#    parser.add_argument('--extract_from', type=str, help='Property to extract data from (e.g., predicted_z_scores)')
#    parser.add_argument('--prop_save_as', type=str, help='Name to save the property as (e.g., adopt)')
#    parser.add_argument('--all', action='store_true', help='Calculate all features (mean, max, min, iqr, top3_mean, bottom3_mean)')
#    parser.add_argument('--mean', action='store_true', help='Calculate mean feature')
#    parser.add_argument('--max', action='store_true', help='Calculate max feature')
#    parser.add_argument('--min', action='store_true', help='Calculate min feature')
#    parser.add_argument('--iqr', action='store_true', help='Calculate IQR feature')
#    parser.add_argument('--top3_mean', action='store_true', help='Calculate mean of top 3 highest values')
#    parser.add_argument('--bottom3_mean', action='store_true', help='Calculate mean of bottom 3 lowest values')
#    
#    args = parser.parse_args()
#
#    if args.all:
#        args.mean = args.max = args.min = args.iqr = args.top3_mean = args.bottom3_mean = True
#
#    df_features_path = args.features_path
#    df_properties_path = args.properties_path
#    output_dir = args.output_dir
#    extract_from = args.extract_from
#    prop_save_as = args.prop_save_as
#
#    df_merged = load_and_merge_dataframes(df_features_path, df_properties_path)
#    process_local_features(df=df_merged, extract_from=extract_from, prop_save_as=prop_save_as, output_dir=output_dir, calculate_mean=args.mean, calculate_max=args.max, calculate_min=args.min, calculate_iqr=args.iqr, calculate_top3_mean=args.top3_mean, calculate_bottom3_mean=args.bottom3_mean)
