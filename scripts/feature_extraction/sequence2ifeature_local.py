import argparse
import pandas as pd
import iFeatureOmegaCLI
import re
import pickle as pk
import sys
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

    print(f"Topology mismatch. Early termination.")
    sys.exit(1)
    

def get_prop(prop, segment):
    i, j = segment
    return prop[i:j]

def get_substrings(x):
    name = x["name"]
    topology = x["topology"]
    sequence = x["sequence"]
    dssp = x["dssp"]
    segments = find_topn_segments(dssp, topology)
    l = []
    for i, char in enumerate(topology):
        new_name = f"{name}_{char}{i+1}"
        substr = get_prop(sequence, segments[i])
        l.append([new_name, substr])
    return l

def df_to_fasta(df, output, keys=["name", "sequence"], unique=False, l=[]):
    
    if unique:
        l = l
        with open(output, "w") as f:
            for _, line in df.iterrows(): 
                if line[keys[0]] not in l:
                    l.append(line[keys[0]])
                    f.write(f'>{line[keys[0]]}\n')
                    f.write(f'{line[keys[1]]}\n')
    else:
        with open(output, "w") as f:
            for _, line in df.iterrows():                
                f.write(f'>{line[keys[0]]}\n')
                f.write(f'{line[keys[1]]}\n')


def process_row(df, row):
        
    prefix = row.name.split("_")[-1]
    name = "_".join(row.name.split("_")[:-1])
    
    #print(f"Processing {name}...")
    
    cols = [f"{prefix}_{col}" for col in row.index]
    
    values = row.values
    
    df.loc[df.name == name, cols] = values 


def main(fasta, dssp, output):
    df_dssp = pd.read_json(dssp)
    protein = iFeatureOmegaCLI.iProtein(fasta)

    l = df_dssp.apply(lambda x: get_substrings(x), axis=1)
    df_substrings = pd.DataFrame([i for j in l for i in j], columns=["name", "sequence"])
    df_to_fasta(df_substrings, output=output + ".substrings.fasta", keys=["name", "sequence"], unique=True, l=[])
    

    protein = iFeatureOmegaCLI.iProtein(output + ".substrings.fasta")

    descriptors = ["AAC", "GAAC", "CTDC"]
    ls = []
    for descriptor in descriptors:
        protein.get_descriptor(descriptor)
        suffix_ = "_frac" if "type 1" in descriptor else "_count" if "type 2" in descriptor else ""
        protein.encodings.columns = [col + suffix_ for col in protein.encodings.columns]
        ls.append(protein.encodings)
    
    df_substrings_features = pd.concat(ls, axis=1)
    #print(f"Features DataFrame shape: {df_substrings_features.shape}")

    #print("processing names...")

     
    names = set(["_".join(i.split("_")[:-1]) for i in df_substrings_features.index])
    
    #print("processing columns...")
    cols = df_substrings_features.columns.tolist()
    #print("processing prefix...")
    prefixes = set([i.split("_")[-1] for i in df_substrings_features.index])
    #print("processing new columns...")
    new_cols = [f"{prefix}_{col}" for prefix in prefixes for col in cols]
    
    #print("generating new dataframe...")
    new_df = pd.DataFrame(columns=new_cols)
    new_df["name"] = pd.Series(list(names))
    
    new_df = pd.merge(new_df, df_dssp[["name", "sequence", "topology"]], on="name", how="left")


    df_substrings_features.apply(lambda x: process_row(new_df, x), axis=1)

    new_df.to_json(output)
    print(f"Output saved to {output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process sequence data and extract features.')
    parser.add_argument('--fasta', type=str, required=True, help="Fasta with sequence to be processed")
    parser.add_argument('--dssp', type=str, required=True, help='Path to input JSON file containing sequence and dssp')
    parser.add_argument('--output', type=str, required=True, help='Path to save the processed data and extracted features in JSON format.')
    args = parser.parse_args()
    main(args.fasta, args.dssp, args.output)

