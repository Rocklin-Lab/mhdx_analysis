import argparse
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pyrosetta
from pyrosetta import rosetta
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet

# Initialize PyRosetta
pyrosetta.init("-mute all")


def get_hb_df(pdb_file_path):
    """Extract hydrogen bond information from a PDB file."""
    pose = rosetta.core.import_pose.pose_from_file(pdb_file_path)
    hbond_set = HBondSet(pose, bb_only=False)

    hbonds_data = []
    for hb in hbond_set.hbonds():
        residue_acc = pose.residue(hb.acc_res())
        residue_don = pose.residue(hb.don_res())
        hbonds_data.append([
            hb.acc_res(), hb.don_res(),
            residue_acc.name1(), residue_don.name1(),
            residue_acc.atom_name(hb.acc_atm()), residue_don.atom_name(hb.don_hatm()),
            hb.acc_atm_is_protein_backbone(), hb.don_hatm_is_backbone(),
            hb.energy()
        ])

    columns = ["res_num_acc", "res_num_don", "res_name_acc", "res_name_don", "atom_acc", "atom_don", "acc_is_bb", "don_is_bb", "hb_energy"]
    return pd.DataFrame(hbonds_data, columns=columns)

def get_protein_sequence(pdb_file_path):
    """Extract the protein sequence from a PDB file."""
    pose = rosetta.core.import_pose.pose_from_file(pdb_file_path)
    return pose.sequence()

def get_hb_df(pdb_file_path):
    """Extract hydrogen bond information from a PDB file."""
    pose = rosetta.core.import_pose.pose_from_file(pdb_file_path)
    hbond_set = HBondSet(pose, bb_only=False)

    hbonds_data = []
    for hb in hbond_set.hbonds():
        residue_acc = pose.residue(hb.acc_res())
        residue_don = pose.residue(hb.don_res())
        hbonds_data.append([
            hb.acc_res(), hb.don_res(),
            residue_acc.name1(), residue_don.name1(),
            residue_acc.atom_name(hb.acc_atm()), residue_don.atom_name(hb.don_hatm()),
            hb.acc_atm_is_protein_backbone(), hb.don_hatm_is_backbone(),
            hb.energy()
        ])

    columns = ["res_num_acc", "res_num_don", "res_name_acc", "res_name_don", "atom_acc", "atom_don", "acc_is_bb", "don_is_bb", "hb_energy"]
    return pd.DataFrame(hbonds_data, columns=columns)

def get_protein_sequence(pdb_file_path):
    """Extract the protein sequence from a PDB file."""
    pose = rosetta.core.import_pose.pose_from_file(pdb_file_path)
    return pose.sequence()




def analyze_hbonds(input_file, burial_file=None, output_file=None):
    """Analyze hydrogen bonds in a PDB file and save the results to a JSON file."""
    name = input_file.split("/")[-1]
    sequence = get_protein_sequence(input_file)
    df = get_hb_df(input_file)

    # Your analysis and data aggregation logic here
    n_hb_bb_bb_don = len(set(df.query("(don_is_bb == True) & (acc_is_bb == True) & (res_num_don != 1) & (res_num_don != 2)")["res_num_don"].values))
    n_hb_bb_bb_acc = len(set(df.query("(don_is_bb == True) & (acc_is_bb == True) & (res_num_don != 1) & (res_num_don != 2)")["res_num_acc"].values))
    n_hb_bb_sc_don = len(set(df.query("(don_is_bb == True) & (acc_is_bb == False) & (res_num_don != 1) & (res_num_don != 2)")["res_num_don"].values))
    n_hb_bb_sc_acc = len(set(df.query("(don_is_bb == True) & (acc_is_bb == False) & (res_num_don != 1) & (res_num_don != 2)")["res_num_acc"].values))
    n_hb_sc_sc_don = len(set(df.query("(don_is_bb == False) & (acc_is_bb == False) & (res_num_don != 1) & (res_num_don != 2)")["res_num_don"].values))
    n_hb_sc_sc_acc = len(set(df.query("(don_is_bb == False) & (acc_is_bb == False) & (res_num_don != 1) & (res_num_don != 2)")["res_num_acc"].values))

    n_hb_all = len(df)

    array_bb_bb_don = df.query("(don_is_bb == True) & (acc_is_bb == True)")["res_num_don"].values.tolist()
    array_bb_bb_acc = df.query("(don_is_bb == True) & (acc_is_bb == True)")["res_num_acc"].values.tolist()
    array_bb_sc_don = df.query("(don_is_bb == True) & (acc_is_bb == False)")["res_num_don"].values.tolist()
    array_bb_sc_acc = df.query("(don_is_bb == True) & (acc_is_bb == False)")["res_num_acc"].values.tolist()
    array_sc_sc_don = df.query("(don_is_bb == False) & (acc_is_bb == False)")["res_num_don"].values.tolist()
    array_sc_sc_acc = df.query("(don_is_bb == False) & (acc_is_bb == False)")["res_num_acc"].values.tolist()

    array_hb_energy = df["hb_energy"].values.tolist()
    array_acc = df["res_num_acc"].values.tolist()
    array_don = df["res_num_don"].values.tolist()

    n_hb_bb_all = len([i for i in set(array_bb_bb_don + array_bb_sc_don) if (i!=1 and i!=2)])

    l = []
    l.append([name, sequence, n_hb_all, n_hb_bb_all, n_hb_bb_bb_don,n_hb_bb_bb_acc,n_hb_sc_sc_don,n_hb_sc_sc_acc, 
              n_hb_bb_sc_don, n_hb_bb_sc_acc, array_bb_bb_don, array_bb_bb_acc, array_bb_sc_don, 
              array_bb_sc_acc, array_sc_sc_don, array_sc_sc_acc, array_acc,array_don, array_hb_energy])

    result_df = pd.DataFrame(l, columns=["name", "sequence", "n_hb_all", "n_hb_bb_all", "n_hb_bb_bb_don", 
                                         "n_hb_bb_bb_acc", "n_hb_sc_sc_don", "n_hb_sc_sc_acc", 
                                         "n_hb_bb_sc_don", "n_hb_bb_sc_acc", "array_bb_bb_don",
                                         "array_bb_bb_acc", "array_bb_sc_don", "array_bb_sc_acc", 
                                         "array_sc_sc_don", "array_sc_sc_acc", 
                                         "array_acc", "array_don", "array_hb_energy"])
    
    
    return result_df


def sparse_to_full_array(row, full_col, sparse_col):
    
    sequence = row['sequence']
    full_values = row[full_col]
    sparse_values = row[sparse_col]
    
    # Initialize the third array with zeros
    array = np.zeros(len(sequence))

    # Iterate over the values and indexes arrays
    for value, index in zip(full_values, sparse_values):
        # Adjusting for 1-indexed array, and check if index is within bounds
        if 0 < index <= len(sequence):
            # Subtract 1 from index to convert from 1-indexed to 0-indexed
            adjusted_index = index - 1
            # Add the value from the values array to the corresponding position in the third array
            array[adjusted_index] += value

    return array


def analyze_hbonds(input_file, output_file=None):
    """Analyze hydrogen bonds in a PDB file and save the results to a JSON file."""
    name = input_file.split("/")[-1]
    sequence = get_protein_sequence(input_file)
    df = get_hb_df(input_file)

    # Your analysis and data aggregation logic here
    n_hb_bb_bb_don = len(set(df.query("(don_is_bb == True) & (acc_is_bb == True) & (res_num_don != 1) & (res_num_don != 2)")["res_num_don"].values))
    n_hb_bb_bb_acc = len(set(df.query("(don_is_bb == True) & (acc_is_bb == True) & (res_num_don != 1) & (res_num_don != 2)")["res_num_acc"].values))
    n_hb_bb_sc_don = len(set(df.query("(don_is_bb == True) & (acc_is_bb == False) & (res_num_don != 1) & (res_num_don != 2)")["res_num_don"].values))
    n_hb_bb_sc_acc = len(set(df.query("(don_is_bb == True) & (acc_is_bb == False) & (res_num_don != 1) & (res_num_don != 2)")["res_num_acc"].values))
    n_hb_sc_sc_don = len(set(df.query("(don_is_bb == False) & (acc_is_bb == False) & (res_num_don != 1) & (res_num_don != 2)")["res_num_don"].values))
    n_hb_sc_sc_acc = len(set(df.query("(don_is_bb == False) & (acc_is_bb == False) & (res_num_don != 1) & (res_num_don != 2)")["res_num_acc"].values))

    n_hb_all = len(df)

    array_bb_bb_don = df.query("(don_is_bb == True) & (acc_is_bb == True)")["res_num_don"].values.tolist()
    array_bb_bb_acc = df.query("(don_is_bb == True) & (acc_is_bb == True)")["res_num_acc"].values.tolist()
    array_bb_sc_don = df.query("(don_is_bb == True) & (acc_is_bb == False)")["res_num_don"].values.tolist()
    array_bb_sc_acc = df.query("(don_is_bb == True) & (acc_is_bb == False)")["res_num_acc"].values.tolist()
    array_sc_sc_don = df.query("(don_is_bb == False) & (acc_is_bb == False)")["res_num_don"].values.tolist()
    array_sc_sc_acc = df.query("(don_is_bb == False) & (acc_is_bb == False)")["res_num_acc"].values.tolist()

    array_hb_energy = df["hb_energy"].values.tolist()
    array_acc = df["res_num_acc"].values.tolist()
    array_don = df["res_num_don"].values.tolist()

    l = []
    l.append([name, sequence, n_hb_all, n_hb_bb_bb_don,n_hb_bb_bb_acc,n_hb_sc_sc_don,n_hb_sc_sc_acc, 
              n_hb_bb_sc_don, n_hb_bb_sc_acc, array_bb_bb_don, array_bb_bb_acc, array_bb_sc_don, 
              array_bb_sc_acc, array_sc_sc_don, array_sc_sc_acc, array_acc,array_don, array_hb_energy])

    result_df = pd.DataFrame(l, columns=["name", "sequence", "n_hb_all", "n_hb_bb_bb_don", 
                                         "n_hb_bb_bb_acc", "n_hb_sc_sc_don", "n_hb_sc_sc_acc", 
                                         "n_hb_bb_sc_don", "n_hb_bb_sc_acc", "array_bb_bb_don",
                                         "array_bb_bb_acc", "array_bb_sc_don", "array_bb_sc_acc", 
                                         "array_sc_sc_don", "array_sc_sc_acc", 
                                         "array_acc", "array_don", "array_hb_energy"])
    

    cols = [i for i in result_df if ('array_' in i) and ('array_hb_energy' not in i)]
    for col in cols:
        new_col = col.replace('array_', 'hbond_')
    
        result_df[new_col] = result_df.apply(
            lambda x: sparse_to_full_array(
                x, 'array_hb_energy', col
            ), axis=1).tolist()

    result_df['hb_total_energy'] = result_df.apply(lambda x: np.array(x['hbond_acc']) + np.array(x['hbond_don']), axis=1)

    if output_file is not None:
        result_df.to_json(output_file)
    else:
        return result_df



def main():
    parser = argparse.ArgumentParser(description="Analyze hydrogen bonds in a PDB file and output the results to a JSON file.")
    parser.add_argument('-i', '--input', help="Input PDB file path", required=True)
    parser.add_argument('-o', '--output', help="Output JSON file path", required=True)
    
    args = parser.parse_args()
    
    analyze_hbonds(args.input, args.output)

if __name__ == "__main__":
    main()

