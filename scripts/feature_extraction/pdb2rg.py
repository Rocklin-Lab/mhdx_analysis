import argparse
import pandas as pd
import pyrosetta
from pyrosetta import pose_from_pdb
import numpy as np

# Initialize PyRosetta
pyrosetta.init()

def radius_of_gyration(pose):
    coordinates = []
    for residue in pose.residues:
        for atom_index in range(1, residue.natoms() + 1):
            atom_xyz = residue.xyz(atom_index)
            coordinates.append([atom_xyz.x, atom_xyz.y, atom_xyz.z])
    
    coordinates = np.array(coordinates)
    center_of_mass = np.mean(coordinates, axis=0)
    distances = np.linalg.norm(coordinates - center_of_mass, axis=1)
    rg = np.sqrt(np.mean(distances**2))
    return rg

def get_sequence(pose):
    sequence = ''
    for residue in pose.residues:
        sequence += residue.name1()
    return sequence

def main(pdb_file, output_file):
    protein_pose = pose_from_pdb(pdb_file)
    rg = radius_of_gyration(protein_pose)
    sequence = get_sequence(protein_pose)
    data = {'sequence': [sequence], 'rg': [rg]}
    df = pd.DataFrame(data)
    df.to_json(output_file, orient='records')
    print(f"Data saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute the radius of gyration and save with sequence in JSON.")
    parser.add_argument("-i", "--input", help="Path to the PDB file.")
    parser.add_argument("-o", "--output", help="Path to the output JSON file.")
    args = parser.parse_args()
    
    main(args.input, args.output)

