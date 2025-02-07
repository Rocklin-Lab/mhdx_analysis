import argparse
import Bio
from Bio.PDB import PDBParser, PPBuilder, calc_dihedral
import pandas as pd
from numpy import pi
import os
import pyrosetta
from pyrosetta.toolbox import cleanATOM

# Initialize PyRosetta
pyrosetta.init()

penalty_table_path = '/projects/p30802/grocklin/score_monomeric_designs/scripts/penalty_table'

penalties = {}
aas = 'ACDEFGHIKLMNPQRSTVWY'
with open(penalty_table_path) as file:
    for line in file.readlines():
        abego_aa = tuple(line.split()[0:2])
        penalties[abego_aa] = float(line.split()[-2])


def phi_psi_omega_to_abego(phi, psi, omega):
    if psi is None: return 'O'
    if omega is None: omega = 180
    if phi is None: phi = 90
    phi = 180 * phi / pi
    psi = 180 * psi / pi
    omega = 180 * omega / pi

    if abs(omega) < 90:
        return 'O'
    elif phi > 0:
        if -100.0 <= psi < 100:
            return 'G'
        else:
            return 'E'
    else:
        if -75.0 <= psi < 50:
            return 'A'
        else:
            return 'B'
    return 'X'

def abego_string(phi_psi_omega):
    out = ''
    for x in phi_psi_omega:
        out += phi_psi_omega_to_abego(x[0], x[1], x[2])
    return out


def penalty_to_numeral(value, cutoffs=[-3, -1, -0.5, 0, 0.2, 0.6, 1.0, 1.5]):
    # Returns the numeral based on the given value and cutoffs
    return str(len([x for x in cutoffs if x < value]) + 1)


def file_to_abego(file):
    
    
    phi_psi = []
    nres=[]
    seqs = []
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(file,file)
    for chain in structure:
        polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
        for polypeptide in polypeptides:
            nres.append(len(polypeptide))
            seqs.append(polypeptide.get_sequence())
            if len(nres) > 1:
                nres[-1] = nres[-1] + nres[-2]
            phi_psi += polypeptide.get_phi_psi_list()
    residues = [res for res in structure.get_residues()]
    phi_psi_omega = []
    #print fn, len(phi_psi), len(residues), nres
    for i in range(len(residues)-1):
        if i +1 in nres:
            omega = None
        else:
            a1 = residues[i]['CA'].get_vector()
            a2 = residues[i]['C'].get_vector()
            a3 = residues[i+1]['N'].get_vector()
            a4 = residues[i+1]['CA'].get_vector()
            omega = Bio.PDB.calc_dihedral(a1,a2,a3,a4)
        phi_psi_omega.append((phi_psi[i][0], phi_psi[i][1], omega))
    phi_psi_omega.append((phi_psi[-1][0], phi_psi[-1][1], None))

    temp_abego_string = abego_string(phi_psi_omega)
    my_seq_string = ''
    my_abego_string = ''
    #print abego_string
    #print
    for s in seqs:
        my_seq_string += str(s)
        my_abego_string += temp_abego_string[0:len(str(s))]
        if str(s) != str(seqs[-1]):
            my_seq_string += ' '
            my_abego_string += ' '
            temp_abego_string = temp_abego_string[len(str(s)):]
            #print
            #print abego_string


    my_penalty_string = '55'
    scores = [0,0]
    best_aa_string = my_seq_string[0:2]
    best_aa_scores = '55'
    for i in range(2,len(my_abego_string)-2):
        if ' ' not in my_abego_string[i-2:i+3]:
            if (my_abego_string[i-1:i+2], my_seq_string[i]) in penalties:
                term = penalties[(my_abego_string[i-1:i+2], my_seq_string[i])]
                aa_penalties=[(penalties[my_abego_string[i-1:i+2], aa], aa) for aa in aas]
                aa_penalties.sort()
                #print aa_penalties
                best_aa = aa_penalties[-1][-1]
                best_aa_string += best_aa
                best_aa_scores += penalty_to_numeral(  penalties[(my_abego_string[i-1:i+2], best_aa)] )
                scores.append(-10*term)
                my_penalty_string += penalty_to_numeral(term)

        elif my_abego_string[i] != ' ':
            scores.append(0)
            my_penalty_string += '5'
            best_aa_string += my_seq_string[i]
            best_aa_scores += '5'
        else:
            my_penalty_string += ' '
            best_aa_string += ' '
            best_aa_scores += ' '

                
    return my_abego_string, my_penalty_string


def get_dssp_pyrosetta(f):
    
    # Create a Pose object from a PDB file
    pose = pyrosetta.pose_from_pdb(f)

    # Create a secondary structure object
    secstruct = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose)

    # Extract secondary structure information
    secstruct.insert_ss_into_pose(pose)
    
    # Iterate through the pose and print residue number and corresponding secondary structure
    ss = ''
    for i in range(1, pose.total_residue() + 1):
        ss += pose.secstruct(i)
        
    return ss


def get_sequence_from_pose(f):
    
     
    # Create a Pose object from a PDB file
    pose = pyrosetta.pose_from_pdb(f)

    return pose.sequence()

def main():
    parser = argparse.ArgumentParser(description='Process a PDB file to calculate ABEGO string.')
    parser.add_argument('--input', help='Path to the input PDB file.')
    parser.add_argument('--output', help='Path to the output JSON file.')
    parser.add_argument('--topology', default=None, help="String defining (expected) topology, e.g., HHH")
    args = parser.parse_args()

    sequence = get_sequence_from_pose(args.input)
    dssp = get_dssp_pyrosetta(args.input)
    abego, abego_penalty = file_to_abego(args.input)

    # Prepare data for DataFrame
    data = {
        'name': os.path.basename(args.input),
        'sequence': sequence,
        'dssp': dssp,
        'abego': abego,
        'abego_penalty': abego_penalty,
        'topology': args.topology
    }

    # Save to JSON
    df = pd.DataFrame([data])
    df.to_json(args.output, orient='records', indent=4)


if __name__ == "__main__":
    main()

