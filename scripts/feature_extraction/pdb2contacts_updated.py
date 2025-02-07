import numpy as np
import pandas as pd
from Bio.PDB.Polypeptide import PPBuilder
import Bio.PDB
import argparse

def vtoa(v):
    return np.array([float(x) for x in v])

def unit(v):
    return v / (np.dot(v,v)**0.5)

def norm(v):
    return np.dot(v,v) ** 0.5

def pdb2sasa_burial(file_loc):

    # define basic paramters
    d_param = 9
    d_param_sc = 9
    
    # import pdb file
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(file_loc,file_loc)
    structure = structure[0]
    residues = [res for res in structure.get_residues() if res.id[0] == " "]
    ppb=PPBuilder()
    aa_list = list(ppb.build_peptides(structure)[0].get_sequence())
    hbond_list = list()
    
    
    burial_list,burial_hydrophobic_list,sc_contact_list, local_sc_contact_list, non_local_sc_contact_list,arom_contact_list,local_aromatic_contact_list, non_local_aromatic_contact_list,non_favorable_contacts_DE_list, local_non_favorable_contacts_DE_list, non_local_non_favorable_contacts_DE_list,non_favorable_contacts_RK_list, local_non_favorable_contacts_RK_list, non_local_non_favorable_contacts_RK_list,favorable_contact_list,local_favorable_contacts_list, non_local_favorable_contacts_list=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    
    
    # calculate burial and side chain contact
    for i, res in enumerate(residues):
        burial,burial_hydrophobic,sc_contact,local_sc_contact, non_local_sc_contact,arom_contact, local_aromatic_contact, non_local_aromatic_contact,non_favorable_contacts_DE, local_non_favorable_contacts_DE, non_local_non_favorable_contacts_DE,non_favorable_contacts_RK, local_non_favorable_contacts_RK,non_local_non_favorable_contacts_RK,favorable_contacts, local_favorable_contacts, non_local_favorable_contacts = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0

        # calculate vector from CA to CB of the amino acid of interest
        if res.get_resname() == 'GLY':
            try:
                start = res['2HA']
            except:
                start = res['HA2']
        else:
            start = res['CB']
        my_cb_vector = unit(vtoa(start.get_vector() - res['CA'].get_vector()))
        
        
        for j, res2 in enumerate(residues):
            if res2 != res:
                
                sequence_difference = abs(i-j)

                # calculate burial: # of Cα atoms in a cone projecting out 9 Å away from the Cβ atom on residue X in the direction of the residue X Cα-Cβ vector.
                to_res2 = vtoa(res2['CA'].get_vector() - start.get_vector())
                dist_term = 1.0 / (1.0 + (np.exp(norm(to_res2) - d_param)))
                angle_term = (0.5 + np.dot(unit(to_res2), my_cb_vector))
                if angle_term < 0: angle_term = 0
                angle_term = angle_term ** 2.0
                burial += dist_term * angle_term / 2.25
                
                
#                if aa_list[j] in 'AILVPM' and aa_list[i] in 'AILVPM':
                if aa_list[j] in 'AILVMFWY' and aa_list[i] in 'AILVMFWY':
                    to_res2 = vtoa(res2['CA'].get_vector() - start.get_vector())
                    dist_term = 1.0 / (1.0 + (np.exp(norm(to_res2) - d_param)))
                    angle_term = (0.5 + np.dot(unit(to_res2), my_cb_vector))
                    if angle_term < 0: angle_term = 0
                    angle_term = angle_term ** 2.0
                    burial_hydrophobic += dist_term * angle_term / 2.25
                    
                
                # calculate (total) side chain contact: # of Cβ atoms in a cone projecting out 9 Å away from the Cβ atom on residue X in the direction of the residue X Cα-Cβ vector.
                if aa_list[j] == 'G':
                    try:
                        to_res2_sc = vtoa(res2['2HA'].get_vector() - start.get_vector())
                    except:
                        to_res2_sc = vtoa(res2['HA2'].get_vector() - start.get_vector())
                else:
                    to_res2_sc = vtoa(res2['CB'].get_vector() - start.get_vector())
                dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                if angle_term_sc < 0: angle_term_sc = 0
                angle_term_sc = angle_term_sc ** 2.0
                sc_contact += dist_term_sc * angle_term_sc / 2.25

                if sequence_difference <= 4:
                    if aa_list[j] == 'G':
                        try:
                            to_res2_sc = vtoa(res2['2HA'].get_vector() - start.get_vector())
                        except:
                            to_res2_sc = vtoa(res2['HA2'].get_vector() - start.get_vector())
                    else:
                        to_res2_sc = vtoa(res2['CB'].get_vector() - start.get_vector())
                    dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                    angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                    if angle_term_sc < 0: angle_term_sc = 0
                    angle_term_sc = angle_term_sc ** 2.0                    
                    local_sc_contact += dist_term_sc * angle_term_sc / 2.25

                if sequence_difference > 4:
                    non_local_sc_contact += dist_term_sc * angle_term_sc / 2.25
          
                
                if aa_list[j] in 'FYW' and aa_list[i] in 'FYW':
                    to_res2_sc = vtoa(res2['CE2'].get_vector() - start.get_vector())
                    dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                    angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                    if angle_term_sc < 0: angle_term_sc = 0
                    angle_term_sc = angle_term_sc ** 2.0
                    arom_contact+= dist_term_sc * angle_term_sc / 2.25
                    
                if sequence_difference <= 4:
                    if aa_list[j] in 'FYW' and aa_list[i] in 'FYW':
                        to_res2_sc = vtoa(res2['CE2'].get_vector() - start.get_vector())
                        dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                        angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                        if angle_term_sc < 0: angle_term_sc = 0
                        angle_term_sc = angle_term_sc ** 2.0                    
                        local_aromatic_contact += dist_term_sc * angle_term_sc / 2.25
                if sequence_difference > 4:
                    if aa_list[j] in 'FYW' and aa_list[i] in 'FYW':
                        to_res2_sc = vtoa(res2['CE2'].get_vector() - start.get_vector())
                        dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                        angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                        if angle_term_sc < 0: angle_term_sc = 0
                        angle_term_sc = angle_term_sc ** 2.0                                       
                        non_local_aromatic_contact += dist_term_sc * angle_term_sc / 2.25
                    
                    
                # calculate acidic side chain contact: # of OE1 atoms of Glu + OD1 atoms of Asp in a cone projecting out 9 Å away from the Cβ atom on residue X in the direction of the residue X Cα-Cβ vector.    
                if aa_list[j] in 'DE' and aa_list[i] in 'DE':
                    try:
                        to_res2_sc = vtoa(res2['OE1'].get_vector() - start.get_vector())
                    except:
                        to_res2_sc = vtoa(res2['OD1'].get_vector() - start.get_vector())
                    dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                    angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                    if angle_term_sc < 0: angle_term_sc = 0
                    angle_term_sc = angle_term_sc ** 2.0
                    non_favorable_contacts_DE+= dist_term_sc * angle_term_sc / 2.25
                    
                if sequence_difference <= 4:
                    if aa_list[j] in 'DE' and aa_list[i] in 'DE':
                        try:
                            to_res2_sc = vtoa(res2['OE1'].get_vector() - start.get_vector())
                        except:
                            to_res2_sc = vtoa(res2['OD1'].get_vector() - start.get_vector())
                        dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                        angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                        if angle_term_sc < 0: angle_term_sc = 0
                        angle_term_sc = angle_term_sc ** 2.0                    
                        local_non_favorable_contacts_DE += dist_term_sc * angle_term_sc / 2.25
                if sequence_difference > 4:
                    if aa_list[j] in 'DE' and aa_list[i] in 'DE':
                        try:
                            to_res2_sc = vtoa(res2['OE1'].get_vector() - start.get_vector())
                        except:
                            to_res2_sc = vtoa(res2['OD1'].get_vector() - start.get_vector())
                        dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                        angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                        if angle_term_sc < 0: angle_term_sc = 0
                        angle_term_sc = angle_term_sc ** 2.0                                       
                        non_local_non_favorable_contacts_DE += dist_term_sc * angle_term_sc / 2.25
                    
                    
                if aa_list[j] in 'DE' and aa_list[i] in 'RK':
                    try:
                        to_res2_sc = vtoa(res2['OE1'].get_vector() - start.get_vector())
                    except:
                        to_res2_sc = vtoa(res2['OD1'].get_vector() - start.get_vector())
                    dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                    angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                    if angle_term_sc < 0: angle_term_sc = 0
                    angle_term_sc = angle_term_sc ** 2.0
                    favorable_contacts+= dist_term_sc * angle_term_sc / 2.25
                
                if sequence_difference <= 4:
                    if aa_list[j] in 'DE' and aa_list[i] in 'RK':
                        try:
                            to_res2_sc = vtoa(res2['OE1'].get_vector() - start.get_vector())
                        except:
                            to_res2_sc = vtoa(res2['OD1'].get_vector() - start.get_vector())
                        dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                        angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                        if angle_term_sc < 0: angle_term_sc = 0
                        angle_term_sc = angle_term_sc ** 2.0                   
                        local_favorable_contacts += dist_term_sc * angle_term_sc / 2.25
                if sequence_difference > 4:
                    if aa_list[j] in 'DE' and aa_list[i] in 'RK':
                        try:
                            to_res2_sc = vtoa(res2['OE1'].get_vector() - start.get_vector())
                        except:
                            to_res2_sc = vtoa(res2['OD1'].get_vector() - start.get_vector())
                        dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                        angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                        if angle_term_sc < 0: angle_term_sc = 0
                        angle_term_sc = angle_term_sc ** 2.0                         
                        non_local_favorable_contacts += dist_term_sc * angle_term_sc / 2.25                  
                
                
                # calculate basic side chain contact: # of NZ atoms of Lys + NE atoms of Arg in a cone projecting out 9 Å away from the Cβ atom on residue X in the direction of the residue X Cα-Cβ vector.    
                if aa_list[j] in 'RK' and aa_list[i] in 'DE':
                    try:
                        to_res2_sc = vtoa(res2['NZ'].get_vector() - start.get_vector())
                    except:
                        to_res2_sc = vtoa(res2['NE'].get_vector() - start.get_vector())
                    dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                    angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                    if angle_term_sc < 0: angle_term_sc = 0
                    angle_term_sc = angle_term_sc ** 2.0
                    favorable_contacts += dist_term_sc * angle_term_sc / 2.25 
                    
                if sequence_difference <= 4:
                    if aa_list[j] in 'RK' and aa_list[i] in 'DE':
                        try:
                            to_res2_sc = vtoa(res2['NZ'].get_vector() - start.get_vector())
                        except:
                            to_res2_sc = vtoa(res2['NE'].get_vector() - start.get_vector())
                        dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                        angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                        if angle_term_sc < 0: angle_term_sc = 0
                        angle_term_sc = angle_term_sc ** 2.0                    
                        local_favorable_contacts += dist_term_sc * angle_term_sc / 2.25
                if sequence_difference > 4:
                    if aa_list[j] in 'RK' and aa_list[i] in 'DE':
                        try:
                            to_res2_sc = vtoa(res2['NZ'].get_vector() - start.get_vector())
                        except:
                            to_res2_sc = vtoa(res2['NE'].get_vector() - start.get_vector())
                        dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                        angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                        if angle_term_sc < 0: angle_term_sc = 0
                        angle_term_sc = angle_term_sc ** 2.0                      
                        non_local_favorable_contacts += dist_term_sc * angle_term_sc / 2.25 
                    
                if aa_list[j] in 'RK' and aa_list[i] in 'RK':
                    try:
                        to_res2_sc = vtoa(res2['NZ'].get_vector() - start.get_vector())
                    except:
                        to_res2_sc = vtoa(res2['NE'].get_vector() - start.get_vector())
                    dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                    angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                    if angle_term_sc < 0: angle_term_sc = 0
                    angle_term_sc = angle_term_sc ** 2.0
                    non_favorable_contacts_RK += dist_term_sc * angle_term_sc / 2.25
                    
                if sequence_difference <= 4:
                    if aa_list[j] in 'RK' and aa_list[i] in 'RK':
                        try:
                            to_res2_sc = vtoa(res2['NZ'].get_vector() - start.get_vector())
                        except:
                            to_res2_sc = vtoa(res2['NE'].get_vector() - start.get_vector())
                        dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                        angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                        if angle_term_sc < 0: angle_term_sc = 0
                        angle_term_sc = angle_term_sc ** 2.0                    
                        local_non_favorable_contacts_RK += dist_term_sc * angle_term_sc / 2.25
                if sequence_difference > 4:
                    if aa_list[j] in 'RK' and aa_list[i] in 'RK':
                        try:
                            to_res2_sc = vtoa(res2['NZ'].get_vector() - start.get_vector())
                        except:
                            to_res2_sc = vtoa(res2['NE'].get_vector() - start.get_vector())
                        dist_term_sc = 1.0 / (1.0 + (np.exp(norm(to_res2_sc) - d_param_sc)))
                        angle_term_sc = (0.5 + np.dot(unit(to_res2_sc), my_cb_vector))
                        if angle_term_sc < 0: angle_term_sc = 0
                        angle_term_sc = angle_term_sc ** 2.0                      
                        non_local_non_favorable_contacts_RK += dist_term_sc * angle_term_sc / 2.25

        
        burial_list.append(burial)
        burial_hydrophobic_list.append(burial_hydrophobic)
        sc_contact_list.append(sc_contact)
        local_sc_contact_list.append(local_sc_contact)
        non_local_sc_contact_list.append(non_local_sc_contact)
        arom_contact_list.append(arom_contact)
        local_aromatic_contact_list.append(local_aromatic_contact)
        non_local_aromatic_contact_list.append(non_local_aromatic_contact)
        non_favorable_contacts_DE_list.append(non_favorable_contacts_DE)
        local_non_favorable_contacts_DE_list.append(local_non_favorable_contacts_DE)
        non_local_non_favorable_contacts_DE_list.append(non_local_non_favorable_contacts_DE)
        non_favorable_contacts_RK_list.append(non_favorable_contacts_RK)
        local_non_favorable_contacts_RK_list.append(local_non_favorable_contacts_RK)
        non_local_non_favorable_contacts_RK_list.append(non_local_non_favorable_contacts_RK)
        favorable_contact_list.append(favorable_contacts)
        local_favorable_contacts_list.append(local_favorable_contacts)
        non_local_favorable_contacts_list.append(non_local_favorable_contacts)


    burial_sasa_df = pd.DataFrame({
        'sequence': ["".join(aa_list)],
        'burial': [burial_list],
        'burial_hydrophobic': [burial_hydrophobic_list],
        'sc_contact': [sc_contact_list],
        'local_sc_contact': [local_sc_contact_list],
        'non_local_sc_contact': [non_local_sc_contact_list],
        'arom_contact': [arom_contact_list],
        'local_arom_contact': [local_aromatic_contact_list],
        'non_local_arom_contact': [non_local_aromatic_contact_list],
        'non_favorable_contacts_DE': [non_favorable_contacts_DE_list],
        'local_non_favorable_contacts_DE': [local_non_favorable_contacts_DE_list],
        'non_local_non_favorable_contacts_DE': [non_local_non_favorable_contacts_DE_list],
        'non_favorable_contacts_RK': [non_favorable_contacts_RK_list],
        'local_non_favorable_contacts_RK': [local_non_favorable_contacts_RK_list],
        'non_local_non_favorable_contacts_RK': [non_local_non_favorable_contacts_RK_list],
        'favorable_contacts_electrostatics': [favorable_contact_list],
        'local_favorable_contacts_electrostatics': [local_favorable_contacts_list],
        'non_local_favorable_contacts_electrostatics': [non_local_favorable_contacts_list],
        # Add other metrics as needed
    })

    return burial_sasa_df


def main():
    parser = argparse.ArgumentParser(description="Process PDB file and output metrics as JSON.")
    parser.add_argument("--input", required=True, help="Path to the input PDB file.")
    parser.add_argument("--output", required=True,  help="Path to the output JSON file.")
    args = parser.parse_args()

    # Process the PDB file
    df = pdb2sasa_burial(args.input)

    # Optionally, add a column for the file name if needed
    df['name'] = [args.input.split("/")[-1]]

    # Save the DataFrame to a JSON file
    df.to_json(args.output, orient='records')

if __name__ == "__main__":
    main()
