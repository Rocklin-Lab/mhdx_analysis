import argparse
import pandas as pd
from io import StringIO

def extract_sequence_from_pdb(pdb_file_path):
    sequence = ""
    seen_residues = set()

    # Read the lines from the PDB file
    with open(pdb_file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract information for CA atoms
                atom_name = line[12:16].strip()
                if atom_name == "CA":
                    # Extract residue information
                    residue_name = line[17:20].strip()
                    residue_number = int(line[22:26].strip())
                    
                    # Create a unique identifier for the residue position
                    residue_position = f"{residue_name}_{residue_number}"
                    
                    # Check if we have already seen this residue
                    if residue_position not in seen_residues:
                        seen_residues.add(residue_position)
                        
                        # Map three-letter amino acid codes to single characters
                        aa_mapping = {
                            'HIS': 'H', 'MET': 'M', 'SER': 'S', 'GLU': 'E', 'ALA': 'A',
                            'LYS': 'K', 'ARG': 'R', 'VAL': 'V', 'LEU': 'L', 'ASN': 'N',
                            'GLY': 'G', 'ILE': 'I', 'PRO': 'P', 'ASP': 'D', 'GLN': 'Q',
                            'PHE': 'F', 'TRP': 'W', 'TYR':'Y', 'THR':'T', 'CYS':'C'
                        }
                        
                        # If the residue is in the mapping, use the mapped single character
                        # Otherwise, use the original three-letter code
                        sequence += aa_mapping.get(residue_name, residue_name)

    return sequence


def process_pdb(pdb_file_path):

    # Read the lines from the file
    with open(pdb_file_path, 'r') as file:
        lines = file.readlines()

    # Identify the line number where the table starts
    start_line = None
    for i, line in enumerate(lines):
        if line.startswith("label"):
            start_line = i
            break

    columns = lines[i].split()
#     print(columns)
    # Check if the start_line is found
    if start_line is not None:
        # Extract the table lines
        table_lines = lines[start_line + 1:]

        # Create a DataFrame from the table lines
        df = pd.read_csv(
            StringIO('\n'.join(table_lines)),
            delim_whitespace=True,
            comment='#', names=columns
        )

        # Map three-letter amino acid codes to single characters
        aa_mapping = {
                        'HIS': 'H', 'MET': 'M', 'SER': 'S', 'GLU': 'E', 'ALA': 'A',
                        'LYS': 'K', 'ARG': 'R', 'VAL': 'V', 'LEU': 'L', 'ASN': 'N',
                        'GLY': 'G', 'ILE': 'I', 'PRO': 'P', 'ASP': 'D', 'GLN': 'Q',
                        'PHE': 'F', 'TRP': 'W', 'TYR':'Y', 'THR':'T', 'CYS':'C' 
                    }
                    

        # Extract amino acid label and number
        df['aa'] = df['label'].apply(lambda x: x.split('_')[0].split(':')[0])
        df['aa'] = df['aa'].apply(lambda x: aa_mapping.get(x, x))

        # Extract amino acid number
        df['Amino Acid Number'] = df['label'].apply(lambda x: int(x.split('_')[-1]) if x.split('_')[-1].isdigit() else None)

        return df

    else:
        print("Table not found in the PDB file.")
        
        
def get_list_of_values(pdb_file_path):
    
    name = pdb_file_path.split("/")[-1].strip(".pdb")
    
    df = process_pdb(pdb_file_path)
    
    sequence = ''.join(df['aa'].values[2:-1])
    
    values = df.values[2:-1, 1:-2].T\
    
    print(name)
    
    return [[name, sequence] + list(values)]
    

def process_files_and_output_json(inputs, output):

    l = []
    for f in inputs:
        l.append(get_list_of_values(f))

    columns = ["name", "sequence"] + [
        'fa_atr', 'fa_rep', 'fa_sol', 'fa_in tra_atr_xover4', 'fa_intra_rep_xover4', 'fa_intra_sol_xover4', 
        'lk_ball', 'lk_ball_iso', 'lk_ball_bridge', 'lk_ball_bridge_uncpl', 'fa_elec',
        'fa_intra_elec', 'pro_close','hbond_sr_bb','hbond_lr_bb', 'hbond_bb_sc', 'hbond_sc', 'dslf_fa13',
        'coordinate_constraint','omega', 'fa_dun_dev', 'fa_dun_rot', 'fa_dun_semi', 'p_aa_pp', 'hxl_tors',
        'ref', 'rama_prepro', 'gen_bonded', 'total']
    df = pd.DataFrame([i[0] for i in l], columns=columns)
    df.to_json(output)

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Process PDB files and output to a JSON file.")
    parser.add_argument('--input', type=str, nargs='+', help='Path to input PDB files')
    parser.add_argument('--output', type=str, help='Path to output JSON file')

    # Parse arguments
    args = parser.parse_args()

    # Process files and output JSON
    process_files_and_output_json(args.input, args.output)

if __name__ == "__main__":
    main()

