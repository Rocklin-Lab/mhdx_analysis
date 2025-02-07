import json
import argparse
from pyrosetta import init, pose_from_file

# Define the order of amino acids for mapping
AMINO_ACID_ORDER = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 
                    'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def analyze_prolines(pdb_files):
    """
    Analyze cis/trans conformations of prolines in PDB structures.

    Parameters:
        pdb_files (list): List of PDB file paths.

    Returns:
        dict: A list of dictionaries containing the sequence, number of cis prolines,
              number of trans prolines, and cis/trans arrays for each structure.
    """
    results = []

    for pdb in pdb_files:
        try:
            # Load the structure
            pose = pose_from_file(pdb)
            sequence = pose.sequence()

            # Initialize cis and trans arrays
            cis_array = [0] * 20
            trans_array = [0] * 20

            cis_count = 0
            trans_count = 0

            # Loop over residues to analyze prolines
            for i in range(2, pose.size() + 1 ):  # Start from 2 since Pro-1 has no preceding residue
#                print(i, pose.omega(i))
                if pose.residue(i).name3() == "PRO":
                    # Identify the preceding residue
                    prev_res = pose.residue(i - 1).name1()  # Single-letter code
                    if prev_res in AMINO_ACID_ORDER:
                        index = AMINO_ACID_ORDER.index(prev_res)

                        # Check omega torsion angle to classify cis/trans
                        omega = pose.omega(i-1)
                        print(i, omega)
                        if abs(omega) < 30:  # cis conformation has omega ~ 0°
                            cis_array[index] += 1
                            cis_count += 1
                        else:  # trans conformation has omega ~ 180°
                            trans_array[index] += 1
                            trans_count += 1

            # Append results
            results.append({
                "pdb": pdb,
                "sequence": sequence,
                "cis_prolines": cis_count,
                "trans_prolines": trans_count,
                "cis_array": cis_array,
                "trans_array": trans_array
            })
        except Exception as e:
            print(f"Error processing {pdb}: {e}")
            results.append({
                "pdb": pdb,
                "error": str(e)
            })

    return results

def main():
    # Argument parser
    parser = argparse.ArgumentParser(description="Analyze cis/trans prolines in PDB structures.")
    parser.add_argument("pdb_files", nargs="+", help="List of PDB file paths to analyze.")
    parser.add_argument("-o", "--output", required=True, help="Output JSON file.")
    args = parser.parse_args()

    # Initialize PyRosetta
    init()

    # Print amino acid mapping
    print("Amino Acid Mapping:", AMINO_ACID_ORDER)

    # Analyze prolines
    results = analyze_prolines(args.pdb_files)

    # Write results to JSON
    with open(args.output, "w") as json_file:
        json.dump(results, json_file, indent=4)
    print(f"Results written to {args.output}")

if __name__ == "__main__":
    main()

