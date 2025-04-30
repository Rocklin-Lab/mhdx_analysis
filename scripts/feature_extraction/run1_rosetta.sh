#!/bin/bash

# Start timer
start_time=$(date +%s)

eval "$(conda shell.bash hook)"

# Paths
scripts_path="path-to/mhdx_analysis/scripts"
relax_exe="path-to/Rosetta/main/source/bin/relax.linuxgccrelease"
pdb2hbond_exe="python $scripts_path/pdb2hbond.py"

# Check arguments
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <path-to-af-models-folder>"
    exit 1
fi

pdbs_folder=$1

# Output directories
mkdir -p rosetta_relax/pdbs features/hbonds

# Utility function
check_output() {
    if [ ! -f "$1" ]; then
        echo "Error: Expected output file $1 not found." >> "$error_file"
        echo "Check $error_file for details."
        exit 1
    fi
}

# Run Rosetta Relax
echo "Running Rosetta Relax..."
for pdb in "$pdbs_folder"/*.pdb; do
    name=$(basename "$pdb" .pdb)
    error_file="${name}.err"
    out_pdb="rosetta_relax/pdbs/${name}.pdb"
    if [ ! -f "$out_pdb" ]; then
        $relax_exe -s "$pdb" -constrain_relax_to_start_coords -relax:ramp_constraints false --beta
        relaxed_model=$(readlink -f "${name}"*_0001.pdb)
        mv "$relaxed_model" "$out_pdb"
        check_output "$out_pdb"
    else
        echo "Rosetta Relax output for ${name} already exists."
    fi
done
echo "Rosetta Relax complete."

# Extract HBOND features
echo "Activating PyRosetta environment..."
conda activate pyrosetta

echo "Extracting HBONDS features..."
for pdb in rosetta_relax/pdbs/*.pdb; do
    name=$(basename "$pdb" .pdb)
    error_file="${name}.err"
    out_json="features/hbonds/${name}.json"
    if [ ! -f "$out_json" ]; then
        $pdb2hbond_exe --input "$pdb" --output "$out_json"
        check_output "$out_json"
    else
        echo "HBONDS output for ${name} already exists."
    fi
done
conda deactivate

# Report execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "RosettaRelax + HBONDS completed in ${execution_time} seconds."
