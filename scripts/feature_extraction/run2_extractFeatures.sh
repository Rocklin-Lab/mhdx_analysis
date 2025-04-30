#!/bin/bash

# Measure script execution time
start_time=$(date +%s)

eval "$(conda shell.bash hook)"

scripts_path="path-to/mhdx_analysis/scripts"

# Define paths to executables
relax_exe="path-to/Rosetta/main/source/bin/relax.linuxgccrelease"
pdb2rosetta_local_exe="python $scripts_path/pdb2rosetta_local.py"
plddt_exe="python $scripts_path/plddt.py"
score_designs_AF_exe="python path-to/score_monomeric_designs/scripts/score_designs.py"
pdb2dssp_exe="python $scripts_path/pdb2dssp.py"
pdb2hbond_exe="python $scripts_path/pdb2hbond.py"
sequence2ifeature_global_exe="python $scripts_path/sequence2ifeature_global.py"
sequence2ifeature_local_exe="python $scripts_path/sequence2ifeature_local.py"
pdb2contacts_exe="python $scripts_path/pdb2contacts.py"
pdb2rg_exe="python $scripts_path/pdb2rg.py"
freesasa_exe="python path-to/freesasa-python/run_freesasa.py"

# Input arguments
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <path-to-fasta-file> <topology>"
    exit 1
fi

fasta=$1
topology=$2
name=$(basename "$fasta" | sed "s|.fasta||g" | sed "s|.a3m||g")
error_file="${name}.err"

# Output directories
mkdir -p \
    af_prediction/output af_prediction/pdbs \
    rosetta_relax/pdbs \
    features/{ifeature,rosetta/{haddox,local},sasa,contacts,plddt,abego_dssp,hbonds,concatenated,rg,final}

check_output() {
    if [ ! -f "$1" ]; then
        echo "Error: Expected output file $1 not found." >> "$error_file"
        echo "Check $error_file for details."
        exit 1
    fi
}

# INPUT
rosetta_relax_pdb="rosetta_relax/pdbs/${name}.pdb"

# ----------------------------------------
# Rosetta local features
output_file="features/rosetta/local/${name}.json"
if [ ! -f "$output_file" ]; then
    echo "Extracting Rosetta local features..."
    $pdb2rosetta_local_exe --input "$rosetta_relax_pdb" --output "$output_file"
    check_output "$output_file"
else
    echo "Rosetta local features for ${name} already exist."
fi

# PLDDT/PAE
plddt_output="features/plddt/${name}.json"
if [ ! -f "$plddt_output" ]; then
    echo "Extracting PLDDT/PAE features..."
    score_json=$(readlink -f af_prediction/output/${name}_scores_rank_001_alphafold2_ptm_model*json)
    $plddt_exe --fasta "$fasta" --plddt "$score_json" --output "$plddt_output"
    check_output "$plddt_output"
else
    echo "PLDDT features for ${name} already exist."
fi

# Contact map
contacts_output="features/contacts/${name}.json"
if [ ! -f "$contacts_output" ]; then
    echo "Extracting contact features..."
    $pdb2contacts_exe --input "$rosetta_relax_pdb" --output "$contacts_output"
    check_output "$contacts_output"
else
    echo "Contact features for ${name} already exist."
fi

# ----------------------------------------
# Rosetta scoring (PyRosetta)
conda activate pyrosetta

rosetta_output="features/rosetta/haddox/${name}.pdb.scores.csv"
if [ ! -f "$rosetta_output" ]; then
    echo "Running Rosetta scoring..."
    $score_designs_AF_exe --pdb_path "$rosetta_relax_pdb" --output_dir "features/rosetta/haddox"
    check_output "$rosetta_output"
else
    echo "Rosetta scores for ${name} already exist."
fi

hbonds_output="features/hbonds/${name}.json"
if [ ! -f "$hbonds_output" ]; then
    echo "Extracting HBOND features..."
    $pdb2hbond_exe --input "$rosetta_relax_pdb" --output "$hbonds_output"
    check_output "$hbonds_output"
else
    echo "HBOND features for ${name} already exist."
fi

dssp_output="features/abego_dssp/${name}.json"
if [ ! -f "$dssp_output" ]; then
    echo "Extracting DSSP/ABEGO..."
    $pdb2dssp_exe --input "$rosetta_relax_pdb" --output "$dssp_output" --topology "$topology"
    check_output "$dssp_output"
else
    echo "DSSP/ABEGO features for ${name} already exist."
fi

rg_output="features/rg/${name}.json"
if [ ! -f "$rg_output" ]; then
    echo "Extracting radius of gyration..."
    $pdb2rg_exe --input "$rosetta_relax_pdb" --output "$rg_output"
    check_output "$rg_output"
else
    echo "RG features for ${name} already exist."
fi

conda deactivate

# ----------------------------------------
# iFeature (global and local)
conda activate ifeature

ifeature_global_output="features/ifeature/${name}_ifeature_global.json"
if [ ! -f "$ifeature_global_output" ]; then
    echo "Running iFeature (global)..."
    $sequence2ifeature_global_exe -i "$fasta" -o "$ifeature_global_output"
    check_output "$ifeature_global_output"
else
    echo "iFeature global for ${name} already exist."
fi

ifeature_local_output="features/ifeature/${name}_ifeature_local.json"
if [ ! -f "$ifeature_local_output" ]; then
    echo "Running iFeature (local)..."
    $sequence2ifeature_local_exe --fasta "$fasta" --dssp "$dssp_output" --output "$ifeature_local_output"
    check_output "$ifeature_local_output"
else
    echo "iFeature local for ${name} already exist."
fi

conda deactivate

# ----------------------------------------
# FreeSASA
freesasa_output="features/sasa/${name}.json"
if [ ! -f "$freesasa_output" ]; then
    echo "Extracting SASA..."
    $freesasa_exe --input "$rosetta_relax_pdb" --output "$freesasa_output"
    check_output "$freesasa_output"
else
    echo "SASA features for ${name} already exist."
fi

# ----------------------------------------
# Concatenate features
numeric_output="features/concatenated/${name}_numeric.json"
lists_output="features/concatenated/${name}_lists.json"
if [ ! -f "$numeric_output" ] || [ ! -f "$lists_output" ]; then
    echo "Concatenating features..."
    python $scripts_path/generate_numeric_list_dataframes.py \
        features/*/"${name}.json" \
        features/rosetta/haddox/"${name}.pdb.scores.csv" \
        features/rosetta/local/"${name}.json" \
        features/ifeature/"${name}_ifeature_global.json" \
        features/ifeature/"${name}_ifeature_local.json" \
        --output_numeric "$numeric_output" \
        --output_lists "$lists_output"
    check_output "$numeric_output"
    check_output "$lists_output"
else
    echo "Concatenated features for ${name} already exist."
fi

# Final assembly
final_output="features/final/${name}.json"
if [ ! -f "$final_output" ]; then
    echo "Generating final feature dataframe..."
    python $scripts_path/generate_final_dataframe.py \
        --numeric_path "$numeric_output" \
        --list_path "$lists_output" \
        --topology "$topology" \
        --output "$final_output"
else
    echo "Final feature file for ${name} already exists."
fi

# Done
end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Feature extraction completed in ${execution_time} seconds."
