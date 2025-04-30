#!/bin/bash

# Measure script execution time
start_time=$(date +%s)

eval "$(conda shell.bash hook)"

# Path to ColabFold executable
colabfold_batch_exe="path-to/software/localcolabfold/colabfold-conda/bin/colabfold_batch"

# Check for required input argument
if [ "$#" -lt 1 ]; then
    echo "Use an allocation with GPU for this!"
    echo "Usage: $0 <path-to-fasta-or-a3m-folder>"
    exit 1
fi

fasta_folder=$1

# Create required output folders
mkdir -p af_prediction/{output,pdbs}

# Function to check output existence
check_output() {
    if [ ! -f "$1" ]; then
        echo "Error: Expected output file $1 not found." >> "$error_file"
        echo "Check $error_file for details."
        exit 1
    fi
}

# Run ColabFold batch prediction
echo "Running AlphaFold prediction..."
$colabfold_batch_exe "$fasta_folder" af_prediction/output --amber --num-relax 1

# Collect and rename relaxed models
for fasta in "$fasta_folder"/*; do 
    name=$(basename "$fasta" | sed "s|.fasta||g" | sed "s|.a3m||g")
    error_file="${name}.err"
    final_pdb="af_prediction/pdbs/${name}.pdb"
    if [ ! -f "$final_pdb" ]; then
        af_model=$(readlink -f af_prediction/output/${name}_relaxed_rank_001_alphafold2_ptm_model*)
        cp "$af_model" "$final_pdb"
        check_output "$final_pdb"
    else
        echo "AlphaFold model for ${name} already exists."
    fi
done

# Report execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "AlphaFold structure prediction completed in ${execution_time} seconds."
