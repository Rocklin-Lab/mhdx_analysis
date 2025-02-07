#!/bin/bash

# Measure the entire script execution time
start_time=$(date +%s)

eval "$(conda shell.bash hook)"

scripts_path="/projects/b1107/allan/HDX_analysis/HX_datasets/features/scripts"

# Define paths to executables
colabfold_batch_exe="/projects/b1107/allan/software/localcolabfold/colabfold-conda/bin/colabfold_batch"

# Check if at least two arguments are provided
if [ "$#" -lt 1 ]; then
    echo "Use an allocation with GPU for this!"  
    echo "Usage: $0 <path-to-fasta-or-a3m-folder>"
    exit 1
fi

# Define other variables
fasta_folder=$1

# Setup directories
echo "Setting up directories..."
mkdir -p af_prediction/{input,output,pdbs} \
         rosetta_relax/pdbs \
         features/{ifeature,rosetta/{,haddox,local},sasa,adopt,contacts,fldpnn,loops,plddt,abego_dssp,hbonds,concatenated,rg,final}

# FUNCTION TO CHECK OUTPUT AND LOG ERROR
check_output() {
    if [ ! -f "$1" ]; then
        echo "Error: Expected output file $1 not found." >> "$error_file"
        echo "Check $error_file for details."
        exit 1
    fi
}

# STRUCTURE PREDICTION

# AlphaFold prediction
echo "Running AlphaFold prediction..."
$colabfold_batch_exe ${fasta_folder} af_prediction/output --amber --num-relax 1
for fasta in $(ls ${fasta_folder}); do 
	name=$(basename "$fasta" | sed "s|.fasta||g" | sed "s|.a3m||g")
        error_file="${name}.err"
	if [ ! -f "af_prediction/pdbs/${name}.pdb" ]; then
		af_prediction=$(readlink -f "af_prediction/output/${name}_relaxed_rank_001_alphafold2_ptm_model"*)
    		cp "$af_prediction" "af_prediction/pdbs/${name}.pdb"
    		check_output "af_prediction/pdbs/${name}.pdb"
	else
		echo "Final AlphaFold prediction model for ${name} already exists."
	fi
done

echo "Structure prediction finished!"

# Measure end time and calculate total execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

echo "Structure Prediction Total execution time: ${execution_time} seconds"
echo "Structure Prediction Total execution time: ${execution_time} seconds" >> LOGFILE
