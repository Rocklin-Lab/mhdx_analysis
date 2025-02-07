#!/bin/bash

# Measure the entire script execution time
start_time=$(date +%s)

eval "$(conda shell.bash hook)"

scripts_path="/projects/b1107/allan/software/scripts"

# Define paths to executables
adopt_exe="python /projects/b1107/allan/software/ADOPT/run-adopt.py"
fldpnn_exe="python /projects/b1107/allan/software/fldpnn/run_flDPnn.py"

# Check if at least two arguments are provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <path-to-fasta-folder>"
    exit 1
fi

# Define other variables
fastas_folder=$1

# FUNCTION TO CHECK OUTPUT AND LOG ERROR
check_output() {
    if [ ! -f "$1" ]; then
        echo "Error: Expected output file $1 not found." >> "$error_file"
        echo "Check $error_file for details."
        exit 1
    fi
}

# Setup directories
echo "Setting up directories..."
mkdir -p af_prediction/{input,output,pdbs} \
         rosetta_relax/pdbs \
         features/{ifeature,rosetta/{,haddox,local},sasa,adopt,contacts,fldpnn,loops,plddt,abego_dssp,hbonds,concatenated,rg,final}


# Run ADOPT
echo "Activating ADOPT environment..."
conda activate adopt
echo "Running ADOPT..."
for fasta in $(readlink -f ${fastas_folder}/*); do
	name=$(basename "$fasta" | sed "s|.fasta||g" | sed "s|.a3m||g")
	error_file="${name}.err"
	adopt_output="features/adopt/${name}.json"
	if [ ! -f "$adopt_output" ]; then
    		$adopt_exe --fasta "$fasta" --output "$adopt_output"
    		check_output "$adopt_output"
	else
    		echo "ADOPT output for ${name} already exists."
	fi
done
conda deactivate

# Run fldpnn
echo "Activating flDPnn environment..."
conda activate flDPnn
echo "Running fldpnn..."
for fasta in $(readlink -f ${fastas_folder}/*); do
	name=$(basename "$fasta" | sed "s|.fasta||g" | sed "s|.a3m||g" )
	error_file="${name}.err"
	fldpnn_output="features/fldpnn/${name}.csv"
	if [ ! -f "$fldpnn_output" ]; then
    		echo $fldpnn_exe "$fasta" "features/fldpnn/${name}.csv"
    		$fldpnn_exe "$fasta" "features/fldpnn/${name}.csv"
    		check_output "$fldpnn_output"
	else
    		echo "fldpnn output for ${name} already exists."
	fi
done
conda deactivate


echo "Workflow completed."
# Measure end time and calculate total execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

echo "ADOPT/FLPNN Total execution time: ${execution_time} seconds"
echo "ADOPT/FLPNN Total execution time: ${execution_time} seconds" >> LOGFILE
