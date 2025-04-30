#!/bin/bash

# Measure script execution time
start_time=$(date +%s)

eval "$(conda shell.bash hook)"

# Paths to executables
adopt_exe="python path-to/software/ADOPT/run-adopt.py"
fldpnn_exe="python path-to/projects/b1107/allan/software/fldpnn/run_flDPnn.py"

# Input validation
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <path-to-fasta-folder>"
    exit 1
fi

fastas_folder=$1

# Output directories
mkdir -p features/adopt features/fldpnn

# Utility: Check output and log error
check_output() {
    if [ ! -f "$1" ]; then
        echo "Error: Expected output file $1 not found." >> "$error_file"
        echo "Check $error_file for details."
        exit 1
    fi
}

# -------------------------
# Run ADOPT
# -------------------------
echo "Activating ADOPT environment..."
conda activate adopt
echo "Running ADOPT..."
for fasta in "$fastas_folder"/*; do
    name=$(basename "$fasta" | sed 's|.fasta||g' | sed 's|.a3m||g')
    error_file="${name}.err"
    output_file="features/adopt/${name}.json"
    if [ ! -f "$output_file" ]; then
        $adopt_exe --fasta "$fasta" --output "$output_file"
        check_output "$output_file"
    else
        echo "ADOPT output for ${name} already exists."
    fi
done
conda deactivate

# -------------------------
# Run flDPnn
# -------------------------
echo "Activating flDPnn environment..."
conda activate flDPnn
echo "Running flDPnn..."
for fasta in "$fastas_folder"/*; do
    name=$(basename "$fasta" | sed 's|.fasta||g' | sed 's|.a3m||g')
    error_file="${name}.err"
    output_file="features/fldpnn/${name}.csv"
    if [ ! -f "$output_file" ]; then
        $fldpnn_exe "$fasta" "$output_file"
        check_output "$output_file"
    else
        echo "flDPnn output for ${name} already exists."
    fi
done
conda deactivate

# Report execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Disorder prediction (ADOPT + flDPnn) completed in ${execution_time} seconds."
