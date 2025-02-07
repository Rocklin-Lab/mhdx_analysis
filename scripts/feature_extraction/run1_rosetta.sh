#!/bin/bash

# Measure the entire script execution time
start_time=$(date +%s)


eval "$(conda shell.bash hook)"

scripts_path="/projects/b1107/allan/HDX_analysis/HX_datasets/features/scripts"

# Define paths to executables
relax_exe="/projects/p30802/Rosetta/main/source/bin/relax.linuxgccrelease"
pdb2hbond_exe="python $scripts_path/pdb2hbond.py"

# Check if at least two arguments are provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <path-to-af-models-folder>" 
    exit 1
fi

# Define other variables
pdbs_folder=$1

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

# Rosetta Relax
echo "Running Rosetta Relax..."
for pdb in $(ls ${pdbs_folder}); do 
	name=$(basename "$pdb" | sed "s|\.pdb$||g")
	error_file="${name}.err"
	if [ ! -f "rosetta_relax/pdbs/${name}.pdb" ]; then
		$relax_exe -s "$(readlink -f af_prediction/pdbs/${name}.pdb)" -constrain_relax_to_start_coords -relax:ramp_constraints false --beta
		rosetta_relax=$(readlink -f "${name}"*_0001.pdb)
		mv "$rosetta_relax" "rosetta_relax/pdbs/${name}.pdb"
		check_output "rosetta_relax/pdbs/${name}.pdb"
	else
	echo "Rosetta Relax output for ${name} already exists."
	fi
done

echo "Rosetta Relax finished"

echo "Activating pyrosetta environment"
conda activate pyrosetta
echo "Running HBONDS features..."
for pdb in $(readlink -f "rosetta_relax/pdbs/*pdb"); do
    name=$(basename "$pdb" | sed "s|\.pdb$||g")    
    error_file="${name}.err"
    hbonds_output="features/hbonds/${name}.json"
    if [ ! -f "$hbonds_output" ]; then
    #    echo $pdb2hbond_exe --input "$rosetta_relax" --output "$hbonds_output"
        $pdb2hbond_exe --input "$pdb" --output "$hbonds_output"
        check_output "$hbonds_output"
    else
        echo "HBONDS output for ${name} already exists."
    fi
done

# Measure end time and calculate total execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

echo "RosettaRelax Total execution time: ${execution_time} seconds" 
echo "RosettaRelax Total execution time: ${execution_time} seconds" >> LOGFILE
