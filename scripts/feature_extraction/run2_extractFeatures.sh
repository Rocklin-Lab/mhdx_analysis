#!/bin/bash

# Measure the entire script execution time
start_time=$(date +%s)

eval "$(conda shell.bash hook)"

scripts_path="/projects/b1107/allan/HDX_analysis/HX_datasets/features/scripts"

# Define paths to executables
adopt_exe="python /projects/b1107/allan/software/ADOPT/run-adopt.py"
fldpnn_exe="python /projects/b1107/allan/software/fldpnn/run_flDPnn.py"
freesasa_exe="python /projects/b1107/allan/software/freesasa-python/run_freesasa.py"
colabfold_batch_exe="/projects/b1107/allan/software/localcolabfold/colabfold-conda/bin/colabfold_batch"
relax_exe="/projects/p30802/Rosetta/main/source/bin/relax.linuxgccrelease"
pdb2rosetta_local_exe="python $scripts_path/pdb2rosetta_local.py"
plddt_exe="python $scripts_path/plddt.py"
score_designs_AF_exe="python /projects/p30802/grocklin/score_monomeric_designs/scripts/score_designs_AF.py"
pdb2dssp_exe="python $scripts_path/pdb2dssp.py"
pdb2hbond_exe="python $scripts_path/pdb2hbond.py"
sequence2ifeature_global_exe="python $scripts_path/sequence2ifeature_global.py"
sequence2ifeature_local_exe="python $scripts_path/sequence2ifeature_local.py"
pdb2contacts_exe="python $scripts_path/pdb2contacts.py"
pdb2rg_exe="python $scripts_path/pdb2rg.py"

# Check if at least two arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <path-to-fasta-file> <topology>"
    exit 1
fi

# Define other variables
fasta=$1
topology=$2
name=$(basename "$fasta" | sed "s|.fasta||g" | sed "s|.a3m||g")
error_file="${name}.err"

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

rosetta_relax="$(readlink -f rosetta_relax/pdbs/${name}.pdb)"

# EXTRACTING FEATURES

# Rosetta local features
echo "Extracting Rosetta local features..."
output_file="features/rosetta/local/${name}.json"
if [ ! -f "$output_file" ]; then
    $pdb2rosetta_local_exe --input "rosetta_relax/pdbs/${name}.pdb" --output "$output_file"
    check_output "$output_file"
else
    echo "Rosetta local features for ${name} already exist."
fi

# PLDDT
echo "Extracting PLDDT/PAE features..."
plddt_output="features/plddt/${name}.json"
if [ ! -f "$plddt_output" ]; then
        $plddt_exe --fasta "$fasta" --plddt $(readlink -f "af_prediction/output/${name}_scores_rank_001_alphafold2_ptm_model*json") --output "$plddt_output"
    check_output "$plddt_output"
else
    echo "PLDDT features for ${name} already exist."
fi


# Run in-house contacts script
echo "Running contacts code..."
contacts_output="features/contacts/${name}.json"
if [ ! -f "$contacts_output" ]; then
    $pdb2contacts_exe --input "rosetta_relax/pdbs/${name}.pdb" --output "$contacts_output"
    check_output "$contacts_output"
else
    echo "Contacts output for ${name} already exists."
fi

# Activate PyRosetta environment for Rosetta scoring
echo "Activating PyRosetta environment..."
conda activate pyrosetta
echo "Running Rosetta scoring..."
rosetta_output="features/rosetta/haddox/${name}.pdb.scores.csv"
if [ ! -f "$rosetta_output" ]; then
    $score_designs_AF_exe --pdb_path "rosetta_relax/pdbs/${name}.pdb" --output_dir "features/rosetta/haddox"
    check_output "$rosetta_output"
else
    echo "Rosetta scoring output for ${name} already exists."
fi

echo "Running HBONDS features..."
hbonds_output="features/hbonds/${name}.json"
if [ ! -f "$hbonds_output" ]; then
#    echo $pdb2hbond_exe --input "$rosetta_relax" --output "$hbonds_output"
    $pdb2hbond_exe --input "$rosetta_relax" --output "$hbonds_output"
    check_output "$hbonds_output"
else
    echo "HBONDS output for ${name} already exists."
fi

# Run DSSP/ABEGO
echo "Running DSSP/ABEGO..."
dssp_output="features/abego_dssp/${name}.json"
echo $dssp_output
if [ ! -f "$dssp_output" ]; then
    $pdb2dssp_exe --input "rosetta_relax/pdbs/${name}.pdb" --output "$dssp_output" --topology "$topology"
    check_output "$dssp_output"
else
    echo "DSSP/ABEGO output for ${name} already exists."
fi

# Run RG - radius of gyration
echo "Running Radius of gyration..."
rg_output="features/rg/${name}.json"
if [ ! -f "$rg_output" ]; then
    $pdb2rg_exe --input "rosetta_relax/pdbs/${name}.pdb" --output "$rg_output"
    check_output "$rg_output"
else
    echo "RG output for ${name} already exists."
fi

# Deactivate PyRosetta environment
conda deactivate

# Run iFeature (Global and Local)
echo "Activating iFeature environment..."
conda activate ifeature
echo "Running iFeature for global features..."
ifeature_global_output="features/ifeature/${name}_ifeature_global.json"
if [ ! -f "$ifeature_global_output" ]; then
    $sequence2ifeature_global_exe -i "$fasta" -o "$ifeature_global_output"
    check_output "$ifeature_global_output"
else
    echo "iFeature global features for ${name} already exist."
fi

echo "Running iFeature for local features..."
ifeature_local_output="features/ifeature/${name}_ifeature_local.json"
echo $dssp_output
if [ ! -f "$ifeature_local_output" ]; then
    $sequence2ifeature_local_exe --fasta "$fasta" --dssp "$dssp_output" --output "$ifeature_local_output"
    check_output "$ifeature_local_output"
else
    echo "iFeature local features for ${name} already exist."
fi
conda deactivate

# Run freeSASA
echo "Running freeSASA..."
freesasa_output="features/sasa/${name}.json"
if [ ! -f "$freesasa_output" ]; then
    $freesasa_exe --input "rosetta_relax/pdbs/${name}.pdb" --output "$freesasa_output"
    check_output "$freesasa_output"
else
    echo "freeSASA output for ${name} already exists."
fi


# Split into numeric and list types
echo "Running concatenation code..."
numeric_output="features/concatenated/${name}_numeric.json"
lists_output="features/concatenated/${name}_lists.json"
echo python $scripts_path/generate_numeric_list_dataframes.py $(ls features/*/${name}.json features/rosetta/haddox/${name}.pdb.scores.csv features/rosetta/local/${name}.json features/ifeature/${name}_ifeature_global.json  features/ifeature/${name}_ifeature_local.json) --output_numeric $numeric_output --output_lists $lists_output
if [ ! -f "$numeric_output" ] || [ ! -f "$lists_output" ]; then
   python $scripts_path/generate_numeric_list_dataframes.py $(ls features/*/${name}.json features/rosetta/haddox/${name}.pdb.scores.csv features/rosetta/local/${name}.json features/ifeature/${name}_ifeature_global.json  features/ifeature/${name}_ifeature_local.json) --output_numeric $numeric_output --output_lists $lists_output
   #python generate_numeric_list_dataframes.py $(ls features/*/${name}.json features/rosetta/haddox/${name}.pdb.scores.csv features/rosetta/local/${name}.json) --output_numeric $numeric_output --output_lists $lists_output
else
    echo "Concatenated outputs for ${name} already exist."
fi

echo "Generating final dataframe..."
final_output="features/final/${name}.json"
if [ ! -f "$final_output" ]; then
	python $scripts_path/generate_final_dataframe.py --numeric_path $numeric_output --list_path $lists_output --topology $topology --output $final_output
else
	echo "Final output for ${name} already exist."
fi

echo "Workflow completed."

# Measure end time and calculate total execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))

echo "FeatureExtraction Total execution time: ${execution_time} seconds" >> LOG
