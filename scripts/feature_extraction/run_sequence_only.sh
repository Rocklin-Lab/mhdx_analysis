
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

echo "Processing ${name} ..."

# Setup directories
echo "Setting up directories..."
mkdir -p af_prediction/{input,output,pdbs} \
         rosetta_relax/pdbs \
         features/{ifeature,rosetta/{,haddox,local},sasa,adopt,contacts,fldpnn,loops,plddt,abego_dssp,hbonds,concatenated,rg,final}

# Copy fasta to input directory
if [ ! -f "af_prediction/input/${name}.a3m" ]; then
    cp "$fasta" "af_prediction/input/${name}.a3m"
else
    echo "${name}.a3m already exists, proceeding to next step."
fi

# FUNCTION TO CHECK OUTPUT AND LOG ERROR
check_output() {
    if [ ! -f "$1" ]; then
        echo "Error: Expected output file $1 not found." >> "$error_file"
        echo "Check $error_file for details."
        exit 1
    fi
}

## STRUCTURE PREDICTION
#
## AlphaFold prediction
#echo "Running AlphaFold prediction..."
#if [ ! -f "af_prediction/pdbs/${name}.pdb" ]; then
#    $colabfold_batch_exe af_prediction/input af_prediction/output --amber --num-relax 1
#    af_prediction=$(readlink -f "af_prediction/output/${name}_relaxed_rank_001_alphafold2_ptm_model"*)
#    mv "$af_prediction" "af_prediction/pdbs/${name}.pdb"
#    check_output "af_prediction/pdbs/${name}.pdb"
#else
#    echo "AlphaFold prediction for ${name} already exists."
#fi
#af_prediction="$(readlink -f af_prediction/pdbs/${name}.pdb)"
#
## Rosetta Relax
#echo "Running Rosetta Relax..."
#if [ ! -f "rosetta_relax/pdbs/${name}.pdb" ]; then
#    $relax_exe -s "$(readlink -f af_prediction/pdbs/${name}.pdb)" -constrain_relax_to_start_coords -relax:ramp_constraints false --beta
#    rosetta_relax=$(readlink -f "${name}"*_0001.pdb)
#    mv "$rosetta_relax" "rosetta_relax/pdbs/${name}.pdb"
#    check_output "rosetta_relax/pdbs/${name}.pdb"
#else
#    echo "Rosetta Relax output for ${name} already exists."
#fi
#rosetta_relax="$(readlink -f rosetta_relax/pdbs/${name}.pdb)"

# EXTRACTING FEATURES

## Rosetta local features
#echo "Extracting Rosetta local features..."
#output_file="features/rosetta/local/${name}.json"
#if [ ! -f "$output_file" ]; then
#    $pdb2rosetta_local_exe --input "rosetta_relax/pdbs/${name}.pdb" --output "$output_file"
#    check_output "$output_file"
#else
#    echo "Rosetta local features for ${name} already exist."
#fi
#
## PLDDT
#echo "Extracting PLDDT/PAE features..."
#plddt_output="features/plddt/${name}.json"
#if [ ! -f "$plddt_output" ]; then
#        $plddt_exe --fasta "$fasta" --plddt $(readlink -f "af_prediction/output/${name}_scores_rank_001_alphafold2_ptm_model*json") --output "$plddt_output"
#    check_output "$plddt_output"
#else
#    echo "PLDDT features for ${name} already exist."
#fi
#
#
## Run in-house contacts script
#echo "Running contacts code..."
#contacts_output="features/contacts/${name}.json"
#if [ ! -f "$contacts_output" ]; then
#    $pdb2contacts_exe --input "rosetta_relax/pdbs/${name}.pdb" --output "$contacts_output"
#    check_output "$contacts_output"
#else
#    echo "Contacts output for ${name} already exists."
#fi
#
## Activate PyRosetta environment for Rosetta scoring
#echo "Activating PyRosetta environment..."
#conda activate pyrosetta
#echo "Running Rosetta scoring..."
#rosetta_output="features/rosetta/haddox/${name}.pdb.scores.csv"
#if [ ! -f "$rosetta_output" ]; then
#    $score_designs_AF_exe --pdb_path "rosetta_relax/pdbs/${name}.pdb" --output_dir "features/rosetta/haddox"
#    check_output "$rosetta_output"
#else
#    echo "Rosetta scoring output for ${name} already exists."
#fi
#
#echo "Running HBONDS features..."
#hbonds_output="features/hbonds/${name}.json"
#if [ ! -f "$hbonds_output" ]; then
##    echo $pdb2hbond_exe --input "$rosetta_relax" --output "$hbonds_output"
#    $pdb2hbond_exe --input "$rosetta_relax" --output "$hbonds_output"
#    check_output "$hbonds_output"
#else
#    echo "HBONDS output for ${name} already exists."
#fi
#
## Run DSSP/ABEGO
#echo "Running DSSP/ABEGO..."
#dssp_output="features/abego_dssp/${name}.json"
#if [ ! -f "$dssp_output" ]; then
#    $pdb2dssp_exe --input "rosetta_relax/pdbs/${name}.pdb" --output "$dssp_output" --topology "$topology"
#    check_output "$dssp_output"
#else
#    echo "DSSP/ABEGO output for ${name} already exists."
#fi
#
## Run RG - radius of gyration
#echo "Running Radius of gyration..."
#rg_output="features/rg/${name}.json"
#if [ ! -f "$rg_output" ]; then
#    $pdb2rg_exe --input "rosetta_relax/pdbs/${name}.pdb" --output "$rg_output"
#    check_output "$rg_output"
#else
#    echo "RG output for ${name} already exists."
#fi

# Deactivate PyRosetta environment
conda deactivate

# Run ADOPT
echo "Activating ADOPT environment..."
conda activate adopt
echo "Running ADOPT..."
adopt_output="features/adopt/${name}.json"
if [ ! -f "$adopt_output" ]; then
    $adopt_exe --fasta "$fasta" --output "$adopt_output"
    check_output "$adopt_output"
else
    echo "ADOPT output for ${name} already exists."
fi
conda deactivate

# Run fldpnn
echo "Activating flDPnn environment..."
conda activate flDPnn
echo "Running fldpnn..."
fldpnn_output="features/fldpnn/${name}.csv"
if [ ! -f "$fldpnn_output" ]; then
    $fldpnn_exe "$fasta" "features/fldpnn/${name}.csv"
    check_output "$fldpnn_output"
else
    echo "fldpnn output for ${name} already exists."
fi
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
if [ ! -f "$ifeature_local_output" ]; then
    $sequence2ifeature_local_exe --fasta "$fasta" --dssp "$dssp_output" --output "$ifeature_local_output"
    check_output "$ifeature_local_output"
else
    echo "iFeature local features for ${name} already exist."
fi
conda deactivate

## Run freeSASA
#echo "Running freeSASA..."
#freesasa_output="features/sasa/${name}.json"
#if [ ! -f "$freesasa_output" ]; then
#    $freesasa_exe --input "rosetta_relax/pdbs/${name}.pdb" --output "$freesasa_output"
#    check_output "$freesasa_output"
#else
#    echo "freeSASA output for ${name} already exists."
#fi
#
#
## Split into numeric and list types
#echo "Running concatenation code..."
#numeric_output="features/concatenated/${name}_numeric.json"
#lists_output="features/concatenated/${name}_lists.json"
#if [ ! -f "$numeric_output" ] || [ ! -f "$lists_output" ]; then
#   python $scripts_path/generate_numeric_list_dataframes.py $(ls features/*/${name}.json features/rosetta/haddox/${name}.pdb.scores.csv features/rosetta/local/${name}.json features/ifeature/${name}_ifeature_global.json  features/ifeature/${name}_ifeature_local.json) --output_numeric $numeric_output --output_lists $lists_output
#   #python generate_numeric_list_dataframes.py $(ls features/*/${name}.json features/rosetta/haddox/${name}.pdb.scores.csv features/rosetta/local/${name}.json) --output_numeric $numeric_output --output_lists $lists_output
#else
#    echo "Concatenated outputs for ${name} already exist."
#fi
#
#echo "Generating final dataframe..."
#final_output="features/final/${name}.json"
#if [ ! -f "$final_output" ]; then
#	python $scripts_path/generate_final_dataframe.py --numeric_path $numeric_output --list_path $lists_output --topology $topology --output $final_output
#else
#	echo "Final output for ${name} already exist."
#fi
#
#echo "Workflow completed."

