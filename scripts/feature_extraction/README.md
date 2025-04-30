# Protein Feature Extraction Pipeline

This repository contains a modular pipeline for extracting structural and sequence-based features from protein sequences. The pipeline integrates tools such as AlphaFold (via ColabFold), Rosetta, ADOPT, flDPnn, iFeature, and others. It is designed for high-throughput feature computation and supports downstream applications like protein stability and cooperativity prediction.

---

## ⚙️ Requirements

Each tool must be installed in its own conda environment to avoid dependency conflicts. The table below outlines the required tools, suggested environment names, and official installation links:

| Tool | Environment Name | Link |
|------|------------------|------|
| ADOPT | `adopt` | [ADOPT](https://github.com/PeptoneLtd/ADOPT) |
| flDPnn | `fldpnn` | [flDPnn](https://gitlab.com/sina.ghadermarzi/fldpnn) |
| ColabFold (local) | `colabfold` | [ColabFold](https://github.com/sokrypton/ColabFold) |
| FreeSASA | `freesasa` | [FreeSASA](https://freesasa.github.io/) |
| Rosetta (compiled) | system install | [RosettaCommons](https://docs.rosettacommons.org/demos/latest/tutorials/install_build/install_build) |
| PyRosetta | `pyrosetta` | [PyRosetta](https://www.pyrosetta.org/) |
| score_monomeric_designs | `pyrosetta` | [Haddox/score_monomeric_designs](https://github.com/Haddox/score_monomeric_designs) |
| iFeature | `ifeature` | [iFeature](https://github.com/Superzchen/iFeature) |

> ⚠️ Ensure that each environment is correctly activated within the corresponding script. Paths in the scripts (e.g., to ColabFold, ADOPT, Rosetta, etc.) must be updated according to your local installation.

---


### 📥 Input Parameters

- **`<fastas-folder>`**: Path to a folder containing `.fasta` or `.a3m` input files.
- **`<fasta-file>`**: A single FASTA file used for structure and feature extraction.
- **`<topology>`**: A string representing the expected secondary structure using `H` (helix) and `E` (strand), e.g., `HHEEH`.

Each script automatically checks for existing outputs and skips already completed steps to save time and avoid redundancy.

## 🚀 Usage

The pipeline is split into four main scripts. Run each one sequentially:

### 1. Structure prediction
```
bash run0_structure_prediction.sh <fastas-folder>
```

### 2. Disorder prediction (ADOPT + flDPnn)
```
bash run0_serial_adopt_fldpnn.sh <fastas-folder>
```

### 3. Rosetta relaxation
```
bash run1_rosetta.sh af_prediction/pdbs
```

### 4. Feature extraction
```
bash run2_extractFeatures.sh <fasta-file> <topology>
```

## Output organization

``` 
af_prediction/
├── output/           # ColabFold raw output
├── pdbs/             # Final predicted structures

rosetta_relax/
└── pdbs/             # Relaxed models

features/
├── adopt/            # ADOPT disorder predictions
├── fldpnn/           # flDPnn disorder predictions
├── ifeature/         # Global/local iFeature outputs
├── sasa/             # FreeSASA calculations
├── contacts/         # Contact maps
├── plddt/            # AlphaFold PLDDT and PAE
├── abego_dssp/       # DSSP/ABEGO analysis
├── hbonds/           # PyRosetta hydrogen bond analysis
├── rosetta/
│   ├── haddox/       # Rosetta scoring
│   └── local/        # Residue-level Rosetta terms
├── rg/               # Radius of gyration
├── concatenated/     # Merged features (split numeric/list)
└── final/            # Final unified feature JSON
```

