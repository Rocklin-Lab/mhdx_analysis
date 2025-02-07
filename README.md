# mhdx_analysis

This repository contains scripts and notebooks to analyze results from **mhdx_pipeline** and **hdxrate_pipeline** associated with the *Ferrari 2025* paper.

## Table of Contents
1. [Compute Cooperativities](#compute-cooperativities)
2. [AF2 Prediction](#af2-prediction)
3. [Rosetta Relaxation](#rosetta-relaxation)
4. [Hydrogen Bond Extraction](#hydrogen-bond-extraction)
5. [Dependencies](#dependencies)
6. [How to Cite](#how-to-cite)

---

## Compute Cooperativities
The code to compute cooperativities (*normalized cooperativity* and *family-normalized cooperativity*) from scratch can be found in:

**Path:** `notebooks/DeriveCooperativities.ipynb`

### Requirements:
- **Processed Data:** Output from `hdxrate_pipeline` (`{hdxrate_output}/consolidated_results/deduplicated.json`).
- **Structural Data:** Dataframe with protein names and the number of expected protected amides based on their 3D structure.

### Structural Modeling Details:
- Structures from *Ferrari 2025* were generated using **AlphaFold2** via **ColabFold** (version X), followed by **RosettaRelax** applied to the top-scoring model.
- Hydrogen bond information was extracted using a custom **PyRosetta** script.

*Scripts for all these steps are provided. Ensure ColabFold and Rosetta are installed. Scripts include the exact parameters used in our study.*

---

## AF2 Prediction
To perform structure prediction using AlphaFold2:

```bash
bash mhdx_analysis/scripts/feature_extraction/run0_structure_prediction.sh {folder-to-fasta-or-a3m-files}
```

---

## Rosetta Relaxation
To relax the AlphaFold2 predicted structures using Rosetta:

```bash
bash mhdx_analysis/scripts/feature_extraction/run1_rosetta.sh {af_prediction/pdbs/}
```

---

## Hydrogen Bond Extraction
Extract hydrogen bond information from the relaxed structures.

**Run for each structure:**

```bash
python mhdx_analysis/scripts/feature_extraction/pdb2hbond.py --input "$rosetta_relax" --output "$hbonds_output"
```

After extraction, concatenate all H-bond JSON files and populate a column `PF` with the corresponding protein family.

---

## Dependencies
Ensure the following tools and libraries are installed:
- **ColabFold** (version 1.5.2): https://github.com/sokrypton/ColabFold
- **Rosetta** (version 2019.11) : https://rosettacommons.org/software/
- **PyRosetta** (version 4 2022.33) : https://www.pyrosetta.org/

---

## How to Cite
If you use this repository, please cite the following paper:

*Ferrari, A. et al. (2025). Title of the Paper. *Journal Name*. DOI: XXXXXX*

---

For any questions or issues, please contact `ajrferrari@gmail.com`.

