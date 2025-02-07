#!/bin/bash

# Check if seqkit is available
if ! command -v seqkit &> /dev/null; then
    echo "Error: 'seqkit' is not installed or not in PATH. Please install seqkit and try again." >&2
    exit 1
fi

# Base directory where the fasta files are located
BASE_DIR="."  # Update this path
OUTPUT_DIR="mmseqs_outputs_dynamic"
TMP_DIR="tmp_mmseqs"

# Ensure the output and temporary directories exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$TMP_DIR"

# List of FASTA files
FASTA_FILES=("Cold-Shock.fasta" "EEHEE.fasta" "EHEE.fasta" "HEEH.fasta" "HHH.fasta"
             "LysM.fasta" "PASTA.fasta" "PDB.fasta" "SH3.fasta" "WW.fasta")

#FASTA_FILES=("HHH.fasta")

# Iterate over each fasta file
for FASTA_FILE in "${FASTA_FILES[@]}"
do
    # Extract the base name (without extension)
    BASE_NAME=$(basename "$FASTA_FILE" .fasta)

    # Create a folder for each fasta file
    FASTA_FOLDER="$OUTPUT_DIR/$BASE_NAME"
    mkdir -p "$FASTA_FOLDER"

    # Define paths
    FASTA_PATH="$BASE_DIR/$FASTA_FILE"
    CLUSTER_OUTPUT_DIR="$FASTA_FOLDER"
    PAIRWISE_OUTPUT_DIR="$FASTA_FOLDER/pairwise_search"
    mkdir -p "$CLUSTER_OUTPUT_DIR"  # Create cluster directory for each file
    mkdir -p "$PAIRWISE_OUTPUT_DIR"  # Create pairwise directory

    echo "Processing $FASTA_FILE ..." >> outlog

    # Count the total number of sequences in the fasta file
    TOTAL_SEQUENCES=$(grep -c "^>" "$FASTA_PATH")
    if [ "${TOTAL_SEQUENCES}" -eq 0 ]; then
        echo "Error: No sequences found in $FASTA_FILE. Skipping."
        continue
    fi

    echo "Total sequences in $FASTA_FILE: ${TOTAL_SEQUENCES}" >> outlog

    # Dynamic clustering thresholds
    MIN_SEQ_IDS=(0.1 0.2 0.3 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75)

    for MIN_SEQ_ID in "${MIN_SEQ_IDS[@]}"
    do
        # Step 1: Cluster the sequences using easy-cluster with the current min-seq-id
        echo "Clustering sequences with min-seq-id ${MIN_SEQ_ID}..." >> outlog
        mmseqs easy-cluster "${FASTA_PATH}" "${CLUSTER_OUTPUT_DIR}/${BASE_NAME}_${MIN_SEQ_ID}" "${TMP_DIR}" --min-seq-id "${MIN_SEQ_ID}" -e 0.001

        # Check if cluster output exists
        CLUSTER_TSV="$CLUSTER_OUTPUT_DIR/${BASE_NAME}_${MIN_SEQ_ID}_cluster.tsv"
        if [ ! -f "$CLUSTER_TSV" ]; then
            echo "Error: Cluster file ${CLUSTER_TSV} not found. Skipping ${MIN_SEQ_ID}." >> outlog
            continue
        fi

        # Analyze the cluster sizes
        echo "Analyzing cluster sizes for min-seq-id $MIN_SEQ_ID..." >> outlog
        awk '{print $1}' "$CLUSTER_TSV" | sort | uniq -c | sort -nr > "$CLUSTER_OUTPUT_DIR/${BASE_NAME}_${MIN_SEQ_ID}_cluster_sizes.txt"

        # Get the largest cluster size
        LARGEST_CLUSTER_SIZE=$(head -n 1 "$CLUSTER_OUTPUT_DIR/${BASE_NAME}_${MIN_SEQ_ID}_cluster_sizes.txt" | awk '{print $1}')

        echo "Largest cluster size: $LARGEST_CLUSTER_SIZE" >> outlog

        # Check stop criteria: if the largest cluster is less than 10% of total sequences, stop
        if [ "$LARGEST_CLUSTER_SIZE" -le $((TOTAL_SEQUENCES / 10)) ]; then
            echo "Stop criteria met for min-seq-id $MIN_SEQ_ID: Largest cluster is <= 10% of total sequences." >> outlog
            cp ${CLUSTER_TSV} $CLUSTER_OUTPUT_DIR/${BASE_NAME}_final.tsv
            break
        fi
    done

    echo "Finished processing $FASTA_FILE with selected min-seq-id $MIN_SEQ_ID." >> outlog

    # Step 2: Group members by representative from the cluster.tsv
    echo "Grouping members by representatives..." >> outlog
    sort -k1,1 "$CLUSTER_OUTPUT_DIR/${BASE_NAME}_final.tsv" | uniq | awk '{print $1}' | uniq > "$CLUSTER_OUTPUT_DIR/${BASE_NAME}_representatives.txt"  # Get all unique representatives

    # Temporary file to store concatenated results
    FINAL_OUTPUT_TSV="$FASTA_FOLDER/${BASE_NAME}_pairwise.tsv"
    > "$FINAL_OUTPUT_TSV"  # Empty the file if it exists

    # Add header to the output file

    # Iterate through each representative in the cluster file
    while read -r rep; do
        # Extract all members for this representative
        awk -v rep="$rep" '$1 == rep {print $2}' $CLUSTER_OUTPUT_DIR/${BASE_NAME}_final.tsv  > "$CLUSTER_OUTPUT_DIR/${rep}_members.txt"
        #awk -v rep="$rep" '$1 == rep {print $2}' "$CLUSTER_OUTPUT_DIR/${rep}_members.txt"

        # Get the expected member count from members.txt
        expected_member_count=$(wc -l < "$CLUSTER_OUTPUT_DIR/${rep}_members.txt")

        # Create separate FASTA files for the representative and the members
        REP_FASTA="$PAIRWISE_OUTPUT_DIR/${rep}_rep.fasta"
        MEMBERS_FASTA="$PAIRWISE_OUTPUT_DIR/${rep}_members.fasta"

        # Add the representative sequence to the rep FASTA
        seqkit grep -p "$rep" "$FASTA_PATH" > "$REP_FASTA"

        # Add all members to the members FASTA
        seqkit grep -f "$CLUSTER_OUTPUT_DIR/${rep}_members.txt" "$FASTA_PATH" > "$MEMBERS_FASTA"

        # Create MMseqs2 databases for the representative and members
        REP_DB="$PAIRWISE_OUTPUT_DIR/${rep}_rep_db"
        MEMBER_DB="$PAIRWISE_OUTPUT_DIR/${rep}_members_db"

        mmseqs createdb "$REP_FASTA" "$REP_DB"
        mmseqs createdb "$MEMBERS_FASTA" "$MEMBER_DB"

        # Export the member database to FASTA using convert2fasta and count the number of sequences
        mmseqs convert2fasta "$MEMBER_DB" "$TMP_DIR/member_check.fasta"
        member_count=$(grep -c "^>" "$TMP_DIR/member_check.fasta")  # Count the number of sequences in the fasta

        echo "Running pairwise search for $rep (1 rep, $member_count members)..." >> outlog

        # Perform pairwise identity search between the representative and its members
        CLUSTER_OUTPUT="$PAIRWISE_OUTPUT_DIR/${rep}_pairwise"
        mmseqs search "$REP_DB" "$MEMBER_DB" "$CLUSTER_OUTPUT" "$TMP_DIR" --min-seq-id 0.0 --cov-mode 0 --min-aln-len 0 --prefilter-mode 2 --alignment-mode 3 --remove-tmp-files -e 1000

        # Convert the search results to a readable format (TSV) and append to the final output
        mmseqs convertalis "$REP_DB" "$MEMBER_DB" "$CLUSTER_OUTPUT" "$PAIRWISE_OUTPUT_DIR/${rep}_pairwise.tsv"

        # Append the results to the final output TSV with required formatting (rep, member, identity, coverage, e-value, bit score)
        awk -v rep="$rep" 'BEGIN{FS="\t"; OFS="\t"} {print rep, $2, $3, $4, $11, $12}' "$PAIRWISE_OUTPUT_DIR/${rep}_pairwise.tsv" >> "$FINAL_OUTPUT_TSV"

    done < "$CLUSTER_OUTPUT_DIR/${BASE_NAME}_representatives.txt"

    echo "Finished processing $FASTA_FILE. Final results saved to $FINAL_OUTPUT_TSV" >> outlog
done

echo "All tasks completed." >> outlog

