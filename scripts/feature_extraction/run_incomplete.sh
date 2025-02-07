for i in $(ls *fasta); do name=$(echo $i | sed 's|.fasta||g'); if [[ ! -f features/final/${name}.json ]]; then echo $name; sbatch $(ls extract_features_${name}.sh); sleep 10; fi; done
