#!/bin/bash
#SBATCH --account p30802    
#SBATCH --partition short 
#SBATCH --nodes=1								
#SBATCH --ntasks-per-node=1						
#SBATCH --time=1:00:00							
#SBATCH --mem-per-cpu=12G						
#SBATCH --job-name=template						
#SBATCH --output=outlog_template
#SBATCH --error=errlog_template

bash /projects/b1107/allan/HDX_analysis/HX_datasets/features/scripts/run_updated.sh template.fasta topology
