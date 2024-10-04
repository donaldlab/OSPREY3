#!/bin/csh -f
#SBATCH --output=submit.out
#SBATCH --mem=750G
#SBATCH --cpus-per-task=48
#SBATCH -A grisman --partition=grisman
~/dlab/henry/projects/TLDR/efxn_tune/osprey3-3.3/bin/osprey3 kstar --complex-confspace ./complex.ccsx --target-confspace ./target.ccsx --design-confspace ./peptide.ccsx --ensemble-dir ./ensembles --max-simultaneous-mutations 10000
rm *confdb
echo "completed"
