#!/bin/bash

# Generates template coordinates from a PDB file's atom records.
# Can be used to generate template coordinates for, e.g., a small molecule ligand.
#
# ./gen-templ-coords.sh sch7.pdb 38Z
# (script) (structure) (residue ID in structure)

if [[ $# -ne 2 ]]; then
	echo "Expecting 2 parameters: <structure-file> and <residue ID>" >&2
	exit 1
fi

echo $2 $(wc --lines < $1)
tr --squeeze-repeats ' ' < $1 | cut --delimiter=' ' --fields=3,7,8,9
