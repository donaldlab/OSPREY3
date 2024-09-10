# this file renames each chain for all PDBs in a directory
# for DexDesign (when using analyze.py), it is recommended to use:
# target = z, peptide = y
# this avoids issues with MASTER returns
# NOTE: this also auto-removes PDB REMARKS and ANISOU records

from Bio.PDB import PDBParser, PDBIO
import os
import warnings

# hide atom name warnings
warnings.simplefilter('ignore')


def chain_renamer(complex_pdb: str, desired_target_id: str, desired_peptide_id: str):

    print("\nRenaming %s IDs. Target ID: %s, Peptide ID: %s.\n" % (complex_pdb, desired_target_id, desired_peptide_id))

    # get the old IDs
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("complex", complex_pdb)
    model = structure[0]

    chain_ids = []
    for chain in model:
        chain_ids.append(chain.id)
    old_target_chain = chain_ids[0]
    old_peptide_chain = chain_ids[1]

    # change both IDs
    for chain in model:
        if chain.id == old_target_chain:
            chain.id = desired_target_id
        elif chain.id == old_peptide_chain:
            chain.id = desired_peptide_id

    # write out the new PDB with suffix -renamed.pdb
    io = PDBIO()
    io.set_structure(structure)
    filename = complex_pdb[:-4] + "-renamed.pdb"
    io.save(filename)


# example call
# chain_renamer("8gaj-trim.pdb", "z", "y")
