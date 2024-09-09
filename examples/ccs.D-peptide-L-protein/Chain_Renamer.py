# this file renames each chain in a PDB file
# for DexDesign (when using analyze.py), it is recommended to use:
# target = z, peptide = y
# this avoids issues with MASTER returns
# NOTE: this also auto-removes PDB REMARKS and ANISOU records

from Bio.PDB import PDBParser, PDBIO

def chain_renamer(pdb_filepath: str, desired_target_chain: str, desired_peptide_chain: str):

    # get the old IDs
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("complex", pdb_filepath)
    model = structure[0]

    chain_ids = []
    for chain in model:
        chain_ids.append(chain.id)
    old_target_chain = chain_ids[0]
    old_peptide_chain = chain_ids[1]

    # change both IDs
    for chain in model:
        if chain.id == old_target_chain:
            chain.id = desired_target_chain
        elif chain.id == old_peptide_chain:
            chain.id = desired_peptide_chain

    # write out the new PDB with suffix -renamed.pdb
    io = PDBIO()
    io.set_structure(structure)
    filename = pdb_filepath[:-4] + "-renamed.pdb"
    io.save(filename)
