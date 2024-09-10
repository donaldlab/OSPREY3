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


def chain_renamer(input_directory: str, output_directory: str, desired_target_id: str, desired_peptide_id: str):

    for f in os.listdir(input_directory):
        input_file = os.path.join(input_directory, f)

        # get the old IDs
        parser = PDBParser(PERMISSIVE=1)
        structure = parser.get_structure("complex", input_file)
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
        filename = f[:-4] + "-renamed.pdb"
        outfile = os.path.join(output_directory, filename)
        io.save(outfile)


# example call
# chain_renamer("test", "test2", "z", "y")