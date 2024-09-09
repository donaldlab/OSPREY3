# this file generates D:L complexes from MASTER returns
# assumes L-protein is chain z, L-peptide is chain y
# MASTER-returned L-peptide can be any id
import warnings

from Bio.PDB import *

import os

# target_pdb is the L-L complex
# input_directory is the MASTER matches returns (all PDBs)
# output_directory is where to put the D-L complexes

# hide warnings
warnings.simplefilter('ignore')

def scaffold_generator(complex_pdb: str, input_directory: str, output_directory: str):

    for pdb_path in os.listdir(input_directory):

        # get the L-complex
        parser = PDBParser(PERMISSIVE=1)
        complex_structure = parser.get_structure('complex', complex_pdb)
        complex_model = complex_structure[0]

        # get the old (L-space) peptide
        old_peptide = complex_model["y"]

        # flip each MASTER L-match to D-space, and change ID if needed
        full_pdb_path = os.path.join(input_directory, pdb_path)
        M_pep_structure = parser.get_structure('match', full_pdb_path)

        M_peptide_model = M_pep_structure[0]
        for chain in M_peptide_model:
            if chain.id == "z" or chain.id == "y":
                chain.id = "A"
                print("%s had the same chain ID as the complex, so changed to A" % pdb_path)
            for res in chain:
                for a in res:
                    a.coord[2] = a.coord[2] * -1

        # align the old L-pep to the MASTER L-pep
        # get the CA from both
        old_atoms = []
        M_atoms = []
        for res in old_peptide:
            for a in res:
                if a.name == 'CA':
                    old_atoms.append(a)
        for chain in M_peptide_model:
            for res in chain:
                for a in res:
                    if a.name == 'CA':
                        M_atoms.append(a)

        imposer = Superimposer()
        imposer.set_atoms(old_atoms, M_atoms)
        imposer.apply(M_peptide_model.get_atoms())

        # add to new D-pep to the PDB
        for chain in M_peptide_model:
            complex_model.add(chain)

        # remove the old peptide from the complex
        complex_model.detach_child("y")

        # print new PDB to output_directory
        io = PDBIO()
        io.set_structure(complex_model)
        pdb_name = pdb_path[:-4] + "-complex.pdb"
        savename = os.path.join(output_directory, pdb_name)
        io.save(savename)
