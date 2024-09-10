# this file generates D:L complexes from MASTER returns
# assumes L-protein is chain z, L-peptide is chain y
# MASTER-returned L-peptide can be any id
import warnings

from Bio.PDB import *

import os

# target_pdb is the L-L complex
# input_directory is the MASTER matches returns (all PDBs)
# output_directory is where to put the D-L complexes

# hide warnings - throws for a lot of atom names, but ok
warnings.simplefilter('ignore')


def scaffold_generator(complex_pdb: str, input_directory: str, output_directory: str):

    for pdb_path in os.listdir(input_directory):

        # get the L-complex
        parser = PDBParser(PERMISSIVE=1)
        complex_structure = parser.get_structure('complex', complex_pdb)
        complex_model = complex_structure[0]

        # get the old (L-space) peptide
        old_peptide = complex_model["y"]

        # get the MASTER L-peptide
        full_pdb_path = os.path.join(input_directory, pdb_path)
        M_pep_structure = parser.get_structure('match', full_pdb_path)
        M_peptide_model = M_pep_structure[0]

        # skip the match if it has a disjoint segment
        has_disjoint = False
        last_res = 0
        for chain in M_peptide_model:
            for res in chain:
                if last_res == 0:
                    last_res = res.id[1]
                elif (res.id[1] - last_res) != 1:
                    has_disjoint = True
                    print("SKIP: %s has a disjoint segment and will not be prepared" % pdb_path)
                    break
                else:
                    last_res = res.id[1]
        if has_disjoint:
            continue

        # flip each MASTER L-match to D-space, and change ID if needed
        for chain in M_peptide_model:
            if chain.id == "z" or chain.id == "y":
                chain.id = "A"
                print("ALTERATION: %s had the same chain ID as the complex, so changed ID to A" % pdb_path)
            for res in chain:
                for a in res:
                    a.coord[2] = a.coord[2] * -1

        # MASTER already aligns, so just add the D-pep
        for chain in M_peptide_model:
            complex_model.add(chain)

        # remove the old L-peptide from the complex
        complex_model.detach_child("y")

        # print new PDB to output_directory
        io = PDBIO()
        io.set_structure(complex_model)
        pdb_name = pdb_path[:-4] + "-complex.pdb"
        savename = os.path.join(output_directory, pdb_name)
        io.save(savename)

        # update terminal
        print("------ completed %s ------" % pdb_path)


# example call
# scaffold_generator("8gaj-trim-renamed.pdb", "matches", "scaffolds")