import osprey
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO

# this script is to be run before DL.py
# this handles common PDB preparation outside of OSPREY including atomic labeling and minimization

# design preprocess step one: change atomic labelling
# commonly, PDB atom labels for carbon don't match OSPREY templates. Let's change this manually to avoid
# any issues with python-automated molecule prep in DL.py
# f = open('match4.pdb', 'r')
# filedata = f.read()
# f.close()
#
# # MOST COMMON: change CD -> CD1 labelling
# # affected residues: ILE, LEU, PHE, TRP, TYR
# ILE = filedata.replace("CD  ILE", "CD1 ILE")
# LEU = ILE.replace("CD  LEU", "CD1 LEU")
# PHE = LEU.replace("CD  PHE", "CD1 PHE")
# TRP = PHE.replace("CD  TRP", "CD1 TRP")
# TYR = TRP.replace("CD  TYR", "CD1 TYR")
#
# f = open("match4-ready.pdb", 'w')
# f.write(TYR)
# f.close()

# design preprocess step two: minimize the D-ligand wrt the target



# # design preprocess step three: extract and flip the D-ligand to L-space
# # save the resulting D-ligand and L-target, now minimized and with correct atomic labels
# osprey.start()
#
# # import the prep module after starting Osprey
# import osprey.prep
#
# # load the complex PDB
# pdb_path = 'match4.pdb'
# pdb = osprey.prep.loadPDB(open(pdb_path, 'r').read())
#
# # separate the two chains
# target = pdb[0]
# ligand = pdb[1]
#
# with osprey.prep.LocalService():
#     with open('L-target.pdb', 'w') as file:
#         file.write(osprey.prep.savePDB(target))
#
#     with open('D-peptide.pdb', 'w') as file:
#         file.write(osprey.prep.savePDB(ligand))

# using biopython, create the necessary objects to access atomic coordinates of the D-peptide
parser = PDBParser(PERMISSIVE=1)
structure_id = 'peptide'
filename = "D-peptide.pdb"
structure = parser.get_structure(structure_id, filename)
model = structure[0]
chain = model["A"]

# loop over the atoms in each residue and flip over z axis
for r in chain:
    for a in r:
        a.coord[2] = a.coord[2] * -1

# we want to preserve the peptide binding location, so align to L-peptide to the D-peptide in complex

# save the inverted (now L-space + aligned) peptide
io = PDBIO()
io.set_structure(structure)
io.save("L-peptide.pdb")