import osprey
from Bio.PDB import PDBParser, PDBIO

# this script is to be run before DL.py
# this handles PDB prep including atomic labelling, bonds, protonation, minimization, and inversions
# this supposes a D-peptide in complex with an L-target, with the L-target listed 1st in the PDB (index 0)
# modify the index for chain_ids in design process step 3 to change these assumptions

# design preprocess step one: change atomic labelling
# commonly, PDB atom labels for carbon don't match ambertools templates. Let's change this manually to avoid
# any issues with adding missing atoms
parser = PDBParser(PERMISSIVE=1)
structure_id = 'complex'
filename = "match1.pdb"
structure = parser.get_structure(structure_id, filename)
model = structure[0]

# MOST COMMON: change CD -> CD1 labelling
# affected residues: ILE, LEU, PHE, TRP, TYR
# loop over both chains + change CD labelling
for chain in model:
    for res in chain:
        if res.get_resname() in ("ILE", "LEU", "PHE", "TRP", "TYR"):
            for a in res:
                if a.fullname == ' CD ':
                    a.fullname = ' CD1'
                    print("changed residue " + str(res) + " to" + a.fullname)

# save the CD-corrected complex
io = PDBIO()
io.set_structure(structure)
io.save("match1-CD.pdb")

# design preprocess step two: minimize the D-ligand wrt the target



# # design preprocess step three: extract and flip the D-ligand to L-space
# # up until now, we've needed to preserve the D:L in the same PDB in order to generate
# # a correct mask for SANDER minimization. Now we can separate the ligand from the target.
parser = PDBParser(PERMISSIVE=1)
structure_id = 'complex'
# TODO: changed to post-min PDB
filename = "match1-CD.pdb"
structure = parser.get_structure(structure_id, filename)
model = structure[0]

# get the chain names (note: assumes L-target is listed 1st in the PDB)
ids = []
for chain in model:
    ids.append(chain.id)
target_id = ids[0]
peptide_id = ids[1]

# assign chains
target = model[target_id]
peptide = model[peptide_id]

# save the L-target
io = PDBIO()
io.set_structure(target)
io.save("match1-target.pdb")

# flip the D-peptide to a L-peptide
for r in peptide:
    for a in r:
        a.coord[2] = a.coord[2] * -1

# save the inverted (now L-space) peptide
io = PDBIO()
io.set_structure(peptide)
io.save("L-peptide.pdb")

# design preprocess step 4: prepare the L-ligand and L-target
# ambertools has trouble with D-peptides, so we'll prepare everything in L-space, then flip back
# remove duplicate atoms

# add missing heavy atoms

# add bonds

# strip Hs. We don't know if the PDB depositor used reliable protonation software, so let's remove and use our own.

# add Hs