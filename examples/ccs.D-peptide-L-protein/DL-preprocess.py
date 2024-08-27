import os
from Bio.PDB import PDBParser, PDBIO

# this script is to be run before DL.py
# this handles PDB prep including atomic labelling, bonds, protonation, minimization, and inversions
# this supposes a D-peptide in complex with an L-target, with the L-target listed 1st in the PDB (index 0)
# modify the index for the ids list to change these assumptions

# design preprocess step one: change atomic labelling
# commonly, PDB atom labels for carbon don't match ambertools templates. Let's change this manually to avoid
# any issues with adding missing atoms

# this is the directory where the MASTER-returned D:L complexes are stored. Forward slash is required.
directory = "scaffolds/"

for f in os.listdir(directory):
    file = os.path.join(directory, f)

    # create the required BioPython modules so we can work with the PDB
    parser = PDBParser(PERMISSIVE=1)
    structure_id = 'complex'
    filename = file
    structure = parser.get_structure(structure_id, filename)

    # get 1st PDB in the ensemble
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

    # design preprocess step two: flip the D-peptide to L-space, so we can work in ambertools w/ a homochiral system
    # get the chain names (note: assumes L-target is listed 1st in the PDB)
    ids = []
    for chain in model:
        ids.append(chain.id)
    target_id = ids[0]
    peptide_id = ids[1]

    # assign chains
    target = model[target_id]
    peptide = model[peptide_id]

    # flip the D-peptide to a L-peptide
    # this separates the chains in space, so we don't need to worry about chiral-induced clashes when prepping (step 4)
    for r in peptide:
        for a in r:
            a.coord[2] = a.coord[2] * -1

    # save the L:L complex with a different suffix
    io = PDBIO()
    io.set_structure(structure)
    L_filename = directory + f[:-4] + "-L.pdb"
    # io.save(L_filename)

    # design preprocess step 4: prepare the L-ligand and L-target
    # ambertools has trouble with D-peptide protonation, so we'll prepare everything in L-space, then flip back
    # PDB Parser from BioPython automatically checks for duplicate atoms, so we can skip this step

    # add missing heavy atoms

    # add bonds

    # strip Hs. We don't know if the PDB depositor used reliable protonation software, so let's remove and use our own.
    # it's also very possible no H are present, but it's a good idea to check

    # add Hs

    # flip the L-peptide back to D-space

    # save the D-L complex

    # design preprocess step two: minimize the D-ligand wrt the L-target
    # SANDER requires the chemical context of the D-L context for minimization, so we must pass the complex

    # save the minimized D:L complex