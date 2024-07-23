
import osprey
osprey.start()

# import the prep module after starting Osprey
import osprey.prep

# load PDB 6ov7 (kcal01, all L-amino acids)
# this has been manually trimmed in PyMOL to remove waters, unwanted dimer, and kcal01 to 6 residues
# note: filenames ending in "-GUI" indicate that operation was performed in the GUI interface
# note: filenames ending in "-python" indicate that operation was performed via the python scripts (like this one)
pdb_path = '6ov7.pdb'
pdb = osprey.prep.loadPDB(open(pdb_path, 'r').read())

# print confirm we have chains A and C (what we trimmed to in PyMOL)
print('Loaded %d molecules:' % len(pdb))
for mol in pdb:
    print('\t%s: %s' % (mol, osprey.prep.molTypes(mol)))

# looks like the PDB file has two protein chains (all L-space)
# chain is CALP, ligand is kCAL01
# optional: this can also be used to ignore water molecules instead of trim in PYMOL
CALP = pdb[0]
kCAL01 = pdb[1]
mols = [CALP, kCAL01]

# optional: give our molecules better names
# must be done manually on the omol file if using the GUI
CALP.setName('CALP')
kCAL01.setName('kCAL01')

# start the local service that calls AmberTools for us (Note: ONLY WORKS ON LINUX)
with osprey.prep.LocalService():

    # Molecule Preparation Step 1: remove duplicate atoms
    # For this PDB example, none should be present
    for mol in mols:
        # remove all but the first duplicated atom from each group
        for group in osprey.prep.duplicateAtoms(mol):
            for atomi in range(1, len(group.getAtoms())):
                group.remove(atomi)
                print('removed duplicate atom %s' % group)


    # Molecule Preparation Step 2: add missing heavy atoms
    # For this PDB, none should be added
    for mol in mols:
        for missing_atom in osprey.prep.inferMissingAtoms(mol):
            missing_atom.add()
            print('added missing atom: %s' % missing_atom)

    # Molecule Preparation Step 3: add bonds
    for mol in mols:
        bonds = osprey.prep.inferBonds(mol)
        for bond in bonds:
            mol.getBonds().add(bond)
        print('added %d bonds to %s' % (len(bonds), mol))

    # Molecule Preparation Step 4: add Hs
    # protonate both chains
    for mol in mols:
        protonated_atoms = osprey.prep.inferProtonation(mol)
        for protonated_atom in protonated_atoms:
            protonated_atom.add()
        print('added %d hydrogens to %s' % (len(protonated_atoms), mol))

    # Molecule Preparation Step 5: minimize the ligand
    # restrain the heavy atoms, but let H wander freely
    # change numSteps for more minimization (GUI files are done w/ 100 steps)
    # print('minimizing ...')
    #
    # def heavy_atoms(mol):
    #     return [atom for atom in mol.getAtoms() if atom.getElement().getSymbol() != 'H']
    #
    # osprey.prep.minimize(
    #     [osprey.prep.minimizerInfo(mol, heavy_atoms(mol)) for mol in mols],
    #     numSteps = 100
    # )
    #
    # print('minimization complete!')

    # Molecule Preparation Step 6: save the cleaned OMOL
    omol_path = '6ov7-python.omol'
    open(omol_path, 'w').write(osprey.prep.saveOMOL(mols))
    print('saved prepared OMOL to %s' % omol_path)

    # optional: save the cleaned PDB
    # pdb_path = '6ov7-python.pdb'
    # open(pdb_path, 'w').write(osprey.prep.savePDB(mols))
    # print('saved prepared PDB to %s' % pdb_path)

print('Molecule preparation complete!')