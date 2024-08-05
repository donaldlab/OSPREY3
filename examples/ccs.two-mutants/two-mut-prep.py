
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
    omol_path = '6ov7.omol'
    open(omol_path, 'w').write(osprey.prep.saveOMOL(mols))
    print('saved prepared OMOL to %s' % omol_path)

    # optional: save the cleaned PDB
    # pdb_path = '6ov7-python.pdb'
    # open(pdb_path, 'w').write(osprey.prep.savePDB(mols))
    # print('saved prepared PDB to %s' % pdb_path)

print('Molecule preparation complete!')


# let's create a new conformation space from an OMOL file created by a previous molecule preparation python script
mols_path = '6ov7.omol'
mols = osprey.prep.loadOMOL(open(mols_path, 'r').read())

# the molecules we prepped are CALP and kCAL01
CALP = mols[0]
kCAL01 = mols[1]

# create a new conformation space out of the molecules
# we'll define all of our sequence mutations in the conformation space,
# as well as all of the desired conformations for each sequence,
# desired flexibilities, etc
conf_space = osprey.prep.ConfSpace(mols)

# choose a name for your conformation space
conf_space.setName('Conformation Space')

# Conformation Space Preparation Step 1: load the conformation libraries
# The conformation libraries are collections of mutations and conformations for molecules.
# In this case, the Lovell et al 2000 library describes common rotamers for protein natual amino acids.
# We'll load it from Osprey's built-in database of conformation libraries by scanning for its id.
lovell2000 = next(lib for lib in osprey.prep.confLibs if lib.getId() == 'lovell2000-osprey3').load()
conf_space.getConflibs().add(lovell2000)


# Conformation Space Preparation Step 2: define WT mutations
# Create a design position for that residue, then add the mutations.
# The identifiers for the mutations must match the fragment types in the loaded conformation libraries.

# add kCAL01 mutants (Chain C w/ mutations)
pos6gln = conf_space.addPosition(osprey.prep.ProteinDesignPosition(kCAL01, 'C', '6'))
conf_space.addMutations(pos6gln, 'GLN', 'ALA', 'GLU')

pos10val = conf_space.addPosition(osprey.prep.ProteinDesignPosition(kCAL01, 'C', '10'))
conf_space.addMutations(pos10val, 'VAL')

# add CALP mutants (Chain A, only WT)
pos295ile = conf_space.addPosition(osprey.prep.ProteinDesignPosition(CALP, 'A', '295'))
conf_space.addMutations(pos295ile, 'ILE')

pos341his = conf_space.addPosition(osprey.prep.ProteinDesignPosition(CALP, 'A', '341'))
conf_space.addMutations(pos341his, 'HIS')

pos345val = conf_space.addPosition(osprey.prep.ProteinDesignPosition(CALP, 'A', '345'))
conf_space.addMutations(pos345val, 'VAL')

# now we have defined a conformation space
print('conformation space describes %s sequences:' % conf_space.countSequences())
for pos in conf_space.positions():
    print('\t%6s mutations: %s' % (pos.getName(), conf_space.getMutations(pos)))


# Conformation Space Preparation Step 3: define discrete flexibility using conformations
# Conformations tell osprey how to model structural flexibility when evaluating sequences
for pos in conf_space.positions():

    # First, add conformations from the conformation libraries for our mutations.
    for mutation in conf_space.getMutations(pos):
        conf_space.addConformationsFromLibraries(pos, mutation)

    # Next, if the design position allows the wild-type "mutation"
    if pos.getType() in conf_space.getMutations(pos):
        # Create a conformation out of the PDB structure (the "wild-type" input structure conformation) and add it
        # Often these "wild-type" conformations will have lower energies than the library conformations,
        # since the "wild-type" conformations are specifically adapted to the input structures.
        conf_space.addWildTypeConformation(pos)

# now we have a larger conformation space of 32,000 conformations
print('conformation space describes %s conformations:' % conf_space.countConformations())
for pos in conf_space.positions():
    print('\t%6s conformations:' % pos.getName())
    for frag in conf_space.getFragments(pos):
        conf_ids = [conf_info.getConf().getId() for conf_info in conf_space.getConformations(pos, frag)]
        print('\t\tfragment %10s: %s' % (frag.getId(), conf_ids))


# Conformation Space Preparation Step 4: define continuous flexiblilty and trans/rot
# add dihedral angle motions to all the conformations at each design position
dihedral_settings = osprey.prep.DihedralAngleSettings()
for pos in conf_space.positions():
    for mutation in conf_space.getMutations(pos):
        for conf_info in conf_space.getConformations(pos, mutation):
            for motion in osprey.prep.conformationDihedralAngles(pos, conf_info, dihedral_settings):
                conf_info.getMotions().add(motion)

# add a translation/rotation motion to kCAL01 (adding to CALP is redundant)
conf_space.addMotion(osprey.prep.moleculeTranslationRotation(kCAL01))

# Conformation Space Preparation Step 5: save the conformation space
# Save the complex
path = 'complex.confspace'
open(path, 'w').write(osprey.prep.saveConfSpace(conf_space))
print('saved complex conformation space to %s' % path)

# Save the ligand
conf_space_kCAL01 = conf_space.copy(kCAL01)
# conf_space_kCAL01.setName('kCAL01')
path = 'kCAL01.confspace'
open(path, 'w').write(osprey.prep.saveConfSpace(conf_space_kCAL01))
print('saved kCAL01 conformation space to %s' % path)

# save the target protein
conf_space_CALP = conf_space.copy(CALP)
# conf_space_CALP.setName('CALP')
path = 'CALP.confspace'
open(path, 'w').write(osprey.prep.saveConfSpace(conf_space_CALP))
print('saved CALP conformation space to %s' % path)

# Conformation Space Preparation Step 6: compile the conformation spaces
# "Compilation", in this case, is the process by which Osprey transforms the
# conformation space(s) defined here into an optimized binary file format that
# will can be efficiently used later by a design algorithm.
# This "compilation" process does many things, but most importantly it uses
# AmberTools to gather the forcefield parameters for your molecules and all the mutations
# and conformations in the conformation space.
def compile(cs, name):

    compiler = osprey.prep.ConfSpaceCompiler(cs)

    # add the forcefields to use
    compiler.getForcefields().add(osprey.prep.Forcefield.Amber96)
    compiler.getForcefields().add(osprey.prep.Forcefield.EEF1)

    # run the compiler and wait for it to finish
    print('compiling %s ...' % name)
    progress = compiler.compile()
    progress.printUntilFinish(5000)
    report = progress.getReport()

    # check for compiler errors
    if report.getError() is not None:
        raise Exception('Compilation failed', report.getError())

    # save the compiled conf space
    path = '%s.ccsx' % name
    print('saving %s to %s ...' % (name, path))
    open(path, 'wb').write(osprey.prep.saveCompiledConfSpace(report.getCompiled()))


# start the local service that calls AmberTools for us
# NOTE: this will only work on Linux machines
with osprey.prep.LocalService():

    # compile all three conformation spaces
    compile(conf_space, 'complex')
    compile(conf_space_kCAL01, 'kCAL01')
    compile(conf_space_CALP, 'CALP')


print('Conformation space preparation complete!')