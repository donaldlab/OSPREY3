
import osprey
osprey.start()

# import the prep module after starting Osprey
import osprey.prep

# let's create a new conformation space from an OMOL file created from the GUI
# OMOL only contains the D-ligand
mols_path = 'GUImatch.omol'
mols = osprey.prep.loadOMOL(open(mols_path, 'r').read())

# get the D-ligand
match = mols[0]

# create a new conformation space out of the protein(s)
conf_space = osprey.prep.ConfSpace(mols)

# choose a name for your conformation space
conf_space.setName('Conformation Space')

# Conformation Space Preparation Step 1: load the conformation libraries
# We'll load it from Osprey's built-in database of conformation libraries by scanning for its id.
Dlovell2000 = next(lib for lib in osprey.prep.confLibs if lib.getId() == 'D-lovell2000-osprey3').load()
conf_space.getConflibs().add(Dlovell2000)

# Conformation Space Preparation Step 2: define WT mutations
# Create a design position for that residue, then add the mutations.
# The identifiers for the mutations must match the fragment types in the loaded conformation libraries!

# add ligand mutants (Chain A w/ mutations)
res48 = conf_space.addPosition(osprey.prep.ProteinDesignPosition(match, 'A', '48'))
conf_space.addMutations(res48, 'ILE', 'ARG')

res49 = conf_space.addPosition(osprey.prep.ProteinDesignPosition(match, 'A', '49'))
conf_space.addMutations(res49, 'SER', 'HIS')

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

# now we have a larger conformation space
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
conf_space.addMotion(osprey.prep.moleculeTranslationRotation(match))

# Conformation Space Preparation Step 5: save the conformation space
# Save the D-ligand
# conf_space.setName('Match4')
path = 'PINTmatch.confspace'
open(path, 'w').write(osprey.prep.saveConfSpace(conf_space))
print('saved complex conformation space to %s' % path)

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
    compile(conf_space, 'PINTmatch')


print('Conformation space preparation complete!')