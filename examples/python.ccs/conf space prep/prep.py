
import osprey
osprey.start()

# import the prep module after starting Osprey
import osprey.prep


# let's create a new conformation space from an OMOL file created by a previous molecule preparation
mols_path = '1dg9.omol'
mols = osprey.prep.loadOMOL(open(mols_path, 'r').read())

# the molecules we prepped are Bovine PTPase and HEPES
ptpase = mols[0]
hepes = mols[1]

# create a new conformation space out of the molecules
# we'll define all of our sequence mutations in the conformation space,
# as well as all of the desired conformations for each sequence,
# desired flexibilities, etc
conf_space = osprey.prep.ConfSpace(mols)

# choose a name for your conformation space
conf_space.setName('BPTPase and HEPES')


# Conformation Space Preparation Step 1: load the conformation libraries
# The conformation libraries are collections of mutations and conformations for molecules.
# In this case, the Lovell et al 2000 library describes common rotamers for protein natual amino acids.
# We'll load it from Osprey's built-in database of conformation libraries by scanning for its id.
lovell2000 = next(lib for lib in osprey.prep.confLibs if lib.getId() == 'lovell2000-osprey3').load()
conf_space.getConflibs().add(lovell2000)


# Conformation Space Preparation Step 2: define mutations

# allow the 13 LEU residue to either stay LEU, or mutate to a VAL
# Create a design position for that residue, the add the mutations.
# The identifiers for the mutations must match the fragment types in the loaded conformation libraries.
pos13leu = conf_space.addPosition(osprey.prep.ProteinDesignPosition(ptpase, 'A', '13'))
conf_space.addMutations(pos13leu, 'LEU', 'VAL')

# force the 16 ILE residue to mutate from ILE to either LEU or VAL
pos16ile = conf_space.addPosition(osprey.prep.ProteinDesignPosition(ptpase, 'A', '16'))
conf_space.addMutations(pos16ile, 'LEU', 'VAL')

# allow the 1 ALA residue to either stay ALA, or mutate to VAL
pos1ala = conf_space.addPosition(osprey.prep.ProteinDesignPosition(ptpase, 'A', '1'))
conf_space.addMutations(pos1ala, 'ALA', 'VAL')

# now we have defined a conformation space with 2x2x2=8 sequences in it
print('conformation space describes %s sequences:' % conf_space.countSequences())
for pos in conf_space.positions():
    print('\t%6s mutations: %s' % (pos.getName(), conf_space.getMutations(pos)))


# Conformation Space Preparation Step 3: define discrete flexiblilty using conformations
# Conformations tell osprey how to model structural flexibility when evaluating sequences
for pos in conf_space.positions():

    # First, add conformations from the conformation libraries for our mutations.
    for mutation in conf_space.getMutations(pos):
        conf_space.addConformationsFromLibraries(pos, mutation)

    # Next, if the design position allows the wild-type "mutation"
    if pos.getType() in conf_space.getMutations(pos):
        # Create a conformation out of the PDB structure (the "wild-type" conformation) and add it to the space.
        # Often these "wild-type" conformations will have lower energies than the library conformations,
        # since the "wild-type" conformations are specifically adapted to the input structures.
        conf_space.addWildTypeConformation(pos)

# now we have a conformation space with 9x8x5=360 conformations in it!
print('conformation space describes %s conformations:' % conf_space.countConformations())
for pos in conf_space.positions():
    print('\t%6s conformations:' % pos.getName())
    for frag in conf_space.getFragments(pos):
        conf_ids = [conf_info.getConf().getId() for conf_info in conf_space.getConformations(pos, frag)]
        print('\t\tfragment %10s: %s' % (frag.getId(), conf_ids))


# Conformation Space Preparation Step 4: define continuous flexiblilty using motions

# add dihedral angle motions to all the conformations at each design position
dihedral_settings = osprey.prep.DihedralAngleSettings()
for pos in conf_space.positions():
    for mutation in conf_space.getMutations(pos):
        for conf_info in conf_space.getConformations(pos, mutation):
            for motion in osprey.prep.conformationDihedralAngles(pos, conf_info, dihedral_settings):
                conf_info.getMotions().add(motion)

# add dihedral angles to HEPES too
conf_space.addMotion(osprey.prep.moleculeDihedralAngle(hepes, 'C10', 'S', 'O3S', 'H103', dihedral_settings))
conf_space.addMotion(osprey.prep.moleculeDihedralAngle(hepes, 'C9', 'C10', 'S', 'O3S', dihedral_settings))
conf_space.addMotion(osprey.prep.moleculeDihedralAngle(hepes, 'N1', 'C9', 'C10', 'S', dihedral_settings))
conf_space.addMotion(osprey.prep.moleculeDihedralAngle(hepes, 'C2', 'N1', 'C9', 'C10', dihedral_settings))
conf_space.addMotion(osprey.prep.moleculeDihedralAngle(hepes, 'C3', 'N4', 'C7', 'C8', dihedral_settings))
conf_space.addMotion(osprey.prep.moleculeDihedralAngle(hepes, 'N4', 'C7', 'C8', 'O8', dihedral_settings))
conf_space.addMotion(osprey.prep.moleculeDihedralAngle(hepes, 'C7', 'C8', 'O8', 'H81', dihedral_settings))

# add a translation/rotation motion to HEPES, to optimize the fit of HEPES to PTPase in different conformations
conf_space.addMotion(osprey.prep.moleculeTranslationRotation(hepes))


# Conformation Space Preparation Step 5: save the conformation space
# This saves all the work we just did to the conformation space to a file.
# Later, the file can be re-opened and changed to do something else.
# The conf space save file can even be opened and viewed/edited in the GUI.
path = 'complex.confspace'
open(path, 'w').write(osprey.prep.saveConfSpace(conf_space))
print('saved complex conformation space to %s' % path)

# If we're doing a K* type or other multi-state design,
# we'll need to split the complex conformation space into separate unbound states,
# assuming that we want to use the same PDB files for the bound and unbound states.
# But if you'd rather use a separate struture for the unbound states,
# you'll need to prepare the new structures and conformation spaces separately.
conf_space_ptpase = conf_space.copy(ptpase)
conf_space_ptpase.setName('BPTPase')
path = 'ptpase.confspace'
open(path, 'w').write(osprey.prep.saveConfSpace(conf_space_ptpase))
print('saved PTPase conformation space to %s' % path)

conf_space_hepes = conf_space.copy(hepes)
conf_space_hepes.setName('HEPES')
path = 'hepes.confspace'
open(path, 'w').write(osprey.prep.saveConfSpace(conf_space_hepes))
print('saved HEPES conformation space to %s' % path)


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
    compile(conf_space_ptpase, 'ptpase')
    compile(conf_space_hepes, 'hepes')


print('Conformation space preparation complete!')
