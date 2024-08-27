from Bio.PDB import PDBParser, PDBIO
import osprey
osprey.start()
import osprey.prep

# TODO: for loop over prepared files
pdb_path = "match1-D-L-complex.pdb"

# get target and peptide chain ids for later confspace specification
parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure("complex", pdb_path)
model = structure[0]
chain_ids = []
for chain in model:
    chain_ids.append(chain.id)
target_chain_id = chain_ids[0]
peptide_chain_id = chain_ids[1]

# TODO: create IAS mutable residues dictionaries for target
# need: AA identity, residue number
# https://biopython.org/docs/1.76/api/Bio.PDB.NeighborSearch.html
thisdict = {
    'HIS': '301',
    'VAL': '303'
}

# we have our L-target and D-peptide from DL-preprocessing
pdb = osprey.prep.loadPDB(open(pdb_path, 'r').read())

# get both protein chains
target = pdb[0]
peptide = pdb[1]
mols = [target, peptide]

# save OMOL for the L-target
target_omol_path = 'match1-target.omol'
open(target_omol_path, 'w').write(osprey.prep.saveOMOL([target]))
print('saved prepared OMOL to %s' % target_omol_path)

# save OMOL for the D-peptide
peptide_omol_path = 'match1-peptide.omol'
open(peptide_omol_path, 'w').write(osprey.prep.saveOMOL([peptide]))
print('saved prepared OMOL to %s' % peptide_omol_path)

# # rename the D-peptide.omol to be molecule.1
# # this will be important later when we paste together the files
# f = open('match1-peptide.omol', 'r')
# filedata = f.read()
# f.close()
#
# newatom = filedata.replace("molecule.0", "molecule.1")
# newpoly = newatom.replace("molecule.0.polymer", "molecule.1.polymer")
#
# f = open("match1-peptide.omol", 'w')
# f.write(newpoly)
# f.close()

# let's create the target conformation space
target_omol = osprey.prep.loadOMOL(open(target_omol_path, 'r').read())
target_conf_space = osprey.prep.ConfSpace(target_omol)
lovell2000 = next(lib for lib in osprey.prep.confLibs if lib.getId() == 'lovell2000-osprey3').load()
target_conf_space.getConflibs().add(lovell2000)

# define mutations
# TODO: make this a for loop for length of IAS flex dict
for res in thisdict:
    res_type = res
    res_id = thisdict[res]
    print("target chain id is " + target_chain_id)
    print("residue type is " + res_type)
    print("id is " + res_id)
    # target_mut = target_conf_space.addPosition(osprey.prep.ProteinDesignPosition(target_omol, target_chain_id, res_id))
    # target_conf_space.addMutations(target_mut, res_type)

target_mut = target_conf_space.addPosition(osprey.prep.ProteinDesignPosition(target_omol, 'y', '301'))
target_conf_space.addMutations(target_mut, 'HIS')

# now we have defined a conformation space
print('conformation space describes %s sequences:' % target_conf_space.countSequences())
for pos in target_conf_space.positions():
    print('\t%6s mutations: %s' % (pos.getName(), target_conf_space.getMutations(pos)))
















# # Conformation Space Preparation Step 3: define discrete flexibility using conformations
# # Conformations tell osprey how to model structural flexibility when evaluating sequences
# for pos in conf_space.positions():
#
#     # First, add conformations from the conformation libraries for our mutations.
#     for mutation in conf_space.getMutations(pos):
#         conf_space.addConformationsFromLibraries(pos, mutation)
#
#     # Next, if the design position allows the wild-type "mutation"
#     if pos.getType() in conf_space.getMutations(pos):
#         # Create a conformation out of the PDB structure (the "wild-type" input structure conformation) and add it
#         # Often these "wild-type" conformations will have lower energies than the library conformations,
#         # since the "wild-type" conformations are specifically adapted to the input structures.
#         conf_space.addWildTypeConformation(pos)
#
# # now we have a larger conformation space of 32,000 conformations
# print('conformation space describes %s conformations:' % conf_space.countConformations())
# for pos in conf_space.positions():
#     print('\t%6s conformations:' % pos.getName())
#     for frag in conf_space.getFragments(pos):
#         conf_ids = [conf_info.getConf().getId() for conf_info in conf_space.getConformations(pos, frag)]
#         print('\t\tfragment %10s: %s' % (frag.getId(), conf_ids))
#
#
# # Conformation Space Preparation Step 4: define continuous flexiblilty and trans/rot
# # add dihedral angle motions to all the conformations at each design position
# dihedral_settings = osprey.prep.DihedralAngleSettings()
# for pos in conf_space.positions():
#     for mutation in conf_space.getMutations(pos):
#         for conf_info in conf_space.getConformations(pos, mutation):
#             for motion in osprey.prep.conformationDihedralAngles(pos, conf_info, dihedral_settings):
#                 conf_info.getMotions().add(motion)
#
# # add a translation/rotation motion to kCAL01 (adding to CALP is redundant)
# conf_space.addMotion(osprey.prep.moleculeTranslationRotation(kCAL01))
#
# # Conformation Space Preparation Step 5: save the conformation space
# # Save the complex
# path = 'complex.confspace'
# open(path, 'w').write(osprey.prep.saveConfSpace(conf_space))
# print('saved complex conformation space to %s' % path)
#
# # Save the ligand
# conf_space_kCAL01 = conf_space.copy(kCAL01)
# # conf_space_kCAL01.setName('kCAL01')
# path = 'kCAL01.confspace'
# open(path, 'w').write(osprey.prep.saveConfSpace(conf_space_kCAL01))
# print('saved kCAL01 conformation space to %s' % path)
#
# # save the target protein
# conf_space_CALP = conf_space.copy(CALP)
# # conf_space_CALP.setName('CALP')
# path = 'CALP.confspace'
# open(path, 'w').write(osprey.prep.saveConfSpace(conf_space_CALP))
# print('saved CALP conformation space to %s' % path)
#
# # Conformation Space Preparation Step 6: compile the conformation spaces
# # "Compilation", in this case, is the process by which Osprey transforms the
# # conformation space(s) defined here into an optimized binary file format that
# # will can be efficiently used later by a design algorithm.
# # This "compilation" process does many things, but most importantly it uses
# # AmberTools to gather the forcefield parameters for your molecules and all the mutations
# # and conformations in the conformation space.
# def compile(cs, name):
#
#     compiler = osprey.prep.ConfSpaceCompiler(cs)
#
#     # add the forcefields to use
#     compiler.getForcefields().add(osprey.prep.Forcefield.Amber96)
#     compiler.getForcefields().add(osprey.prep.Forcefield.EEF1)
#
#     # run the compiler and wait for it to finish
#     print('compiling %s ...' % name)
#     progress = compiler.compile()
#     progress.printUntilFinish(5000)
#     report = progress.getReport()
#
#     # check for compiler errors
#     if report.getError() is not None:
#         raise Exception('Compilation failed', report.getError())
#
#     # save the compiled conf space
#     path = '%s.ccsx' % name
#     print('saving %s to %s ...' % (name, path))
#     open(path, 'wb').write(osprey.prep.saveCompiledConfSpace(report.getCompiled()))
#
#
# # start the local service that calls AmberTools for us
# # NOTE: this will only work on Linux machines
# with osprey.prep.LocalService():
#
#     # compile all three conformation spaces
#     compile(conf_space, 'complex')
#     compile(conf_space_kCAL01, 'kCAL01')
#     compile(conf_space_CALP, 'CALP')
#
#
# print('Conformation space preparation complete!')