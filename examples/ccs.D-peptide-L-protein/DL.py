from Bio.PDB import PDBParser, NeighborSearch
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

# create a list of dict of flexible residues for each IAS round
flexible_set = []
for r in model[peptide_chain_id].get_residues():
    # get the res type, id for pep res nearest neighbors
    # we'll search within 6A of the CA of the current peptide residue
    target_atoms = [atom for atom in model[target_chain_id].get_atoms()]
    search_from = r['CA']
    ns = NeighborSearch(target_atoms)
    nearby_res = ns.search(search_from.coord, 6.0, level='R')

    # create a dict list of target's flexible residues
    # key = AA type, value = res #
    flexing = []
    for n in nearby_res:
        kv = {}
        kv[n.get_resname()] = str(n.id[1])
        flexing.append(kv)

    # output the search results
    print('Found flexible set for peptide residue %s %s:' % (r.get_resname(), r.id[1]))
    print(flexing)

    # add to flexible set list
    flexible_set.append(flexing)

# keep track of peptide atom index for later file naming
peptide_residues = []
for p in model[peptide_chain_id].get_residues():
    name = p.get_resname()
    index = str(p.id[1])
    pep_res = name + index
    peptide_residues.append(pep_res)

# calculate the # of IAS rounds required. -1 so we can index peptide_residues
curr_IAS_round = 0
num_IAS_rounds = len(flexible_set) - 1

# we have our L-target and D-peptide from DL-preprocessing
pdb = osprey.prep.loadPDB(open(pdb_path, 'r').read())

# get both protein chains
target = pdb[0]
peptide = pdb[1]
mols = [target, peptide]

# PDB doesn't save bonds in OMOL format, so insert here
with osprey.prep.LocalService():
    for mol in mols:
        bonds = osprey.prep.inferBonds(mol)
        for bond in bonds:
            mol.getBonds().add(bond)
        print('added %d bonds to %s' % (len(bonds), mol))

# save OMOL for the L-target
target_omol_path = pdb_path[:-15] + 'target.omol'
open(target_omol_path, 'w').write(osprey.prep.saveOMOL([target]))
print('saved prepared OMOL to %s' % target_omol_path)

# save OMOL for the D-peptide
peptide_omol_path = pdb_path[:-15] + 'peptide.omol'
open(peptide_omol_path, 'w').write(osprey.prep.saveOMOL([peptide]))
print('saved prepared OMOL to %s' % peptide_omol_path)

# rename the D-peptide.omol to be molecule.1
# this will be important later when we paste together the files
f = open(peptide_omol_path, 'r')
filedata = f.read()
f.close()

newatom = filedata.replace("molecule.0", "molecule.1")
newpoly = newatom.replace("molecule.0.polymer", "molecule.1.polymer")

f = open(peptide_omol_path, 'w')
f.write(newpoly)
f.close()

for s in flexible_set:

    # open the omol and prep the confspace for specification
    target_omol = osprey.prep.loadOMOL(open(target_omol_path, 'r').read())
    target_conf = target_omol[0]
    target_conf_space = osprey.prep.ConfSpace(target_omol)
    lovell2000 = next(lib for lib in osprey.prep.confLibs if lib.getId() == 'lovell2000-osprey3').load()
    target_conf_space.getConflibs().add(lovell2000)

    # define mutations for all flexible residues in the mut scan dictionary list
    for r in s:
        for identity, id in r.items():
            target_mut = target_conf_space.addPosition(osprey.prep.ProteinDesignPosition(target_conf, target_chain_id, id))
            target_conf_space.addMutations(target_mut, identity)

    # report the target conf space for this IAS round
    print('conformation space describes %s sequences:' % target_conf_space.countSequences())
    for pos in target_conf_space.positions():
        print('\t%6s mutations: %s' % (pos.getName(), target_conf_space.getMutations(pos)))

    # add discrete flexibility to the round
    for pos in target_conf_space.positions():

        # add conformations from library
        for mutation in target_conf_space.getMutations(pos):
            target_conf_space.addConformationsFromLibraries(pos, mutation)

        # add WT flexibility
        if pos.getType() in target_conf_space.getMutations(pos):
            target_conf_space.addWildTypeConformation(pos)

    # print current confspace
    print('conformation space describes %s conformations:' % target_conf_space.countConformations())
    for pos in target_conf_space.positions():
        print('\t%6s conformations:' % pos.getName())
        for frag in target_conf_space.getFragments(pos):
            conf_ids = [conf_info.getConf().getId() for conf_info in target_conf_space.getConformations(pos, frag)]
            print('\t\tfragment %10s: %s' % (frag.getId(), conf_ids))


    # define continuous flexiblilty and trans/rot
    dihedral_settings = osprey.prep.DihedralAngleSettings()
    for pos in target_conf_space.positions():
        for mutation in target_conf_space.getMutations(pos):
            for conf_info in target_conf_space.getConformations(pos, mutation):
                for motion in osprey.prep.conformationDihedralAngles(pos, conf_info, dihedral_settings):
                    conf_info.getMotions().add(motion)


    # for each round, save the target confspace w/ a descriptive filename
    curr_residue = peptide_residues[curr_IAS_round]
    print("curr res is " + curr_residue)
    save_confspace_path = pdb_path[:-15] + 'target-' + curr_residue + '.confspace'
    open(save_confspace_path, 'w').write(osprey.prep.saveConfSpace(target_conf_space))
    print('saved complex conformation space to %s' % save_confspace_path)

    # compile the target confspace files for each IAS round
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
        save_compile_name = pdb_path[:-15] + 'target-' + curr_residue
        compile(target_conf_space, save_compile_name)

        # index the current IAS round
        curr_IAS_round += 1

    print('Conformation space preparation complete!')

# TODO: prepare peptide confspace