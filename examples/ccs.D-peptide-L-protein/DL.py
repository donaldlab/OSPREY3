from Bio.PDB import PDBParser, NeighborSearch
from Confspace_Combiner import *
from Slurm_Maker import *
from pathlib import Path
import os
import osprey
osprey.start()
import osprey.prep

# this files expects a properly-prepared D-L complex PDB file. Use DL-preprocess.py for rapid preprocessing.
# this prepares OSPREY confspace/ccsx files for 1 round of IAS for each residue on the D-peptide.
# a new directory is created for each D-L complex PDB, with subdirectories for each IAS round.

# TODO: user must provide the following settings:
# provide the directory (with slash) where ONLY prepared D-L complex PDBs are stored
input_directory = 'prepared-PDBs/'

# provide slurm settings
slurm_mem = 750
slurm_cpus = 48
slurm_partition = "grisman"
kstar_epsilon = 0.05



for pdb_path in os.listdir(input_directory):

    original_pdb_path = os.path.join(input_directory, pdb_path)

    # get the match num + create the directory for saving outputs
    match_name = ""
    for s in pdb_path:
        if s != "-":
            match_name += s
        elif s == "-":
            break

    print("\n\n----------- now preparing " + match_name + " -----------\n\n")

    output_directory = match_name + '/'
    Path(output_directory).mkdir(exist_ok=True)
    print("Created new output parent directory: " + output_directory)

    # get target and peptide chain ids for later confspace specification
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("complex", original_pdb_path)
    model = structure[0]
    chain_ids = []
    for chain in model:
        chain_ids.append(chain.id)
    target_chain_id = chain_ids[0]
    peptide_chain_id = chain_ids[1]

    # get a list of the pep residue numbers
    peptide_set = []
    for r in model[peptide_chain_id]:
        peptide_set.append(str(r.id[1]))

    # create a list of dict of flexible residues for each IAS round
    flexible_set = []

    for r in model[peptide_chain_id].get_residues():

        # get the res type, id for pep res nearest neighbors
        # we'll search within 6A of the CA of the current peptide residue
        tatoms = [atom for atom in model[target_chain_id].get_atoms()]
        search_from = r['CA']
        ns = NeighborSearch(tatoms)
        nearby_res = ns.search(search_from.coord, 6.0, level='R')

        # create a dict list of target's flexible residues
        # key = AA type, value = res #
        flexing = []
        for n in nearby_res:

            # GLY and ALA are small residues, so exclude from flexible set
            if n.get_resname() == "GLY" or n.get_resname() == "ALA":
                continue

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

    # IAS round tracker. curr is very important for file naming and opening.
    curr_IAS_round = 0
    num_IAS_rounds = len(flexible_set)

    # we have our prepared L-target and D-peptide from DL-preprocessing
    pdb = osprey.prep.loadPDB(open(original_pdb_path, 'r').read())

    # get both protein chains
    target = pdb[0]
    peptide = pdb[1]
    mols = [target, peptide]

    # PDB doesn't save bonds in OMOL format (can't do in preprocessing), so insert here
    with osprey.prep.LocalService():
        for mol in mols:
            bonds = osprey.prep.inferBonds(mol)
            for bond in bonds:
                mol.getBonds().add(bond)
            print('added %d bonds to %s' % (len(bonds), mol))

    # save OMOL for the L-target
    target_omol_path = output_directory + pdb_path[:-15] + 'target.omol'
    open(target_omol_path, 'w').write(osprey.prep.saveOMOL([target]))
    print('saved prepared OMOL to %s' % target_omol_path)

    # save OMOL for the D-peptide
    peptide_omol_path = output_directory + pdb_path[:-15] + 'peptide.omol'
    open(peptide_omol_path, 'w').write(osprey.prep.saveOMOL([peptide]))
    print('saved prepared OMOL to %s' % peptide_omol_path)

    # define the confspace compiler for the target, peptide, and complex
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
        path = '%s' % name
        print('saving %s to %s ...' % (name, path))
        open(path, 'wb').write(osprey.prep.saveCompiledConfSpace(report.getCompiled()))

    # prepare each round of IAS for the target (L-space)
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

            # add WT flexibility (this comes from the crystal structure conformation)
            if pos.getType() in target_conf_space.getMutations(pos):
                target_conf_space.addWildTypeConformation(pos)


        # define continuous flexiblilty
        dihedral_settings = osprey.prep.DihedralAngleSettings()
        for pos in target_conf_space.positions():
            for mutation in target_conf_space.getMutations(pos):
                for conf_info in target_conf_space.getConformations(pos, mutation):
                    for motion in osprey.prep.conformationDihedralAngles(pos, conf_info, dihedral_settings):
                        conf_info.getMotions().add(motion)


        # for each round, save the target confspace in a new subdirectory w/ a descriptive filename
        # subdirectories are the name of the corresponding pep residue
        curr_residue = peptide_residues[curr_IAS_round]
        print("current residue is " + curr_residue)
        new_subdirectory = output_directory + curr_residue + "/"
        Path(new_subdirectory).mkdir(exist_ok=True)
        target_confspace_path = new_subdirectory + pdb_path[:-15] + 'target' + '.confspace'
        open(target_confspace_path, 'w').write(osprey.prep.saveConfSpace(target_conf_space))
        print('saved target conformation space to %s' % target_confspace_path)

        # compile the confspace
        # NOTE: this only works on linux systems
        with osprey.prep.LocalService():

            # compile the target confspace
            target_ccsx_path = target_confspace_path[:-9] + 'ccsx'
            compile(target_conf_space, target_ccsx_path)

            # index the current IAS round
            curr_IAS_round += 1

        print('\n\n' + match_name + ': Completed round %s / %s for target\n\n' % (curr_IAS_round, num_IAS_rounds))

    # prepare the peptide conformation spaces
    # for each residue in the peptide, create a confspace that assigns all 20 AA to the residue

    # reset the IAS rounder counter
    curr_IAS_round = 0

    # declare list of 19 AA for reference
    # excludes PRO for lack of N/C-term compatibility
    # preprocessing already filters peptides with WT of N/C-PRO
    amino_acids = ['GLY', 'ALA', 'VAL', 'CYS', 'LEU', 'ILE', 'MET', 'TRP', 'PHE',
                   'LYS', 'ARG', 'HIS',
                   'SER', 'THR', 'TYR', 'ASN', 'GLN',
                   'ASP', 'GLU']

    # prepare the confspace for each peptide residue, then create the complex confspace files
    for s in peptide_set:

        # open the omol and prep the confspace for specification
        peptide_omol = osprey.prep.loadOMOL(open(peptide_omol_path, 'r').read())
        peptide_conf = peptide_omol[0]
        peptide_conf_space = osprey.prep.ConfSpace(peptide_omol)

        # since this is a D-peptide, add the D AA library
        Dlovell2000 = next(lib for lib in osprey.prep.confLibs if lib.getId() == 'D-lovell2000-osprey3').load()
        peptide_conf_space.getConflibs().add(Dlovell2000)

        # add the 19 AA identities to the 'mutations' confspace (WT is not really a mutant)
        peptide_mut = peptide_conf_space.addPosition(osprey.prep.ProteinDesignPosition(peptide_conf, peptide_chain_id, s))
        peptide_conf_space.addMutations(peptide_mut, amino_acids)

        # report the peptide conf space for this IAS round
        print('conformation space describes %s sequences:' % peptide_conf_space.countSequences())
        for pos in peptide_conf_space.positions():
            print('\t%6s mutations: %s' % (pos.getName(), peptide_conf_space.getMutations(pos)))

        # add discrete flexibility to the round
        for pos in peptide_conf_space.positions():

            # add conformations from library
            for mutation in peptide_conf_space.getMutations(pos):
                peptide_conf_space.addConformationsFromLibraries(pos, mutation)

            # add WT flexibility
            if pos.getType() in peptide_conf_space.getMutations(pos):
                peptide_conf_space.addWildTypeConformation(pos)

        # define continuous flexiblilty
        dihedral_settings = osprey.prep.DihedralAngleSettings()
        for pos in peptide_conf_space.positions():
            for mutation in peptide_conf_space.getMutations(pos):
                for conf_info in peptide_conf_space.getConformations(pos, mutation):
                    for motion in osprey.prep.conformationDihedralAngles(pos, conf_info, dihedral_settings):
                        conf_info.getMotions().add(motion)

        # add translation/rotation motion
        # adding this to target would be redundant (and more comp expensive)
        peptide_conf_space.addMotion(osprey.prep.moleculeTranslationRotation(peptide_conf))

        # for each round, save the peptide confspace w/ a descriptive filename in same subdirectory as target
        curr_residue = peptide_residues[curr_IAS_round]
        print("current residue is " + curr_residue)
        target_subdirectory = output_directory + curr_residue + '/'
        peptide_confspace_path = target_subdirectory + pdb_path[:-15] + 'peptide' + '.confspace'
        open(peptide_confspace_path, 'w').write(osprey.prep.saveConfSpace(peptide_conf_space))
        print('saved peptide conformation space to %s' % peptide_confspace_path)

        # compile the peptide residue confspace
        with osprey.prep.LocalService():
            peptide_ccsx_path = peptide_confspace_path[:-9] + 'ccsx'
            compile(peptide_conf_space, peptide_ccsx_path)
        print('\n\n' + match_name + ': Completed round %s / %s for peptide\n\n' % ((curr_IAS_round + 1), num_IAS_rounds))

        # combine the pep residue with its respective target IAS round
        # curr_residue determines corresponding confspace
        target_confspace_path = target_subdirectory + pdb_path[:-15] + 'target' + '.confspace'
        complex_confspace_path = target_subdirectory + pdb_path[:-15] + 'complex' + '.confspace'
        complex_info = combine_Confspaces(target_confspace_path, peptide_confspace_path, complex_confspace_path)

        # load confspace
        complex_confspace = osprey.prep.loadConfSpace(open(complex_confspace_path, 'r').read())

        # compile the complex confspace
        with osprey.prep.LocalService():

            # compile the complex
            complex_ccsx_path = complex_confspace_path[:-9] + 'ccsx'
            compile(complex_confspace, complex_ccsx_path)

        # insert the slurm script and an empty ensemble directory (for Kstar outputs)
        ensemble_dir = target_subdirectory + 'ensembles/'
        Path(ensemble_dir).mkdir(exist_ok=True)
        slurm_info = make_slurm(match_name, slurm_mem, slurm_cpus, slurm_partition, kstar_epsilon)
        slurm_filename = target_subdirectory + match_name + '-' + curr_residue + '.sh'
        with open(slurm_filename, 'w') as file:
            file.write(slurm_info)

        # index the current IAS round
        curr_IAS_round += 1

        print('\n\n' + match_name + ': Completed round %s / %s for complex\n\n' % (curr_IAS_round, num_IAS_rounds))

    print('completed IAS preparations successfully')