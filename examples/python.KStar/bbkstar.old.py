
import osprey
osprey.start()

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb('../../test-resources/4npd.A_DomainA_noH_trim_his_clean.min.pdb')

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# define the protein strand
protein = osprey.Strand(mol, templateLib=templateLib, residues=['A1', 'A37'])
protein.flexibility['A8'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# define the ligand strand
ligand = osprey.Strand(mol, templateLib=templateLib, residues=['A38', 'A58'])
ligand.flexibility['A41'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand.flexibility['A42'].setLibraryRotamers(osprey.WILD_TYPE,'THR','LYS').addWildTypeRotamers().setContinuous()
ligand.flexibility['A45'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand.flexibility['A46'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# make the conf space for the protein
proteinConfSpace = osprey.ConfSpace(protein)

# make the conf space for the ligand
ligandConfSpace = osprey.ConfSpace(ligand)

# make the conf space for the protein+ligand complex
complexConfSpace = osprey.ConfSpace([protein, ligand])

# how should we compute energies of molecules?
# (give the complex conf space to the ecalc since it knows about all the templates and degrees of freedom)
parallelism = osprey.Parallelism(cpuCores=4)
minimizingEcalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism, isMinimizing=True)

# BBK* needs a rigid energy calculator too, for multi-sequence bounds on K*
rigidEcalc = osprey.SharedEnergyCalculator(minimizingEcalc, isMinimizing=False)

# how should we define energies of conformations?
def confEcalcFactory(confSpace, ecalc):
    eref = osprey.ReferenceEnergies(confSpace, ecalc)
    return osprey.ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=eref)

epsilon = 0.99

# configure BBK*
bbkstar = osprey.BBKStar(
    proteinConfSpace,
    ligandConfSpace,
    complexConfSpace,
    numBestSequences=7,
    epsilon=0.68, # you proabably want something more precise in your real designs
    writeSequencesToConsole=True,
    writeSequencesToFile='bbkstar.results.tsv'
)

# configure BBK* inputs for each conf space
for info in bbkstar.confSpaceInfos():

    # how should we define energies of conformations?
    eref = osprey.ReferenceEnergies(info.confSpace, minimizingEcalc)
    info.confEcalcMinimized = osprey.ConfEnergyCalculator(info.confSpace, minimizingEcalc, referenceEnergies=eref)

    # compute the energy matrix
    emat = osprey.EnergyMatrix(info.confEcalcMinimized, cacheFile='emat.%s.dat' % info.id)

    # how should confs be ordered and searched? (don't forget to capture emat by using a defaulted argument)
    def makeAStar(rcs, emat=emat):
        return osprey.AStarTraditional(emat, rcs, showProgress=False)
    info.confSearchFactoryMinimized = osprey.KStar.ConfSearchFactory(makeAStar)

    # BBK* needs rigid energies and conformation search too
    rigidConfEcalc = osprey.ConfEnergyCalculatorCopy(info.confEcalcMinimized, rigidEcalc)
    rigidEmat = osprey.EnergyMatrix(rigidConfEcalc, cacheFile='emat.%s.rigid.dat' % info.id)
    def makeAStarRigid(rcs, emat=rigidEmat):
        return osprey.AStarTraditional(emat, rcs, showProgress=False)
    info.confSearchFactoryMinimized = osprey.KStar.ConfSearchFactory(makeAStarRigid)

    # Specify the input for the partition functions. Providing the confUpperBoundcalc turns on MARK*
    info.pfuncFactory = osprey.PartitionFunctionFactory(info.confSpace, info.confEcalcMinimized, info.id, confUpperBoundcalc=rigidConfEcalc)

# run BBK*
scoredSequences = bbkstar.run()

# use results
for scoredSequence in scoredSequences:
    print("result:")
    print("\tsequence: %s" % scoredSequence.sequence)
    print("\tscore: %s" % scoredSequence.score)
