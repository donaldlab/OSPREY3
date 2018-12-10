
import osprey
osprey.start()

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb('2RL0.min.reduce.pdb')

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# define the protein strand
proteinStrand = osprey.Strand(mol, templateLib=templateLib, residues=['G648', 'G654'])
proteinStrand.flexibility['G649'].setLibraryRotamers(osprey.WILD_TYPE, 'TYR', 'ALA', 'VAL', 'ILE', 'LEU').addWildTypeRotamers().setContinuous()
proteinStrand.flexibility['G650'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
proteinStrand.flexibility['G651'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
proteinStrand.flexibility['G654'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# define the ligand strand
ligandStrand = osprey.Strand(mol, templateLib=templateLib, residues=['A155', 'A194'])
ligandStrand.flexibility['A156'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligandStrand.flexibility['A172'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligandStrand.flexibility['A192'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligandStrand.flexibility['A193'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# define the MSK* states
protein = osprey.MSKStar_State('Protein', osprey.ConfSpace(proteinStrand))
ligand = osprey.MSKStar_State('Ligand', osprey.ConfSpace(ligandStrand))
complex = osprey.MSKStar_State('Complex', osprey.ConfSpace([proteinStrand, ligandStrand]))

# configure COMETS
mskstar = osprey.MSKStar(
	objective = osprey.MSKStar_LMFE({
		complex: 1.0,
		protein: -1.0,
		ligand: -1.0
	}),
	logFile = 'mskstar.tsv',

	# how precisely should we estimate partition function values?
	# (in range (0,1], lower is more precise)
	epsilon = 0.68,

	# how many mutations should we consider at once?
	# (keeping this low can make the search run much faster)
	maxSimultaneousMutations = 1
)

# how should we compute energies of molecules?
# (give all the conf spaces to the ecalc)
confSpaces = [state.confSpace for state in mskstar.states]
parallelism = osprey.Parallelism(cpuCores=4)
ecalc = osprey.EnergyCalculator(confSpaces, ffparams, parallelism=parallelism)

# configure MSK* states
for state in mskstar.states:

	# how should we define energies of conformations?
	eref = osprey.ReferenceEnergies(state.confSpace, ecalc)
	state.confEcalc = osprey.ConfEnergyCalculator(state.confSpace, ecalc, referenceEnergies=eref)

	# compute the energy matrix
	emat = osprey.EnergyMatrix(state.confEcalc, cacheFile='emat.%s.dat' % state.name)
	state.fragmentEnergies = emat

	# how should we compute partition functions?
	state.pfuncFactory = osprey.PartitionFunctionFactory(state.confSpace, state.confEcalc, state.name)

	# how should confs be ordered and searched? (don't forget to capture emat by using a defaulted argument)
	def makeAStar(rcs, emat=emat):
		return osprey.AStarTraditional(emat, rcs, showProgress=False)
	state.confTreeFactory = osprey.MSKStar_ConfSearchFactory(makeAStar)

# finally, run MSK*
mskstar.findBestSequences(5)
