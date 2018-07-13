
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

# configure MSK*
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
	maxSimultaneousMutations = 1,

	# Let MSK* save memory by guaranteeing only a few conf trees will be kept in memory at once.
	# Extra conf trees get deleted when we need more memory, but they can be recreated when needed again.
	# (this setting works in tandem with the maxNumNodes setting down below)
	minNumConfTrees = 5
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

	# how should confs be ordered and searched? (don't forget to capture emat by using a defaulted argument)
	def makeAStar(rcs, emat=emat):
		return osprey.AStarTraditional(
			emat,
			rcs,
			showProgress = False,

			# Set a limit on the size of the conf tree so we don't use unlimited memory.
			# The minimum possible value here is (num design positions + 1), but setting this too
			# small can make A* search take much much longer. So pick something big enough to keep
			# A* running efficiently, but not so big that all the nodes don't fit in memory.
			# Keep in mind, COMETS will keep multiple conf trees in memory at the same time.
			# Nodes use a variable amount of memory though, so there's no hard correspondence
			# between numbers of node and amount of memory used. You might need to experiment
			# to see what values work best here
			# (this setting works in tandem with the minNumConfTrees setting above)
			maxNumNodes = 1000000
		)
	state.confTreeFactory = osprey.MSKStar_ConfSearchFactory(makeAStar)

# finally, run MSK*
mskstar.findBestSequences(5)
