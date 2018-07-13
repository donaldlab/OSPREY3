
import osprey
osprey.start()

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb('2RL0.min.reduce.pdb')

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# define the protein strand
protein = osprey.Strand(mol, templateLib=templateLib, residues=['G648', 'G654'])
protein.flexibility['G649'].setLibraryRotamers(osprey.WILD_TYPE, 'TYR', 'ALA', 'VAL', 'ILE', 'LEU').addWildTypeRotamers().setContinuous()
protein.flexibility['G650'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['G651'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['G654'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# define the ligand strand
ligand = osprey.Strand(mol, templateLib=templateLib, residues=['A155', 'A194'])
ligand.flexibility['A156'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand.flexibility['A172'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand.flexibility['A192'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand.flexibility['A193'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# define the COMETS states
bound = osprey.COMETS_State('Bound', osprey.ConfSpace([protein, ligand]))
unbound = osprey.COMETS_State('Unbound', osprey.ConfSpace(protein))

# configure COMETS
comets = osprey.COMETS(
	objective = osprey.COMETS_LME({
		bound: 1.0,
		unbound: -1.0
	}),
	logFile = 'comets.tsv'
)

# how should we compute energies of molecules?
# (give all the conf spaces to the ecalc)
confSpaces = [state.confSpace for state in comets.states]
parallelism = osprey.Parallelism(cpuCores=4)
ecalc = osprey.EnergyCalculator(confSpaces, ffparams, parallelism=parallelism)

# configure COMETS states
for state in comets.states:

	# how should we define energies of conformations?
	eref = osprey.ReferenceEnergies(state.confSpace, ecalc)
	state.confEcalc = osprey.ConfEnergyCalculator(state.confSpace, ecalc, referenceEnergies=eref)

	# compute the energy matrix
	emat = osprey.EnergyMatrix(state.confEcalc, cacheFile='emat.%s.dat' % state.name)
	state.fragmentEnergies = emat

	# how should confs be ordered and searched? (don't forget to capture emat by using a defaulted argument)
	def makeAStar(rcs, emat=emat):
		return osprey.AStarTraditional(emat, rcs, showProgress=False)
	state.confTreeFactory = osprey.COMETS_ConfSearchFactory(makeAStar)

# finally, run COMETS
comets.findBestSequences(5)
