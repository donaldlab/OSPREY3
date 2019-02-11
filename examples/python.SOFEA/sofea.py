
import osprey
osprey.start()

# let's go more faster
parallelism = osprey.Parallelism(cpuCores=2)

# if you have GPUs, designs can go much faster
#parallelism = osprey.Parallelism(cpuCores=2, gpus=1, streamsPerGpu=16)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb('2RL0.min.reduce.pdb')

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# define the design strand
design = osprey.Strand(mol, templateLib=templateLib, residues=['G648', 'G654'])
design.flexibility['G649'].setLibraryRotamers(osprey.WILD_TYPE, 'ALA', 'VAL').addWildTypeRotamers().setContinuous()
design.flexibility['G650'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
design.flexibility['G651'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# define the ligand strand
ligand = osprey.Strand(mol, templateLib=templateLib, residues=['A155', 'A194'])
ligand.flexibility['A172'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand.flexibility['A192'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand.flexibility['A193'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# make a multi-state conf space
confSpace = osprey.MultiStateConfSpace([
	osprey.StateMutable('complex', osprey.ConfSpace([design, ligand])),
	osprey.StateMutable('design', osprey.ConfSpace(design)),
	osprey.StateUnmutable('ligand', osprey.ConfSpace(ligand))
])

# how should we compute energies of molecules?
ecalc = osprey.EnergyCalculator(confSpace, ffparams, parallelism=parallelism)

# make a function that makes a SOFEA config object for a state
def config(state):

	# how should we define energies of conformations?
	confEcalc = osprey.ConfEnergyCalculator(
		state.confSpace,
		ecalc,
		referenceEnergies = osprey.ReferenceEnergies(state.confSpace, ecalc),

		# use the "AllOnPairs" setting to get tighter energy bounds
		# tighter energy bounds make SOFEA run much faster!
		energyPartition = osprey.EnergyPartition.AllOnPairs
	)

	# make the SOFEA config
	return osprey.SOFEA_StateConfig(

		# calculate the energy matrix
		emat = osprey.EnergyMatrix(
			confEcalc,

			# save the emat to disk so we don't have to calculate it again later
			cacheFile = 'sofea.%s.emat' % state.name,

			# correct overly optimistic energy bounds with triples to get even tighter energy bounds!
			tripleCorrectionThreshold = 0.0
		),

		confEcalc = confEcalc
	)

# configure the SOFEA algorithm
sofea = osprey.SOFEA(
	confSpace = confSpace,
	configFunc = config,
	fringeDBPath = 'sofea.fringedb',
	fringeDBSizeMiB = 100,
	seqDBPath = 'sofea.seqdb'
)

# start new databases if this is the first time running this script
if not osprey.jvm.toFile('sofea.fringedb').exists():
	sofea.init()

# otherwise, keep the existing databases
# if you're resuming a previous design, don't call init().
# SOFEA will try to delete your databases if you do!

# choose a termination criterion for SOFEA
# ie, what are we trying to compute?
criterion = osprey.SOFEA_MinLMFE(

	# for this example, minimize a linear multi state free energy (LMFE) that encodes binding affinity
	lmfe = confSpace.lmfe()
		.addPositive('complex')
		.addNegative('design')
		.addNegative('ligand')
		.build(),

	# get the 2 sequences with the lowest LMFE
	numSequences = 2,

	# how precise do we want each free energy to be?
	minFreeEnergyWidth = 0.1,

	# free energy calculations require high-precision arithmetic
	# copy the same Boltzmann calculator used by SOFEA itself
	bcalc = sofea.bcalc
)

# finally, run the design!
sofea.refine(criterion)

# when the design is finished, write out the results
with sofea.openSeqDB() as seqdb:
	criterion.makeResultDoc(seqdb, osprey.jvm.toFile('sofea.md'))
