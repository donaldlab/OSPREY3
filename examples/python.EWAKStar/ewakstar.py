
import osprey
osprey.start()

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb('2RL0.min.reduce.pdb')

# make sure all strands share the same template library (including wild-type rotamers)
templateLib = osprey.TemplateLibrary(ffparams.forcefld, moleculesForWildTypeRotamers=[mol])

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

complex = osprey.EwakstarDoer_State('PL', osprey.ConfSpace([proteinStrand, ligandStrand]))

# number of CPUs
cpus = 4;

# configure EWAKStar
ewakstarDoer = osprey.EwakstarDoer(

    state = complex,

    smaNodes = 0,

    useSMA = False,

    printPDBs = True,

    # do we want to set our baseline at the wild-type sequence? as in: do we only want sequences better than wild-type?
    useWtBenchmark = False,

	logFile = 'ewakstar.sequences.tsv',

	# how precisely should we estimate partition function values?
	# (in range (0,1], lower is more precise)
	epsilon = 0.01,

	# energy window limit for partition function calculation
	pfEw = 1.0,

	# energy window within the wild-type for finding sequences in the "sequence filter" portion of ewakstarDoer
	eW = 30.0,

	# order of magnitude worse in partition function we want to keep sequences relative to the wild-type
	orderOfMag = 10,

	# how many conformations do we want to limit our partition functions to?
	numPfConfs = 5000,

	# how many top sequences do we want K* estimates for in the end?
	numTopSeqs = 5,

	# what is the mutable type? exact, all, or max
	mutableType = 'exact',

	# how many mutable?
	numMutable = 1,

	# do you only want to do the sequence filter portion?
	seqFilterOnly = False,

	# how many CPUs paralellism?
	numCPUs = cpus

)

# how should we compute energies of molecules?
# (give the complex conf space to the ecalc since it knows about all the templates and degrees of freedom)
parallelism = osprey.Parallelism(cpuCores=cpus)
ecalc = osprey.EnergyCalculator(complex.confSpace, ffparams, parallelism=parallelism)
rigidEcalc = osprey.EnergyCalculator(complex.confSpace, ffparams, parallelism=parallelism, isMinimizing=False);

# how should we define energies of conformations?
eref = osprey.ReferenceEnergies(complex.confSpace, ecalc)
erefRigid = osprey.ReferenceEnergies(complex.confSpace, rigidEcalc)

# how should we define energies of conformations?
eref = osprey.ReferenceEnergies(complex.confSpace, ecalc)
complex.confEcalc = osprey.ConfEnergyCalculator(complex.confSpace, ecalc, referenceEnergies=eref)
complex.confRigidEcalc = osprey.ConfEnergyCalculator(complex.confSpace, rigidEcalc, referenceEnergies=erefRigid)

# compute the energy matrix
complex.emat = osprey.EnergyMatrix(complex.confEcalc, cacheFile='ewakstar.%s.emat' % complex.name)
complex.fragmentEnergies = complex.emat
complex.ematRigid = osprey.EnergyMatrix(complex.confRigidEcalc, cacheFile='ewakstar.%s.ematRigid' % complex.name)


# how should confs be ordered and searched? (don't forget to capture emat by using a defaulted argument)
useSMA = False #do you want to use memory bounded A*?
if(useSMA):
    def makeAStar(rcs, emat=complex.emat):
	    return osprey.AStarTraditional(complex.emat, rcs, showProgress=False, maxNumNodes=2000000)
    def makeAStarRigid(rcs, emat=complex.ematRigid):
	    return osprey.AStarTraditional(complex.ematRigid, rcs, showProgress=False, maxNumNodes=2000000)
else:
    def makeAStar(rcs, emat=complex.emat):
    	return osprey.AStarTraditional(complex.emat, rcs, showProgress=False)
    def makeAStarRigid(rcs, emat=complex.ematRigid):
    	return osprey.AStarTraditional(complex.ematRigid, rcs, showProgress=False)

complex.confTreeFactoryMin = osprey.EwakstarDoer_ConfSearchFactory(makeAStar)
complex.confTreeFactoryRigid = osprey.EwakstarDoer_ConfSearchFactory(makeAStarRigid)

scoredSequences = ewakstarDoer.run(complex);
