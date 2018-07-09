
import osprey
osprey.start(heapSizeMiB=200000)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# INSERT SCRIPT OUTPUT HERE
mol = osprey.readPdb('../../pdb/1GWC.shell.pdb')
pos =  [0, 1, 2, 3, 4, 5, 6, 7, 8]
posP =  [0, 1, 2, 3, 4]
posL =  [5, 6, 7, 8]
resNumsPL =  ['C53', 'C64', 'C67', 'C71', 'C78', 'B93', 'B100', 'B101', 'B104']
resNumsP =  ['C53', 'C64', 'C67', 'C71', 'C78']
resNumsL =  ['B93', 'B100', 'B101', 'B104']
startResP =  'C7'
endResP =  'C80'
startResL =  'B89'
endResL =  'B161'
AATypeOptions =  [['HID', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'ALA', 'VAL', 'ILE', 'LEU', 'CYS', 'GLY'], ['ALA', 'LEU', 'ILE', 'VAL', 'PHE', 'TYR', 'MET', 'GLU', 'ASP','HID', 'ASN', 'GLN', 'GLY'], ['CYS', 'SER', 'THR', 'ASN', 'GLN', 'HIE', 'ALA', 'VAL', 'ILE', 'LEU', 'PHE', 'TYR', 'GLY'], ['ILE', 'ALA', 'LEU', 'VAL', 'PHE', 'TYR', 'MET','GLU', 'ASP', 'HID', 'ASN', 'GLN', 'GLY'], ['GLU', 'ASP', 'PHE', 'TYR', 'ALA', 'VAL', 'ILE', 'LEU', 'HIE', 'HID', 'ASN', 'GLN', 'GLY'], ['TYR', 'PHE', 'ALA', 'VAL', 'ILE','LEU', 'GLU', 'ASP', 'HID', 'HIE', 'ASN', 'GLN', 'GLY'], ['PHE', 'TYR', 'ALA', 'VAL', 'ILE', 'LEU', 'GLU', 'ASP', 'HID', 'HIE', 'ASN', 'GLN', 'GLY'], ['TRP', 'ALA', 'VAL','ILE', 'LEU', 'PHE', 'TYR', 'MET', 'SER', 'THR', 'ASN', 'GLN', 'GLY'], ['TYR', 'PHE', 'ALA', 'VAL', 'ILE', 'LEU', 'GLU', 'ASP', 'HID', 'HIE', 'ASN', 'GLN', 'GLY']]
# make sure all strands share the same template library (including wild-type rotamers)
templateLib = osprey.TemplateLibrary(ffparams.forcefld, moleculesForWildTypeRotamers=[mol])

# define the protein strand
ligand = osprey.Strand(mol, templateLib=templateLib, residues=[startResL, endResL])
for i in range(0,len(resNumsL)):
    ligand.flexibility[resNumsL[i]].setLibraryRotamers(*AATypeOptions[posL[i]]).addWildTypeRotamers().setContinuous()

# define the ligand strand
protein = osprey.Strand(mol, templateLib=templateLib, residues=[startResP, endResP])
for i in range(0,len(resNumsP)):
    protein.flexibility[resNumsP[i]].setLibraryRotamers(*AATypeOptions[posP[i]]).addWildTypeRotamers().setContinuous()

complex = osprey.EwakstarDoer_State('PL', osprey.ConfSpace([protein, ligand]))

# number of CPUs
cpus = 40;

# configure EWAKStar
ewakstarDoer = osprey.EwakstarDoer(

    state = complex,

    # do we want to set our baseline at the wild-type sequence? as in: do we only want sequences better than wild-type?
    useWtBenchmark = False,

    # number of sequences we want to limit ourselves to during the "sequence filter" portion of ewakstarDoer
    numEWAKStarSeqs = 10000,

        logFile = 'ewakstar.sequences.tsv',

        # how precisely should we estimate partition function values?
        # (in range (0,1], lower is more precise)
        epsilon = 0.68,

        # energy window limit for partition function calculation
        pfEw = 1.0,

        # energy window within the wild-type for finding sequences in the "sequence filter" portion of ewakstarDoer
        eW = 10.0,

        # order of magnitude worse in partition function we want to keep sequences relative to the wild-type
        orderOfMag = 5,

        # how many conformations do we want to limit our partition functions to?
        numPfConfs = 500,

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

complex.pmat = osprey.DEE(complex.confSpace, complex.emat, singlesGoldsteinDiffThreshold=10.0, pairsGoldsteinDiffThreshold=10.0, triplesGoldsteinDiffThreshold=10.0,
typeDependent=True, showProgress=False, parallelism=parallelism, cacheFile='ewakstar.%s.pmat' % complex.name)
# how should confs be ordered and searched? (don't forget to capture emat by using a defaulted argument)
def makeAStar(rcs, emat=complex.emat):
        return osprey.AStarTraditional(complex.emat, rcs, showProgress=False)
complex.confTreeFactory = osprey.EwakstarDoer_ConfSearchFactory(makeAStar)

scoredSequences = ewakstarDoer.run(complex);
