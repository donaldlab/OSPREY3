
import osprey
osprey.start()

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb('2RL0.min.reduce.pdb')

pos = [0, 1, 2, 3, 4, 5, 6, 7]
posL = [4, 5, 6, 7]
posP = [0, 1, 2, 3]

AATypeOptions = [["PHE"],["LYS"],["ILE"],["THR"],["TYR","ALA","VAL","ILE","LEU","PHE"],["ASP"],["GLU"],["THR"]]

startResL = "G648"
endResL = "G654"
startResP = "A155"
endResP = "A194"
resNumsPL = ["A156", "A172", "A192", "A193", "G649", "G650", "G651", "G654"]
resNumsL = ["G649", "G650", "G651", "G654"]
resNumsP = ["A156", "A172", "A192", "A193"]

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

# make the conf space for the protein
confSpaceP = osprey.ConfSpace(protein)

# make the conf space for the ligand
confSpaceL = osprey.ConfSpace(ligand)

# make the conf space for the protein+ligand complex
confSpace = osprey.ConfSpace([protein, ligand])

# how should we compute energies of molecules?
# (give the complex conf space to the ecalc since it knows about all the templates and degrees of freedom)
numCPUs = 4;
parallelism = osprey.Parallelism(cpuCores=numCPUs)
ecalc = osprey.EnergyCalculator(confSpace, ffparams, parallelism=parallelism)
rigidEcalc = osprey.EnergyCalculator(confSpace, ffparams, parallelism=parallelism, isMinimizing=False);

# how should we define energies of conformations?
eref = osprey.ReferenceEnergies(confSpace, ecalc)
erefRigid = osprey.ReferenceEnergies(confSpace, rigidEcalc)
confECalc = osprey.ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=eref)
confRigidECalc = osprey.ConfEnergyCalculator(confSpace, rigidEcalc, referenceEnergies=eref)

# variables for EWAKStar
epsilon = 0.68
numTopSeqs = 5
maxPFConfs = 500
seqFilterOnly = False
wtBenchmark = False
mutableType = "all"
numMutable = 1
numFilteredSeqs = 10000
orderOfMag = 5.0
unboundEw = 20.0
boundEw = 20.0
ewakstarEw = 1.0
Ival = 0.0

PLmatrixName = "ewak.*"

PLematMatrixName = "ewak.PL.emat"
emat = osprey.EnergyMatrix(confECalc, cacheFile=PLematMatrixName)

# run EWAK*
ewakstar = osprey.EWAKStar(
        numCPUs,
        wtBenchmark,
        seqFilterOnly,
        mutableType,
        numMutable,
        numTopSeqs,
        maxPFConfs,
        epsilon,
        confRigidECalc,
        confECalc,
        emat,
        ecalc,
        confSpace,
        confSpaceL,
        confSpaceP,
        pos,
        posL,
        posP,
        numFilteredSeqs,
        orderOfMag,
        unboundEw,
        boundEw,
        ewakstarEw,
        startResL,
        endResL,
        startResP,
        endResP,
        mol,
        resNumsPL,
        resNumsL,
        resNumsP,
        Ival,
        PLmatrixName
)

scoredSequences = ewakstar.run();
