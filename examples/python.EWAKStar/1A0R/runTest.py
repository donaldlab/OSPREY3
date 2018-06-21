
import osprey
osprey.start()

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# INSERT SCRIPT OUTPUT HERE
mol = osprey.readPdb('../../pdb/1A0R.b.shell.pdb')
pos =  [0, 1, 2, 3, 4, 5, 6, 7, 8]
posP =  [0, 1, 2]
posL =  [3, 4, 5, 6, 7, 8]
resNumsPL =  ['0311', '0313', '0332', '0601', '0605', '0696', '0697', '0698', '0729']
resNumsP =  ['0311', '0313', '0332']
resNumsL =  ['0601', '0605', '0696', '0697', '0698', '0729']
startResP =  '039'
endResP =  '0339'
startResL =  '0520'
endResL =  '0729'
AATypeOptions =  [['HID', 'HIE', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'ALA', 'VAL', 'ILE', 'LEU', 'GLY'], ['ASN', 'SER', 'THR', 'GLN', 'HID', 'ALA', 'VAL', 'ILE', 'LEU', 'GLY', 'ALA', 'VAL', 'GLY'], ['TRP', 'ALA', 'VAL', 'ILE', 'LEU', 'PHE', 'TYR', 'MET', 'SER', 'THR', 'ASN', 'GLN', 'GLY'], ['MET', 'ILE', 'ALA', 'VAL', 'LEU', 'PHE', 'TYR', 'GLU', 'ASP', 'HID', 'ASN', 'GLN', 'GLY'], ['LEU', 'ILE', 'ALA', 'VAL', 'PHE', 'TYR', 'MET', 'GLU', 'ASP', 'HID', 'ASN', 'GLN', 'GLY'], ['GLU', 'ASP', 'PHE', 'TYR', 'ALA', 'VAL', 'ILE', 'LEU', 'HIE', 'HID', 'ASN', 'GLN', 'GLY'], ['LEU', 'ILE', 'ALA', 'VAL', 'PHE', 'TYR', 'MET', 'GLU', 'ASP', 'HID', 'ASN', 'GLN', 'GLY'], ['LEU', 'ILE', 'ALA', 'VAL', 'PHE', 'TYR', 'MET', 'GLU', 'ASP', 'HID', 'ASN', 'GLN', 'GLY'], ['GLU', 'ASP', 'PHE', 'TYR', 'ALA', 'VAL', 'ILE', 'LEU', 'HIE', 'HID', 'ASN', 'GLN', 'GLY']]

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
numCPUs = 40;
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
useMaxMutable = True
#want to change maxMutable for various designs
maxMutable = 1
numFilteredSeqs = 10000
orderOfMag = 5.0
unboundEw = 8.0
boundEw = 8.0
ewakstarEw = 1.0
Ival = 0.0

PLmatrixName = "ewak.*"

PLematMatrixName = "ewak.PL.emat"
emat = osprey.EnergyMatrix(confECalc, cacheFile=PLematMatrixName)

# run EWAK*
ewakstarDoer = osprey.EWAKStar(
        numCPUs,
        wtBenchmark,
        seqFilterOnly,
        useMaxMutable,
        maxMutable,
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

scoredSequences = ewakstarDoer.run();
