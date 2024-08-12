

'''
This module handles preparing molecules for Osprey designs,
and building Conformation Spaces from those molecules.
'''


import os, sys

# throw an error for python 2
if sys.version_info[0] <= 2:
    raise Exception('Python v2 or earlier is not supported in this module')


import osprey
import osprey.jvm as jvm


# tragically, we can't access the companion directly from JPype due to an API oversight
# but we can still get it by plain old JVM reflection
def _kotlin_companion(cls):
    return cls.class_.getDeclaredField('Companion').get(None)


class LocalService:
    '''
    Runs an instance of the Osprey Service locally on this machine.

    Only works on Linux.

    Access the service in Python code by using a `with` guard:
    ```python
    from osprey.prep import LocalService
    with LocalService():
    	# run code that uses the service here
    ```

    If a function here needs the Osprey Service, but the Osprey Service is not running,
    that function will raise an Exception.
    '''

    # make these class variables instead of instance variables
    # only only one service instance is allowed at once
    _service = None

    def __init__(self):
        # TODO: try to prevent from running on non-linux platforms?
        pass

    def __enter__(self):

        cls = LocalService

        # don't allow more than once instance at once
        if cls._service is not None:
            raise Exception('only one instance of the local service allowed at once')

        # find the progs dir
        # in the release environment, it should be near this script
        here_dir = os.path.dirname(__file__)
        progs_dir = os.path.join(here_dir, 'progs')
        if not os.path.exists(progs_dir):
            # in the dev environment, it should be at the root of the source tree
            progs_dir = os.path.join(here_dir, '../../../../progs')
        if not os.path.exists(progs_dir):
            raise Exception("can't find programs directory for osprey prep service")

        # start the service
        service_path = jvm.c.java.nio.file.Paths.get(os.path.dirname(progs_dir))
        wait = False
        port = osprey.c.service.OspreyService.defaultPort
        useVersionPrefix = True
        cls._service = osprey.c.service.OspreyService.Instance(service_path, wait, port, useVersionPrefix)

        # switch the provider to the server instance
        provider = osprey.c.gui.io.UserSettings.ServiceProvider('localhost', port, False)
        osprey.c.gui.io.OspreyService.INSTANCE.setProvider(provider)

        print('Osprey prep local service started')

    def __exit__(self, type, value, traceback):

        cls = LocalService

        if cls._service is None:
            return

        # stop the service
        cls._service.close()

        # unregister the provider
        osprey.c.gui.io.OspreyService.INSTANCE.setProvider(None)

        cls._service = None

        print('Osprey prep local service stopped')


def _check_local_service():
    if LocalService._service is None:
        raise Exception('no local service has been started yet. try `with osprey.prep.LocalService():`')


# export enums
Hybridization = osprey.c.gui.forcefield.amber.Hybridization
'''
${enum_kotlin(.gui.forcefield.amber/Hybridization)}
'''


# export statics
confLibs = osprey.c.gui.features.components.ConfLibs.INSTANCE.getInfos()
'''
${prop_kotlin(.gui.features.components/ConfLibs/infos)}
'''


class Forcefield:
    '''
    All the forcefields supported by Osprey for energy calculations

    For example, to reference one of the forcefields in Python code:
    ```python
    from osprey.prep import Forcefield
    ff = Forcefield.Amber96
    ```
    '''

    Amber96 = osprey.c.gui.forcefield.Forcefield.Amber96.INSTANCE
    '''
    ${type_kotlin(.gui.forcefield/Forcefield.Amber96)}
    ${class_kdoc(.gui.forcefield/Forcefield.Amber96)}
    '''

    Amber14SB = osprey.c.gui.forcefield.Forcefield.Amber14SB.INSTANCE
    '''
    ${type_kotlin(.gui.forcefield/Forcefield.Amber14SB)}
    ${class_kdoc(.gui.forcefield/Forcefield.Amber14SB)}
    '''

    EEF1 = osprey.c.gui.forcefield.Forcefield.EEF1.INSTANCE
    '''
    ${type_kotlin(.gui.forcefield/Forcefield.EEF1)}
    ${class_kdoc(.gui.forcefield/Forcefield.EEF1)}
    '''


def loadPDB(pdb):
    '''
    Loads molecules from a string in PDB format.
    Supports reading any kind of molecule, including proteins.

    # Arguments
    pdb `str`: The text of a PDB file

    # Returns
    [`${type_kotlin(.molscope.molecule/Molecule)}`]: The molecules from the PDB file
    '''
    pdb = osprey.c.gui.io.PDBKt.fromPDB(_kotlin_companion(osprey.c.molscope.molecule.Molecule), pdb)
    partition = osprey.c.gui.forcefield.amber.MolPartitioningKt.partition(pdb, True)
    return [pair.getSecond() for pair in partition]


def savePDB(mol):
    '''
    Saves a molecule to a string in PDB format.

    If writing the string to a file, you should use the `w` option to `open()`:
    ```python
    import osprey.prep
    with open('path/to/file.pdb', 'w') as file:
        file.write(osprey.prep.savePDB(mol))
    ```

    {{% notice tip %}}
    PDB files cannot save all the information from this molecule representation.
    Unless you explicitly need to export your molecules to PDB format,
    try saving your molecules in OMOL format, which will save all the information
    from this molecule representation.
    {{% /notice %}}

    # Arguments
    mol ${type_kotlin(.molscope.molecule/Molecule)}: The molecule

    # Returns
    str: The text of a PDB file representing the given molecule
    '''
    return osprey.c.gui.io.PDBKt.toPDB(mol, False, False, None, None, False, False, False)


def loadOMOL(omol):
    '''
    Load molecules from a string in OMOL format.

    # Arguments
    omol `str`: The text of the OMOL file.

    # Returns
    [`${type_kotlin(.molscope.molecule/Molecule)}`]: The molecules from the OMOL file
    '''
    return osprey.c.gui.io.OMOLKt.fromOMOL(_kotlin_companion(osprey.c.molscope.molecule.Molecule), omol, True)


def saveOMOL(mols):
    '''
    Saves molecules to a string in OMOL format, Osprey's native molecule file format.

    If writing the string to a file, you should use the `w` option to `open()`:
    ```python
    import osprey.prep
    with open('path/to/file.omol', 'w') as file:
        file.write(osprey.prep.saveOMOL(mols))
    ```
    The canonical filename extension for an Osprey molecule is `.omol`.

    {{% notice note %}}
    Saving Osprey molecules into OMOL format is lossless.
    All information about this molecule representation will be saved,
    unlike saving into the PDB format, which cannot represent all the information
    about this molecule representation. Prefer to save your molecules into OMOL
    format if you want to load them in Osprey again.
    {{% /notice %}}

    # Arguments
    mols `[`${type_kotlin(.molscope.molecule/Molecule)}`]`: The molecules to save

    # Returns
    str: The text of an OMOL file representing the given molecules
    '''
    return osprey.c.gui.io.OMOLKt.toOMOL(jvm.toArrayList(mols))

def loadMol2(mol2):
    mol = osprey.c.gui.io.MOL2Kt.fromMol2(_kotlin_companion(osprey.c.molscope.molecule.Molecule), mol2)
    return mol

def saveMol2(mols):
    combined = osprey.c.molscope.molecule.MoleculeKt.combine(jvm.toArrayList(mols), "", None, None).getFirst()
    return osprey.c.gui.io.MOL2Kt.toMol2(combined, None)


def loadConfSpace(toml):
    '''
    Loads a conformation space from a string

    # Arguments
    toml `str`: The text of a conformation space file, as generated by #saveConfSpace

    # Returns
    ${type_kotlin(.gui.prep/ConfSpace)}
    '''
    return osprey.c.gui.io.ConfSpaceKt.fromToml(
        _kotlin_companion(osprey.c.gui.prep.ConfSpace),
        toml
    )


def saveConfSpace(confSpace):
    '''
    Saves a conformation space to a string

    If writing the string to a file, you should use the `w` option to `open()`:
    ```python
    import osprey.prep
    with open('path/to/file.confspace', 'w') as file:
        file.write(osprey.prep.saveConfSpace(confSpace))
    ```
    The canonical filename extension for a conformation space is `.confspace`.

    # Arguments
    confSpace ${type_kotlin(.gui.prep/ConfSpace)}: The conformation space to save

    # Returns
    str: The text of a conformation space file representing the given conformation space
    '''
    return osprey.c.gui.io.ConfSpaceKt.toToml(confSpace)


def saveCompiledConfSpace(ccs):
    '''
    Saves a compiled conformation space to bytes

    If writing the bytes to a file, you should use the `wb` option to `open()`:
    ```python
    import osprey.prep
    with open('path/to/file.ccsx', 'wb') as file:
        file.write(osprey.prep.saveCompiledConfSpace(ccs))
    ```
    The canonical filename extension for a compiled conformation space with compression is `.ccsx`.

    {{% notice note %}}
    Due to the large size of the compiled conformation space,
    the saved representation is an efficient compressed binary format,
    rather than a human-readable text format.
    {{% /notice %}}

    # Arguments
    ccs ${type_kotlin(.gui.compiler/CompiledConfSpace)}: The compiled conformation space

    # Returns
    `${returns_method_java(.tools.LZMA2#compress(byte[])byte[])}`: The bytes of a compiled conformation space file representing the given compiled conformation space
    '''
    bytes = osprey.c.gui.io.CompiledConfSpaceKt.toBytes(ccs)
    return osprey.c.tools.LZMA2.compress(bytes)


def molTypes(mol):
    '''
    ${func_kdoc(.gui.forcefield.amber//findTypes)}

    # Arguments
    ${receiver_kotlin(mol, .gui.forcefield.amber//findTypes)}:
    The molecule (or multiple molecules in a single Molecule instance) for which to determine molecule types

    # Returns
    ${returns_func_kotlin(.gui.forcefield.amber//findTypes)}
    '''
    return osprey.c.gui.forcefield.amber.MolPartitioningKt.findTypes(mol)


def duplicateAtoms(mol):
    '''
    ${class_kdoc(.gui.prep/DuplicateAtoms)}

    # Arguments
    ${arg_kotlin(mol, .gui.prep/DuplicateAtoms/DuplicateAtoms)}:
    The molecule for which to search for duplicate atoms

    # Returns
    ${returns_func_kotlin(.gui.prep/DuplicateAtoms/DuplicateAtoms)}
    '''
    return osprey.c.gui.prep.DuplicateAtoms(mol)


def inferMissingAtoms(mol):
    '''
    ${func_kdoc(.gui.forcefield.amber//inferMissingAtomsAmber)}

    # Arguments
    ${receiver_kotlin(mol, .gui.forcefield.amber//inferMissingAtomsAmber)}:
    The molecule for which to infer missing atoms

    # Returns
    ${returns_func_kotlin(.gui.forcefield.amber//inferMissingAtomsAmber)}
    '''
    _check_local_service()
    return osprey.c.gui.forcefield.amber.MissingAtomsKt.inferMissingAtomsAmberBlocking(mol)


def inferBonds(mol):
    '''
    ${func_kdoc(.gui.forcefield.amber//inferBondsAmber)}

    # Arguments
    ${receiver_kotlin(mol, .gui.forcefield.amber//inferBondsAmber)}:
    The molecule for which to infer bonds

    # Returns
    ${returns_func_kotlin(.gui.forcefield.amber//inferBondsAmber)}
    '''
    _check_local_service()
    return osprey.c.gui.forcefield.amber.BondsKt.inferBondsAmberBlocking(mol)


def inferProtonation(mol):
    '''
    ${func_kdoc(.gui.forcefield.amber//inferProtonation)}

    # Arguments
    ${receiver_kotlin(mol, .gui.forcefield.amber//inferProtonation)}:
    The molecule for which to infer protonation

    # Returns
    ${returns_func_kotlin(.gui.forcefield.amber//inferProtonation)}
    '''
    _check_local_service()
    return osprey.c.gui.forcefield.amber.ProtonateKt.inferProtonationBlocking(mol)


def deprotonate(mol, atoms=None):
    '''
    Deprotonates a molecule or certain atoms in a molecule.

    Without supplying an argument for atom, the entire molecule is deprotonated.
    When atom is supplied, it should be supplied as an iterable of atoms.

    # Arguments
    ${receiver_kotlin(mol, .gui.forcefield.amber//deprotonate/.molscope.molecule.Molecule#)}:
    The molecule from which to remove atoms

    atoms `[`${type_kotlin(.molscope.molecule/Atom)}`]`:
    An iterable of atoms from which to remove protons, instead of deprotonating the whole molecule.
    '''
    if atoms is None:
        osprey.c.gui.forcefield.amber.ProtonateKt.deprotonate(mol)
        return

    for atom in atoms:
        osprey.c.gui.forcefield.amber.ProtonateKt.deprotonate(mol, atom)


def protonate(mol, atom, numH, hybridization):
    '''
    ${func_kdoc(.gui.forcefield.amber//protonate)}

    # Arguments
    ${receiver_kotlin(mol, .gui.forcefield.amber//protonate)}: The molecule to which to add Hydrogen atoms
    ${arg_kotlin(atom, .gui.forcefield.amber//findProtonation)}: The atom in the molecule to which to add Hydrogen atoms
    ${arg_kotlin(numH, .gui.forcefield.amber//findProtonation)}: The number of Hydrogen atoms to add
    ${arg_kotlin(hybridization, .gui.forcefield.amber//findProtonation)}: The hybridization state
    '''
    _check_local_service()
    protonation = osprey.c.gui.forcefield.amber.ProtonateKt.findProtonation(mol, atom, numH, hybridization)
    if protonation is None:
        raise Exception('H=%d %s protonation not found for atom %s' % (numH, hybridization, atom))
    osprey.c.gui.forcefield.amber.ProtonateKt.protonateBlocking(mol, atom, protonation)


def minimizerInfo(mol, restrainedAtoms):
    '''
    Holds molecule information for a later call to #minimize

    # Arguments
    ${arg_kotlin(mol, .gui.forcefield.amber/MinimizerInfo/MinimizerInfo)}
    restrainedAtoms `[`${type_kotlin(.molscope.molecule/Atom)}`]`: List of atoms to restrain during the minimization

    # Returns
    ${type_kotlin(.gui.forcefield.amber/MinimizerInfo)}
    '''
    return osprey.c.gui.forcefield.amber.MinimizerInfo(mol, jvm.toArrayList(restrainedAtoms))


def minimize(minimizerInfos, numSteps):
    '''
    ${func_kdoc(.gui.forcefield.amber//minimize)}

    # Arguments
    minimizerInfos `[`${type_kotlin(.gui.forcefield.amber/MinimizerInfo)}`]`: The molecules (and their associated information) to minimize
    ${arg_kotlin(numSteps, .gui.forcefield.amber//minimize)}: The number of steps to minimize
    '''
    _check_local_service()
    osprey.c.gui.forcefield.amber.MinimizationKt.minimizeBlocking(jvm.toArrayList(minimizerInfos), numSteps)


def ConfSpace(mols):
    '''
    Create a new conformation space from a list of molecules

    # Arguments
    ${arg_kotlin(mols, .gui.prep/ConfSpace.Companion/fromMols)}

    # Returns
    ${returns_func_kotlin(.gui.prep/ConfSpace.Companion/fromMols)}
    '''
    return _kotlin_companion(osprey.c.gui.prep.ConfSpace).fromMols(mols)


def ProteinDesignPosition(protein, chainId, resId, name=None):
    '''
    ${func_kdoc(.gui.prep/Proteins/makeDesignPosition)}

    # Arguments
    protein ${type_kotlin(.molscope.molecule/Polymer)}: The protein from which to make the design position
    ${args_kotlin(.molscope.molecule/Polymer/findResidueOrThrow/#kotlin.String#kotlin.String,
        [chainId],
        [resId]
    )}
    ${arg_kotlin(name, .gui.prep/Proteins/makeDesignPosition)}: a recognizable name for the design position

    # Returns
    ${returns_func_kotlin(.gui.prep/Proteins/makeDesignPosition)}
    '''
    res = protein.findResidueOrThrow(chainId, resId)
    if name is None:
        name = '%s %s' % (resId, res.getType())
    return osprey.c.gui.prep.Proteins.INSTANCE.makeDesignPosition(protein, res, name)


def DihedralAngleSettings(radiusDegrees = 9.0, includeHydroxyls = True, includeNonHydroxylHGroups = False):
    '''
    Settings for dihedral angles

    The default values match Osprey's classic settings.

    # Arguments
    radiusDegrees `float`: Amount of rotation allowed in either direction from the starting angle, in degrees
    includeHydroxyls `bool`: True to include dihedral angles for hydroxyl groups
    includeNonHydroxylHGroups `bool`: True to include dihedral angles for other groups, like methyl groups

    # Returns
    ${type_kotlin(.gui.motions/DihedralAngle.LibrarySettings)}
    '''
    return osprey.c.gui.motions.DihedralAngle.LibrarySettings(
    	radiusDegrees,
    	includeHydroxyls,
    	includeNonHydroxylHGroups
    )


def conformationDihedralAngles(pos, confInfo, settings):
    '''
    Create motions from all the dihedral angles in the conformation library,
    so the motions can be added to a conformation at a design position

    # Arguments
    ${arg_kotlin(pos, .gui.motions/DihedralAngle.ConfDescription.Companion/makeFromLibrary)}: the design position
    confInfo ${type_kotlin(.gui.prep/ConfSpace.ConfConfSpace)}: the information of the conformation to which we're about to add the dihedral angles,
    from ${link_func_kotlin(.gui.prep/ConfSpace/getConformations/#edu.duke.cs.osprey.gui.prep.DesignPosition#kotlin.String)}
    ${arg_kotlin(settings, .gui.motions/DihedralAngle.ConfDescription.Companion/makeFromLibrary)}: the dihedral angle settings,
    from #DihedralAngleSettings

    # Returns
    ${returns_func_kotlin(.gui.motions/DihedralAngle.ConfDescription.Companion/makeFromLibrary)}
    '''
    return _kotlin_companion(osprey.c.gui.motions.DihedralAngle.ConfDescription)\
        .makeFromLibrary(pos, confInfo.getFrag(), confInfo.getConf(), settings)


def moleculeDihedralAngle(mol, a, b, c, d, settings):
    '''
    Creates a dihedral angle motion for a general molecule, outside of a design position

    # Arguments
    mol ${type_kotlin(.molscope.molecule/Molecule)}: The molecule for which to create the motion
    a `str`: the name of first atom to define the dihedral angle
    b `str`: the name of second atom to define the dihedral angle
    c `str`: the name of third atom to define the dihedral angle
    d `str`: the name of fourth atom to define the dihedral angle
    settings ${type_kotlin(.gui.motions/DihedralAngle.LibrarySettings)}: dihedral angle settings, from #DihedralAngleSettings

    # Returns
    ${type_kotlin(.gui.motions/DihedralAngle.MolDescription)}
    '''
    return osprey.c.gui.motions.DihedralAngle.MolDescription(
        mol,
        mol.getAtoms().findOrThrow(a),
        mol.getAtoms().findOrThrow(b),
        mol.getAtoms().findOrThrow(c),
        mol.getAtoms().findOrThrow(d),
        settings.getRadiusDegrees()
    )


def moleculeTranslationRotation(mol, dist=1.2, degrees=5.0):
    '''
    Creates a translation/rotation motion for a general molecule

    # Arguments
    mol ${type_kotlin(.molscope.molecule/Molecule)}: The molecule for which to create the motion
    dist `float`: the maximum distance in Angstroms to translate from the original position
    degrees `float`: the maximum angle in degrees to rotate from the original orientation

    # Returns
    ${type_kotlin(.gui.motions/TranslationRotation.MolDescription)}
    '''
    return osprey.c.gui.motions.TranslationRotation.MolDescription(
        mol,
        dist,
        degrees
    )


def ConfSpaceCompiler(confSpace):
    '''
    Creates the compiler for a conformation space

    # Arguments
    ${arg_kotlin(confSpace, .gui.compiler/ConfSpaceCompiler/ConfSpaceCompiler)}: the conformation space to compile

    # Returns
    ${type_kotlin(.gui.compiler/ConfSpaceCompiler)}
    '''
    return osprey.c.gui.compiler.ConfSpaceCompiler(confSpace)
