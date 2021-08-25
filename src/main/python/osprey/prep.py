
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

# export statics
confLibs = osprey.c.gui.features.components.ConfLibs.INSTANCE.getInfos()

class Forcefield:
    Amber96 = osprey.c.gui.forcefield.Forcefield.Amber96.INSTANCE
    Amber14SB = osprey.c.gui.forcefield.Forcefield.Amber14SB.INSTANCE
    EEF1 = osprey.c.gui.forcefield.Forcefield.EEF1.INSTANCE


def loadPDB(pdb):
    pdb = osprey.c.gui.io.PDBKt.fromPDB(_kotlin_companion(osprey.c.molscope.molecule.Molecule), pdb)
    partition = osprey.c.gui.forcefield.amber.MolPartitioningKt.partition(pdb, True)
    return [pair.getSecond() for pair in partition]


def loadOMOL(omol):
    return osprey.c.gui.io.OMOLKt.fromOMOL(_kotlin_companion(osprey.c.molscope.molecule.Molecule), omol, True)


def saveOMOL(mols):
    return osprey.c.gui.io.OMOLKt.toOMOL(jvm.toArrayList(mols))


def loadConfSpace(toml):
    return osprey.c.gui.io.ConfSpaceKt.fromToml(
        _kotlin_companion(osprey.c.gui.prep.ConfSpace),
        toml
    )


def saveConfSpace(confSpace):
    return osprey.c.gui.io.ConfSpaceKt.toToml(confSpace)


def saveCompiledConfSpace(ccs):
    bytes = osprey.c.gui.io.CompiledConfSpaceKt.toBytes(ccs)
    return osprey.c.tools.LZMA2.compress(bytes)


def molTypes(mol):
    return osprey.c.gui.forcefield.amber.MolPartitioningKt.findTypes(mol)


def duplicateAtoms(mol):
    return osprey.c.gui.prep.DuplicateAtoms(mol)


def inferMissingAtoms(mol):
    _check_local_service()
    return osprey.c.gui.forcefield.amber.MissingAtomsKt.inferMissingAtomsAmberBlocking(mol)


def inferBonds(mol):
    _check_local_service()
    return osprey.c.gui.forcefield.amber.BondsKt.inferBondsAmberBlocking(mol)


def inferProtonation(mol):
    _check_local_service()
    return osprey.c.gui.forcefield.amber.ProtonateKt.inferProtonationBlocking(mol)


def deprotonate(mol, atom):
    osprey.c.gui.forcefield.amber.ProtonateKt.deprotonate(mol, atom)


def protonate(mol, atom, numH, hybridization):
    _check_local_service()
    protonation = osprey.c.gui.forcefield.amber.ProtonateKt.findProtonation(mol, atom, numH, hybridization)
    if protonation is None:
        raise Exception('H=%d %s protonation not found for atom %s' % (numH, hybridization, atom))
    osprey.c.gui.forcefield.amber.ProtonateKt.protonateBlocking(mol, atom, protonation)


def minimizerInfo(mol, restrainedAtoms):
    return osprey.c.gui.forcefield.amber.MinimizerInfo(mol, jvm.toArrayList(restrainedAtoms))


def minimize(minimizerInfos, numSteps):
    _check_local_service()
    osprey.c.gui.forcefield.amber.MinimizationKt.minimizeBlocking(jvm.toArrayList(minimizerInfos), numSteps)


def ConfSpace(mols):
    return _kotlin_companion(osprey.c.gui.prep.ConfSpace).fromMols(mols)


def ProteinDesignPosition(protein, chainId, resId, name=None):
    res = protein.findResidueOrThrow(chainId, resId)
    if name is None:
        name = '%s %s' % (resId, res.getType())
    return osprey.c.gui.prep.Proteins.INSTANCE.makeDesignPosition(protein, res, name)


def DihedralAngleSettings():
    return osprey.c.gui.motions.DihedralAngle.LibrarySettings(
        # set some reasonable defaults
        9.0, # radiusDegrees
        True, # includeHydroxyls
        False # includeNonHydroxylHGroups
    )


def conformationDihedralAngles(pos, confInfo, settings):
    return _kotlin_companion(osprey.c.gui.motions.DihedralAngle.ConfDescription)\
        .makeFromLibrary(pos, confInfo.getFrag(), confInfo.getConf(), settings)


def moleculeDihedralAngle(mol, a, b, c, d, settings):
    return osprey.c.gui.motions.DihedralAngle.MolDescription(
        mol,
        mol.getAtoms().findOrThrow(a),
        mol.getAtoms().findOrThrow(b),
        mol.getAtoms().findOrThrow(c),
        mol.getAtoms().findOrThrow(d),
        settings.getRadiusDegrees()
    )


def moleculeTranslationRotation(mol, dist=1.2, degrees=5.0):
    return osprey.c.gui.motions.TranslationRotation.MolDescription(
        mol,
        dist,
        degrees
    )


def ConfSpaceCompiler(confSpace):
    return osprey.c.gui.compiler.ConfSpaceCompiler(confSpace)
