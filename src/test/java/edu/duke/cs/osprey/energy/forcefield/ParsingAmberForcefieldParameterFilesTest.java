package edu.duke.cs.osprey.energy.forcefield;

import static edu.duke.cs.osprey.energy.forcefield.amber.ForcefieldFileParser.UnmatchedAtom;
import static org.assertj.core.api.Assertions.*;

import edu.duke.cs.osprey.energy.forcefield.amber.*;
import one.util.streamex.StreamEx;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.nio.file.Paths;

public class ParsingAmberForcefieldParameterFilesTest {

    String parm96DatPath = "/config/parm96.dat";
    String parm96aDatPath = "/config/parm96a.dat";
    String dabrafenibFrcmodPath = "/dabrafenib.frcmod";

    AtomSymbolAndMass hydrogenBondedToNitrogenAtoms = new AtomSymbolAndMass("H", 1.008f);
    AtomSymbolAndMass hydrogenInWater = new AtomSymbolAndMass("HW", 1.008f);
    AtomSymbolAndMass oxygenInWater = new AtomSymbolAndMass("OW", 16.00f);
    AtomSymbolAndMass sp3AliphaticCarbon = new AtomSymbolAndMass("CT", 12.01f);
    AtomSymbolAndMass sp2AromaticCarbon = new AtomSymbolAndMass("CA", 12.01f);
    AtomSymbolAndMass flourine = new AtomSymbolAndMass("F", 19.00f);
    AtomSymbolAndMass carbon = new AtomSymbolAndMass("C", 12.01f);
    AtomSymbolAndMass sp2AromaticCarbonInTryptophan = new AtomSymbolAndMass("C*", 12.01f);
    AtomSymbolAndMass oxygen = new AtomSymbolAndMass("O2", 16.00f);
    AtomSymbolAndMass carbonylOxygen = new AtomSymbolAndMass("O", 16.00f);
    AtomSymbolAndMass bromine = new AtomSymbolAndMass("BR", 79.90f);
    AtomSymbolAndMass cesium = new AtomSymbolAndMass("Cs", 132.91f);
    AtomSymbolAndMass hydrogenAliphaticBondedWithCarbon = new AtomSymbolAndMass("H1", 1.008f);
    AtomSymbolAndMass oxygenInHydroxylGroup = new AtomSymbolAndMass("OH", 16.00f);
    AtomSymbolAndMass etherAndEsterOxygen = new AtomSymbolAndMass("OS", 16.00f);
    AtomSymbolAndMass sp2Nitrogen = new AtomSymbolAndMass("N*", 14.01f);
    AtomSymbolAndMass sp2NitrogenInAmideGroups = new AtomSymbolAndMass("N", 14.01f);
    AtomSymbolAndMass sp2NitrogenInFiveMemberRing = new AtomSymbolAndMass("NB", 14.01f);
    AtomSymbolAndMass sp2CarbonPyrimidinesinPostions5And6 = new AtomSymbolAndMass("CM", 12.01f);
    AtomSymbolAndMass bigIonWithWater = new AtomSymbolAndMass("IB", 131.0f);
    AtomSymbolAndMass soAtom = new AtomSymbolAndMass("SO", 32.060f);
    AtomSymbolAndMass sp2CarbonInFiveMemberRingOfPurines = new AtomSymbolAndMass("CQ", 12.01f);
    AtomSymbolAndMass sp2NitrogenInSixMemberRing = new AtomSymbolAndMass("NC", 14.01f);

    ForcefieldFileParser readParm96() throws IOException {
        var file = new ForcefieldFileParser(getClass().getResourceAsStream(parm96DatPath));
        file.read();
        return file;
    }

    ForcefieldFileParser readParm96a() throws IOException {
        var file = new ForcefieldFileParser(getClass().getResourceAsStream(parm96aDatPath));
        file.read();
        return file;
    }

    ForcefieldFileParser readWithDabrafenibFrcmod() throws IOException {
        var file = new ForcefieldFileParser(getClass().getResourceAsStream(parm96DatPath), Paths.get(getClass().getResource(dabrafenibFrcmodPath).getFile()));
        file.read();
        return file;
    }

    @Test
    @DisplayName("Parsing vanilla parm96.dat should not throw an exception")
    void testParsingForceFieldFileThrowsNoException() {
        assertThatNoException().isThrownBy(this::readParm96);
    }

    @Test
    @DisplayName("The parsing should skip any parameters with unknown atom types")
    void testForcefieldParsingHasNoUnknownAtomTypes() throws IOException {

        var ff = readParm96();
        var atoms = ff.atomSymbolsAndMasses();
        var hydrophilicAtoms = ff.hydrophilicAtoms();
        var bondLengthAtoms = StreamEx.of(ff.bondLengthParameters())
                .flatCollection(BondLengthParameter::atoms);
        var bondAngleAtoms = StreamEx.of(ff.bondAngleParameters())
                .flatCollection(BondAngleParameter::atoms);
        var dihedralAtoms = StreamEx.of(ff.dihederalParameters())
                .flatCollection(DihederalParameter::atoms);
        var improperDihedralAtoms = StreamEx.of(ff.improperDihederalParameters())
                .flatCollection(ImproperDihederalParameter::atoms);
        var vanDerWaalsAtoms = StreamEx.of(ff.vanDerWaalsRadii())
                .flatCollection(VanDerWaalsRadius::atoms);
        var allAtoms = StreamEx.of(atoms)
                .append(hydrophilicAtoms)
                .append(bondLengthAtoms)
                .append(bondAngleAtoms)
                .append(dihedralAtoms)
                .append(improperDihedralAtoms)
                .append(vanDerWaalsAtoms)
                .toList();

        assertThat(allAtoms).doesNotContain(UnmatchedAtom);
    }

    @Test
    @DisplayName("parm96.dat contains bromine and cesium")
    void testContainsBromineAndCesium() throws IOException {
        var file = readParm96();
        assertThat(file.atomSymbolsAndMasses()).contains(bromine, cesium);
    }

    @Test
    @DisplayName("parm96.dat specifies C and O2 are hydrophilic")
    void testParsesHydrophilic() throws IOException {
        var file = readParm96();
        // NT is present in the hydrophilic atoms, but not present in the atom names

        assertThat(file.hydrophilicAtoms()).contains(carbon, oxygen);
    }

    @Test
    @DisplayName("Bond length parameters are parsed correctly from parm96.dat")
    void testParseBondLengthParameters() throws IOException {
        var file = readParm96();
        var oxygenHydrogenInWaterBond = new BondLengthParameter(oxygenInWater, hydrogenInWater, 553.0f, 0.9572f);
        var sp3CarbonFlorineBond = new BondLengthParameter(sp3AliphaticCarbon, flourine, 367.0f, 1.380f);
        assertThat(file.bondLengthParameters()).contains(oxygenHydrogenInWaterBond, sp3CarbonFlorineBond);
    }

    @Test
    @DisplayName("Bond angle parameters are parsed correctly from parm96.dat")
    void testParseBondAngleParameters() throws IOException {
        var file = readParm96();
        var hydrogenOxygenHydrogenInWaterAngle = new BondAngleParameter(hydrogenInWater, oxygenInWater, hydrogenInWater, 100.f, 104.52f);
        var flourineSp3CarbonHydrogenAngle = new BondAngleParameter(flourine, sp3AliphaticCarbon, hydrogenAliphaticBondedWithCarbon, 35.0f, 109.50f);
        assertThat(file.bondAngleParameters()).contains(hydrogenOxygenHydrogenInWaterAngle, flourineSp3CarbonHydrogenAngle);
    }

    @Test
    @DisplayName("Dihedral parameters are parsed correctly from parm96.dat")
    void testParseDihedralParameters() throws IOException {
        var file = readParm96();
        var carbonSp2CarbonGeneralDihedral = new DihederalParameter(ForcefieldFileParser.WildcardAtom, carbon, sp2AromaticCarbon, ForcefieldFileParser.WildcardAtom, 4, 14.50f, 180.0f, 2.f);
        var esterOxygenSp3CarbonSp2NitrogenSp2CarbonDihedral = new DihederalParameter(etherAndEsterOxygen, sp3AliphaticCarbon, sp2Nitrogen, sp2CarbonPyrimidinesinPostions5And6, 1, 2.50f, 0.0f, 1.f);
        assertThat(file.dihederalParameters()).contains(carbonSp2CarbonGeneralDihedral, esterOxygenSp3CarbonSp2NitrogenSp2CarbonDihedral);
    }

    @Test
    @DisplayName("Improper dihedral parameters are parsed correctly from parm96.dat")
    void testParseImproperDihedralParameters() throws IOException {
        var file = readParm96();
        var carbonCarbonylOxygenGeneralTypeDihedral = new ImproperDihederalParameter(ForcefieldFileParser.WildcardAtom, ForcefieldFileParser.WildcardAtom, carbon, carbonylOxygen, 10.5f, 180.f, 2.f);
        var sp2CarbonSp2CarbonCarbonHydroxylOxygenDihededral = new ImproperDihederalParameter(sp2AromaticCarbon, sp2AromaticCarbon, carbon, oxygenInHydroxylGroup,1.1f, 180.0f, 2.f);
        assertThat(file.improperDihederalParameters()).contains(carbonCarbonylOxygenGeneralTypeDihedral, sp2CarbonSp2CarbonCarbonHydroxylOxygenDihededral);
    }

    @Test
    @DisplayName("H-Bond 10-12 potential parameters are parsed correctly from parm96.dat")
    void testHBond10_12PotentialParameters() throws IOException {
        var file = readParm96();
        var hydrogenOxygenInWater10_12PotentialParams = new HBond10_12PotentialParameter(hydrogenInWater, oxygenInWater, 0000.f, 0000.f);
        assertThat(file.hbond10_12PotentialParameters()).contains(hydrogenOxygenInWater10_12PotentialParams);
    }

    @Test
    @DisplayName("Equivalencing atom symbols for the non-bonded 6-12 potential parameters are parsed correctly from parm96.dat")
    void testNonBonded6_12PotentialParameters() throws IOException {
        var file = readParm96();
        var a = new EquivalencingAtom(sp2NitrogenInAmideGroups, sp2NitrogenInFiveMemberRing);
        var b = new EquivalencingAtom(carbon, sp2AromaticCarbonInTryptophan);
        assertThat(file.equivalencingAtomsForNonBonded6_12PotentialParameters()).contains(a, b);
    }

    @Test
    @DisplayName("Input for 6-12 potential parameters are parsed correctly from parm96.dat")
    void testInputFor6_12PotentialParameters() throws IOException {
        var file = readParm96();
        assertThat(file.slaterKirkwoodParameters()).isNullOrEmpty();
        assertThat(file.six12PotentialCoefficients()).isNullOrEmpty();

        var hydrogenVanDerWaalsParams = new VanDerWaalsRadius(hydrogenBondedToNitrogenAtoms, 0.60000f, 0.0157f);
        var bigIonWithWaterVanDerWaalsParams = new VanDerWaalsRadius(bigIonWithWater, 5.0f, 0.1f);
        assertThat(file.vanDerWaalsRadii())
                .isNotNull()
                .isNotEmpty()
                .contains(hydrogenVanDerWaalsParams, bigIonWithWaterVanDerWaalsParams);
    }

    @Test
    @DisplayName("Extra Atom Names and Masses can be parsed from frcmod file")
    void testMassesCanBeParsedFromFrcmodFile() throws IOException {
        var file = readWithDabrafenibFrcmod();
        assertThat(file.atomSymbolsAndMasses()).contains(soAtom);
    }

    @Test
    @DisplayName("Frcmod bond lengths overwrite base parameter bond lengths")
    void testFrcmodBondLengthReplacesSameBondLength() throws IOException {
        var file = readWithDabrafenibFrcmod();
        var newBondLength = new BondLengthParameter(sp2CarbonInFiveMemberRingOfPurines, sp2NitrogenInSixMemberRing, 441.10f, 1.369f);
        var oldBondLength = new BondLengthParameter(sp2CarbonInFiveMemberRingOfPurines, sp2NitrogenInSixMemberRing, 502.0f, 1.324f);
        assertThat(file.bondLengthParameters()).contains(newBondLength);
        assertThat(file.bondLengthParameters()).doesNotContain(oldBondLength);
    }

    @Test
    @DisplayName("Frcmod file adds new bond lengths to parameters")
    void testBondLengthCanBeParsedFromFrcModFile() throws IOException {
        var file = readWithDabrafenibFrcmod();
        var bondLength = new BondLengthParameter(sp2AromaticCarbon, soAtom, 258.70f, 1.767f);
        assertThat(file.bondLengthParameters()).contains(bondLength);
    }

    @Test
    @DisplayName("Frcmod bond angles overwrite base parameter bond angles")
    void testBondAnglesFromFrcmodFileOverwriteBaseParms() throws IOException {
        var file = readWithDabrafenibFrcmod();
        var oldBondAngle = new BondAngleParameter(sp2NitrogenInSixMemberRing, sp2CarbonInFiveMemberRingOfPurines, sp2NitrogenInSixMemberRing, 70.0f, 129.10f);
        var newBondAngle = new BondAngleParameter(sp2NitrogenInSixMemberRing, sp2CarbonInFiveMemberRingOfPurines, sp2NitrogenInSixMemberRing, 69.800f, 125.700f);
        assertThat(file.bondAngleParameters()).contains(newBondAngle);
        assertThat(file.bondAngleParameters()).doesNotContain(oldBondAngle);
    }

    @Test
    @DisplayName("New bond angles from frcmod are added")
    void testAddingNewBondAnglesFromFrcmod() throws IOException {
        var file = readWithDabrafenibFrcmod();
        var newBondAngle = new BondAngleParameter(carbonylOxygen, soAtom, carbonylOxygen, 73.600f, 120.050f);
        assertThat(file.bondAngleParameters()).contains(newBondAngle);
    }

    @Test
    @DisplayName("New dihedrals can be added from frcmod")
    void testDihedralsAddedFromFrcmod() throws IOException {
        var file = readWithDabrafenibFrcmod();
        var frcmodDihedral = new DihederalParameter(sp2AromaticCarbon, sp2AromaticCarbon, soAtom, carbonylOxygen, 6, 7.800f, 180.000f, 2.000f);
        assertThat(file.dihederalParameters()).contains(frcmodDihedral);
    }

    @Test
    @DisplayName("New improper dihedrals can be added from frcmod")
    void testImproperDihedralsAddedFromFrcmod() throws IOException {
        var file = readWithDabrafenibFrcmod();
        var frcmodImproperDihedral = new ImproperDihederalParameter(sp2AromaticCarbon, sp2AromaticCarbon, sp2AromaticCarbon, soAtom, 1.1f, 180.0f, 2.0f);
        assertThat(file.improperDihederalParameters()).contains(frcmodImproperDihedral);
    }

    @Test
    @DisplayName("New non-bonded parameters can be added from frcmod")
    void testNonBondedParamsAddedFromFrcmod() throws IOException {
        var file = readWithDabrafenibFrcmod();
        var vanDerWaalsRadius = new VanDerWaalsRadius(carbonylOxygen, 1.6612f, 0.2100f);
        assertThat(file.vanDerWaalsRadii()).contains(vanDerWaalsRadius);
    }

    @Test
    @DisplayName("Can get van der Waal params from parm96a.dat")
    void testInputFor6_12PotentialParametersInParm96a() throws IOException {
        var file = readParm96a();
        assertThat(file.slaterKirkwoodParameters()).isNullOrEmpty();
        assertThat(file.six12PotentialCoefficients()).isNullOrEmpty();

        var hydrogenVanDerWaalsParams = new VanDerWaalsRadius(hydrogenBondedToNitrogenAtoms, 0.60000f, 0.0157f);
        var bigIonWithWaterVanDerWaalsParams = new VanDerWaalsRadius(bigIonWithWater, 5.0f, 0.1f);
        assertThat(file.vanDerWaalsRadii())
                .isNotNull()
                .isNotEmpty()
                .contains(hydrogenVanDerWaalsParams, bigIonWithWaterVanDerWaalsParams);
    }

    // TODO: Add HBON test for frcmod
}
