package edu.duke.cs.osprey.confspace.compiled;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class TestConfSpace {

	public static class AffinityClassic {

		public final SimpleConfSpace complex;
		public final SimpleConfSpace chainA;
		public final SimpleConfSpace chainB;
		public final ForcefieldParams ffparams;

		public AffinityClassic(SimpleConfSpace complex, SimpleConfSpace chainA, SimpleConfSpace chainB, ForcefieldParams ffparams) {
			this.complex = complex;
			this.chainA = chainA;
			this.chainB = chainB;
			this.ffparams = ffparams;
		}

		public int[] makeConfComplexWt() {

			SimpleConfSpace confSpace = complex;

			int[] conf = Conf.make(confSpace);
			for (int posi=0; posi<confSpace.numPos(); posi++) {
				SimpleConfSpace.ResidueConf resConf = confSpace.positions.get(posi).resConfs.stream()
					.filter(rc -> rc.type == SimpleConfSpace.ResidueConf.Type.WildType)
					.findFirst()
					.orElseThrow();
				conf[posi] = resConf.index;
			}

			return conf;
		}
	}

	public static class AffinityCompiled {

		public final ConfSpace complex;
		public final ConfSpace chainA;
		public final ConfSpace chainB;

		public AffinityCompiled(ConfSpace complex, ConfSpace chainA, ConfSpace chainB) {
			this.complex = complex;
			this.chainA = chainA;
			this.chainB = chainB;
		}

		public int[] makeConfComplexWt() {

			ConfSpace confSpace = complex;

			int[] conf = Conf.make(confSpace);
			for (int posi=0; posi<confSpace.numPos(); posi++) {
				int fposi = posi;
				conf[posi] = IntStream.range(0, confSpace.numConf(posi))
					.filter(confi -> confSpace.confId(fposi, confi).startsWith("wt-"))
					.findFirst()
					.orElseThrow();
			}

			return conf;
		}
	}

	public static class Design2RL0Interface7Mut {

		public static AffinityClassic makeClassic() {

			// configure the forcefield
			ForcefieldParams ffparams = new ForcefieldParams();

			Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

			// make sure all strands share the same template library
			ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
				.clearTemplateCoords()
				.addTemplateCoords(FileTools.readResource("/config/template_coords_v2.txt"))
				// add a special template so we don't delete the un-protonated terminal residues
				// (they're missing the amide protons for some reason ...)
				.addTemplates(
					"\n" +
					"\n" +
					"GLYCINE un-capped                                               \n" +
					"                                                                \n" +
					" GLY  INT     1                                                 \n" +
					" CORR OMIT DU   BEG                                             \n" +
					"   0.00000                                                      \n" +
					"   1  DUMM  DU    M    0  -1  -2     0.000     0.000     0.000   0.00000\n" +
					"   2  DUMM  DU    M    1   0  -1     1.449     0.000     0.000   0.00000\n" +
					"   3  DUMM  DU    M    2   1   0     1.522   111.100     0.000   0.00000\n" +
					"   4  N     N     M    3   2   1     1.335   116.600   180.000  -0.41570\n" +
					"   5  CA    CT    M    4   3   2     1.449   121.900   180.000  -0.02520\n" +
					"   6  HA2   H1    E    5   4   3     1.090   109.500   300.000   0.06980\n" +
					"   7  HA3   H1    E    5   4   3     1.090   109.500    60.000   0.06980\n" +
					"   8  C     C     M    5   4   3     1.522   110.400   180.000   0.59730\n" +
					"   9  O     O     E    8   5   4     1.229   120.500     0.000  -0.56790\n" +
					"                                                                \n" +
					"IMPROPER                                                        \n" +
					" CA   +M   C    O                                               \n" +
					"                                                                \n" +
					"DONE\n" +
					"GLUTAMIC ACID un-capped                                         \n" +
					"                                                                \n" +
					" GLU  INT     1                                                 \n" +
					" CORR OMIT DU   BEG                                             \n" +
					"   0.00000                                                      \n" +
					"   1  DUMM  DU    M    0  -1  -2     0.000     0.000     0.000   0.00000\n" +
					"   2  DUMM  DU    M    1   0  -1     1.449     0.000     0.000   0.00000\n" +
					"   3  DUMM  DU    M    2   1   0     1.522   111.100     0.000   0.00000\n" +
					"   4  N     N     M    3   2   1     1.335   116.600   180.000  -0.51630\n" +
					"   5  CA    CT    M    4   3   2     1.449   121.900   180.000   0.03970\n" +
					"   6  HA    H1    E    5   4   3     1.090   109.500   300.000   0.11050\n" +
					"   7  CB    CT    3    5   4   3     1.525   111.100    60.000   0.05600\n" +
					"   8  HB2   HC    E    7   5   4     1.090   109.500   300.000  -0.01730\n" +
					"   9  HB3   HC    E    7   5   4     1.090   109.500    60.000  -0.01730\n" +
					"  10  CG    CT    3    7   5   4     1.510   109.470   180.000   0.01360\n" +
					"  11  HG2   HC    E   10   7   5     1.090   109.500   300.000  -0.04250\n" +
					"  12  HG3   HC    E   10   7   5     1.090   109.500    60.000  -0.04250\n" +
					"  13  CD    C     B   10   7   5     1.527   109.470   180.000   0.80540\n" +
					"  14  OE1   O2    E   13  10   7     1.260   117.200    90.000  -0.81880\n" +
					"  15  OE2   O2    E   13  10   7     1.260   117.200   270.000  -0.81880\n" +
					"  16  C     C     M    5   4   3     1.522   111.100   180.000   0.53660\n" +
					"  17  O     O     E   16   5   4     1.229   120.500     0.000  -0.58190\n" +
					"                                                                \n" +
					"IMPROPER                                                        \n" +
					" CA   +M   C    O                                               \n" +
					" CG   OE1  CD   OE2                                             \n" +
					"                                                                \n" +
					"DONE\n"
				)
				.build();

			Strand chainG = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("G638", "G654")
				.build();
			chainG.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
			chainG.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
			chainG.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "ASP").addWildTypeRotamers().setContinuous();

			Strand chainA = new Strand.Builder(mol)
				.setTemplateLibrary(templateLib)
				.setResidues("A153", "A241")
				.build();
			chainA.flexibility.get("A156").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
			chainA.flexibility.get("A172").setLibraryRotamers(Strand.WildType, "ASP", "GLU").addWildTypeRotamers().setContinuous();
			chainA.flexibility.get("A192").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "PHE", "TYR").addWildTypeRotamers().setContinuous();
			chainA.flexibility.get("A193").setLibraryRotamers(Strand.WildType, "SER", "ASN").addWildTypeRotamers().setContinuous();

			return new AffinityClassic(
				new SimpleConfSpace.Builder()
					.addStrands(chainG, chainA)
					.build(),
				new SimpleConfSpace.Builder()
					.addStrand(chainA)
					.build(),
				new SimpleConfSpace.Builder()
					.addStrand(chainG)
					.build(),
				ffparams
			);
		}

		public static AffinityCompiled makeCompiled() {
			return new AffinityCompiled(
				ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.complex.ccsx")),
				ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.A.ccsx")),
				ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.G.ccsx"))
			);
		}
	}

	@Test
	public void check2RL0_compiled() {

		AffinityCompiled design = Design2RL0Interface7Mut.makeCompiled();

		assert2RL0(design.chainB, design.chainA, design.complex);

		// make a molecule of the wild-type complex
		ConfSpace confSpace = design.complex;
		int[] conf = design.makeConfComplexWt();
		AssignedCoords coords = new AssignedCoords(confSpace, conf);
		Molecule mol = coords.toMol();

		assert2RL0WildType(mol);

		// check the rigid energy
		ConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(confSpace);

		List<PosInter> inters = PosInterDist.all(confSpace, null, conf);
		double energy = ecalc.calcEnergy(conf, inters);

		// This energy is off from the classic OSPREY's energy by about 8.5 kcal/mol.
		// As far as I can tell, the difference in energy is entirely due to differences in
		// electrostatic charges (with premultiplied Coulomb factors)
		// between osprey classic and AmberTools19.
		// So this energy value is as correct as we're going to get.
		assertThat(energy, isAbsolutely(-1556.9551045257604, 1e-9));
	}

	@Test
	public void check2RL0_classic() {

		// make a version of the 2RL0 "design" that matches the conf spaces prepped by the GUI

		AffinityClassic design = Design2RL0Interface7Mut.makeClassic();


		assert2RL0(design.chainB, design.chainA, design.complex);

		// make a molecule of the wild-type complex
		SimpleConfSpace confSpace = design.complex;
		ParametricMolecule pmol = confSpace.makeMolecule(design.makeConfComplexWt());

		assert2RL0WildType(pmol.mol);

		// check the rigid energies
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setIsMinimizing(false)
			.build()) {

			ResidueInteractions inters = new ResidueInteractions();
			inters.addComplete(pmol.mol.residues);
			double energy = ecalc.calcEnergy(pmol, inters).energy;

			assertThat(energy, isAbsolutely(-1548.4851738174825, 1e-6));
		}
	}

	private static void assert2RL0(ConfSpaceIteration chainG, ConfSpaceIteration chainA, ConfSpaceIteration complex) {

		Consumer<ConfSpaceIteration> assertG = (confSpace) -> {
			assertPosition(confSpace, Arrays.asList("G649", "649 PHE"), "PHE", 29, Arrays.asList("PHE:5",  "TYR:8",  "ALA:1", "VAL:3", "ILE:7", "LEU:5"));
			assertPosition(confSpace, Arrays.asList("G650", "650 ASP"), "ASP", 14, Arrays.asList("ASP:6",  "GLU:8"));
			assertPosition(confSpace, Arrays.asList("G651", "651 GLU"), "GLU", 14, Arrays.asList("GLU:9",  "ASP:5"));
		};

		Consumer<ConfSpaceIteration> assertA = (confSpace) -> {
			assertPosition(confSpace, Arrays.asList("A156", "156 PHE"), "PHE", 29, Arrays.asList("PHE:5",  "TYR:8",  "ALA:1", "VAL:3", "ILE:7", "LEU:5"));
			assertPosition(confSpace, Arrays.asList("A172", "172 LYS"), "LYS", 41, Arrays.asList("LYS:28", "ASP:5",  "GLU:8"));
			assertPosition(confSpace, Arrays.asList("A192", "192 ILE"), "ILE", 29, Arrays.asList("ILE:8",  "ALA:1",  "VAL:3", "LEU:5", "PHE:4", "TYR:8"));
			assertPosition(confSpace, Arrays.asList("A193", "193 THR"), "THR", 44, Arrays.asList("THR:19", "SER:18", "ASN:7"));
		};

		// check chain G
		assertThat(chainG.numPos(), is(3));
		assertG.accept(chainG);
		assertThat(chainG.countSingles(), is(57));
		assertThat(chainG.countPairs(), is(1008));

		// check chain A
		assertThat(chainA.numPos(), is(4));
		assertA.accept(chainA);
		assertThat(chainA.countSingles(), is(143));
		assertThat(chainA.countPairs(), is(7575));

		// check the complex
		assertThat(complex.numPos(), is(7));
		assertG.accept(complex);
		assertA.accept(complex);
		assertThat(complex.countSingles(), is(200));
		assertThat(complex.countPairs(), is(16734));
	}

	private static void assertPosition(ConfSpaceIteration confSpace, List<String> names, String wildType, int numConfs, List<String> descriptions) {

		// find the design position
		int posi = IntStream.range(0, confSpace.numPos())
			.filter(i -> names.contains(confSpace.name(i)))
			.findFirst()
			.orElse(-1);
		assertThat("No position with names " + names, posi, is(not(-1)));

		// find the sequence position
		SeqSpace.Position seqPos = confSpace.seqSpace().positions.stream()
			.filter(p -> names.contains(p.resNum))
			.findFirst()
			.orElse(null);
		assertThat("No sequence position with names " + names, seqPos, is(not(nullValue())));

		// check the wild type
		assertThat(confSpace.wildType(posi), is(wildType));
		assertThat(seqPos.wildType.name, is(wildType));

		// check the types
		String[] resTypes = descriptions.stream()
			.map(desc -> desc.split(":")[0])
			.toArray(String[]::new);
		assertThat(
			IntStream.range(0, confSpace.numConf(posi))
				.mapToObj(confi -> confSpace.confType(posi, confi))
				.collect(Collectors.toSet()), // collapse duplicates
			containsInAnyOrder(resTypes)
		);
		assertThat(
			seqPos.resTypes.stream()
				.map(resType -> resType.name)
				.collect(Collectors.toList()),
			containsInAnyOrder(resTypes)
		);

		// count the conformations for each type
		int[] numConfsByDesc = descriptions.stream()
			.mapToInt(desc -> Integer.parseInt(desc.split(":")[1]))
			.toArray();
		for (int i=0; i<descriptions.size(); i++) {
			int fi = i; // javac is dumb ...

			int[] confIndices = IntStream.range(0, confSpace.numConf(posi))
				.filter(confi -> confSpace.confType(posi, confi).equals(resTypes[fi]))
				.toArray();

			assertThat(confIndices.length, is(numConfsByDesc[i]));
		}

		// check the total number of conformations
		assertThat(confSpace.numConf(posi), is(numConfs));
	}

	private static void assert2RL0WildType(Molecule mol) {

		assertThat(mol.residues.size(), is(106));

		int numAtoms = mol.residues.stream()
			.mapToInt(res -> res.atoms.size())
			.sum();
		assertThat(numAtoms, is(1602));
	}
}
