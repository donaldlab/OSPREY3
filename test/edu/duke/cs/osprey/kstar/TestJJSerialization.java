package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.FileTools;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

public class TestJJSerialization {

	public static void main(String[] args) {

		File dir = new File("/home/jeff/donaldlab/osprey test cases/jj-serialization");

		ForcefieldParams ffparams = new ForcefieldParams();

		// read the protein
		Molecule mol = PDBIO.readFile(new File(dir, "gp120SRVRC26.09SR.pdb"));

		// make the template library
		ResidueTemplateLibrary templates = new ResidueTemplateLibrary.Builder(ffparams.forcefld)
			.addTemplates(FileTools.readFile(new File(dir, "all_nuc94_and_gr.in")))
			.addTemplateCoords(FileTools.readFile(new File(dir, "all_amino_coords.in")))
			.addRotamers(FileTools.readFile(new File(dir, "GenericRotamers.dat")))
			.addMoleculeForWildTypeRotamers(mol)
			.build();

		Strand protein = new Strand.Builder(mol)
			.setResidues("H1792", "L2250")
			.setTemplateLibrary(templates)
			.build();
		protein.flexibility.get("H1901").setLibraryRotamers(Strand.WildType).setContinuous();
		protein.flexibility.get("H1906").setLibraryRotamers(Strand.WildType).setContinuous();
		protein.flexibility.get("H1907").setLibraryRotamers(Strand.WildType).setContinuous();
		protein.flexibility.get("H1908").setLibraryRotamers(Strand.WildType).setContinuous();

		List<String> mutations = Arrays.asList((Strand.WildType + " ALA VAL LEU ILE PHE TYR TRP CYS MET SER THR LYS ARG HIS ASP GLU ASN GLN GLY").split(" "));
		protein.flexibility.get("H1904").setLibraryRotamers(mutations).setContinuous();
		protein.flexibility.get("H1905").setLibraryRotamers(mutations).setContinuous();

		Strand ligand = new Strand.Builder(mol)
			.setResidues("F379", "J1791")
			.setTemplateLibrary(templates)
			.build();
		ligand.flexibility.get("G973").setLibraryRotamers(Strand.WildType).setContinuous();
		ligand.flexibility.get("G977").setLibraryRotamers(Strand.WildType).setContinuous();
		ligand.flexibility.get("G978").setLibraryRotamers(Strand.WildType).setContinuous();
		ligand.flexibility.get("G979").setLibraryRotamers(Strand.WildType).setContinuous();
		ligand.flexibility.get("G980").setLibraryRotamers(Strand.WildType).setContinuous();
		ligand.flexibility.get("J1448").setLibraryRotamers(Strand.WildType).setContinuous();

		SimpleConfSpace proteinConfSpace = new SimpleConfSpace.Builder()
			.addStrands(protein)
			.setShellDistance(10.0)
			.build();

		SimpleConfSpace ligandConfSpace = new SimpleConfSpace.Builder()
			.addStrands(ligand)
			.setShellDistance(10.0)
			.build();

		SimpleConfSpace complexConfSpace = new SimpleConfSpace.Builder()
			.addStrands(ligand, protein) // NOTE: need ligand first to match confs.txt
			.setShellDistance(10.0)
			.build();

		// load sequences
		Map<KStar.Sequence,Double> lowEnergySequences = new HashMap<>();
		for (String line : FileTools.parseLines(FileTools.readFile(new File(dir, "confs.txt")))) {

			// e.g.
			// 0 CONF: 18 34 5 27 9 34 5 164 168 0 27 9 RESTYPES: THR ARG ASP LYS LYS ARG ASP GLU GLU GLY LYS GLN ROTS: W0 W0 W0 W0 L9 W0 W0 L3 L7 L W0 W0 Score: -91.23010472425628 Energy: -81.16489688321263 Best so far: -81.16489688321263

			Scanner scanner = new Scanner(line);
			scanner.next();

			// skip confs
			assert (scanner.next().equals("CONF:"));
			for (SimpleConfSpace.Position pos : complexConfSpace.positions) {
				scanner.next();
			}

			// read res types
			assert (scanner.next().equals("RESTYPES:"));
			KStar.Sequence sequence = new KStar.Sequence(complexConfSpace.positions.size());
			for (SimpleConfSpace.Position pos : complexConfSpace.positions) {
				sequence.set(pos.index, scanner.next());
			}

			// skip rots
			assert (scanner.next().equals("ROTS:"));
			for (SimpleConfSpace.Position pos : complexConfSpace.positions) {
				scanner.next();
			}

			// skip score, read energy
			assert (scanner.next().equals("Score:"));
			scanner.next();
			assert (scanner.next().equals("Energy:"));
			double energy = scanner.nextDouble();

			// keep the lowest energy
			if (!lowEnergySequences.containsKey(sequence) || lowEnergySequences.get(sequence) > energy) {
				lowEnergySequences.put(sequence, energy);
			}
		}

		// sort sequences by energy
		System.out.println("Collecting sequences...");
		List<KStar.Sequence> sequences = lowEnergySequences.entrySet().stream()
			.sorted((a, b) -> Double.compare(a.getValue(), b.getValue()))
			.map((entry) -> entry.getKey())
			.collect(Collectors.toList());

		// how should we compute energies of molecules?
		new EnergyCalculator.Builder(complexConfSpace, ffparams)
			//.setParallelism(Parallelism.makeCpu(4))
			.setParallelism(Parallelism.make(4, 1, 16))
			.use((ecalc) -> {

				// how should we define energies of conformations?
				KStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
					return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
						.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
							.build()
							.calcReferenceEnergies()
						).build();
				};

				// how should confs be ordered and searched?
				KStar.ConfSearchFactory confSearchFactory = (emat, pmat) -> {
					return new ConfAStarTree.Builder(emat, pmat)
						//.setTraditional()
						.setMPLP(new ConfAStarTree.MPLPBuilder()
							.setUpdater(new EdgeUpdater())
							.setNumIterations(5)
						)
						.setShowProgress(true)
						.useExternalMemory()
						.build();
				};

				// run K*
				KStar.Settings settings = new KStar.Settings.Builder()
					.setEpsilon(0.68)
					.addScoreConsoleWriter()
					.setShowPfuncProgress(true)
					.setEnergyMatrixCachePattern(new File(dir, "emat.*.dat").getAbsolutePath())
					.build();
				SequenceEnsembleAnalyzer analyzer = new SequenceEnsembleAnalyzer(proteinConfSpace, ligandConfSpace, complexConfSpace, ecalc, confEcalcFactory, confSearchFactory, settings);

				ExternalMemory.use(256, () -> {
					ExternalMemory.setTempDir("/tmp", "osprey");

					// analyze a sequence
					KStar.Sequence sequence = sequences.get(0);
					analyzer.analyze(sequence, 10000);
				});
			});
	}
}
