package edu.duke.cs.osprey.gmec;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.Test;

import java.util.*;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;


public class TestComets {

	private static ForcefieldParams ffparams = new ForcefieldParams();

	private static Comets make2RL0Tiny(boolean boundedMemory) {

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("G648", "G654")
			.build();
		protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the COMETS states
		Comets.State unbound = new Comets.State(
			"Unbound",
			new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build()
		);

		Comets.State bound = new Comets.State(
			"Bound",
			new SimpleConfSpace.Builder()
				.addStrand(protein)
				.addStrand(ligand)
				.build()
		);

		// configure COMETS
		Comets.LME objective = new Comets.LME.Builder()
			.addState(bound, 1.0)
			.addState(unbound, -1.0)
			.build();
		Comets comets = new Comets.Builder(objective)
			.setMinNumConfTrees(boundedMemory ? 5 : null)
			.build();

		initStates(comets.states, boundedMemory);

		return comets;
	}

	private static void check2RL0Tiny(Comets comets) {
		List<Comets.SequenceInfo> sequences = comets.findBestSequences(2);
		assertSequence(comets, sequences, "ASP", -15.11711536, new double [] {  });
		assertSequence(comets, sequences, "GLU", -13.56001913, new double [] {  });
		assertThat(sequences.size(), is(2));
		checkSequencesOrder(sequences);
	}

	private static Comets make2RL0Small(boolean boundedMemory) {

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("G648", "G654")
			.build();
		protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "ASP").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType, "SER", "ASN", "GLN").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		// make the COMETS states
		Comets.State unbound = new Comets.State(
			"Unbound",
			new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build()
		);

		Comets.State bound = new Comets.State(
			"Bound",
			new SimpleConfSpace.Builder()
				.addStrand(protein)
				.addStrand(ligand)
				.build()
		);

		// configure COMETS
		Comets.LME objective = new Comets.LME.Builder()
			.addState(bound, 1.0)
			.addState(unbound, -1.0)
			.build();
		Comets comets = new Comets.Builder(objective)
			.setMinNumConfTrees(boundedMemory ? 5 : null)
			.build();

		initStates(comets.states, boundedMemory);

		return comets;
	}

	private static void check2RL0Small(Comets comets) {
		List<Comets.SequenceInfo> sequences = comets.findBestSequences(11);
		assertSequence(comets, sequences, "PHE ASP GLU GLN", -64.51492231, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU ASN", -63.91337074, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU SER", -63.70306142, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU THR", -62.76135586, new double [] {  });
		assertSequence(comets, sequences, "PHE GLU GLU THR", -61.11651357, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP ASP THR", -60.41594521, new double [] {  });
		assertSequence(comets, sequences, "ILE ASP GLU THR", -58.93722354, new double [] {  });
		assertSequence(comets, sequences, "VAL ASP GLU THR", -58.87767041, new double [] {  });
		assertSequence(comets, sequences, "LEU ASP GLU THR", -58.53375330, new double [] {  });
		assertSequence(comets, sequences, "ALA ASP GLU THR", -57.74755836, new double [] {  });
		assertSequence(comets, sequences, "TYR ASP GLU THR", -57.35573519, new double [] {  });
		assertThat(sequences.size(), is(11));
		checkSequencesOrder(sequences);
	}

	private static Comets make2RL0PPI(boolean boundedMemory) {

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("G648", "G654")
			.build();
		protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "ASP").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType, "SER", "ASN", "GLN").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType, "ASP", "GLU").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType, "ALA", "VAL", "LEU", "TYR").addWildTypeRotamers().setContinuous();
		ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType, "SER", "ASN").addWildTypeRotamers().setContinuous();

		// make the COMETS states
		Comets.State unbound = new Comets.State(
			"Unbound",
			new SimpleConfSpace.Builder()
				.addStrand(protein)
				.build()
		);

		Comets.State bound = new Comets.State(
			"Bound",
			new SimpleConfSpace.Builder()
				.addStrand(protein)
				.addStrand(ligand)
				.build()
		);

		// configure COMETS
		Comets.LME objective = new Comets.LME.Builder()
			.addState(bound, 1.0)
			.addState(unbound, -1.0)
			.build();
		Comets comets = new Comets.Builder(objective)
			.setMaxSimultaneousMutations(1)
			.setObjectiveWindowMax(2000) // need a big window to get all the sequences
			.setObjectiveWindowSize(10000)
			.setMinNumConfTrees(boundedMemory ? 5 : null)
			.build();

		initStates(comets.states, boundedMemory);

		return comets;
	}

	private static void check2RL0PPI(Comets comets) {
		List<Comets.SequenceInfo> sequences = comets.findBestSequences(24);
		assertSequence(comets, sequences, "PHE ASP GLU THR PHE LYS ILE THR", -62.76135586, new double [] {  });
		assertSequence(comets, sequences, "TYR ASP GLU THR PHE LYS ILE THR", -57.35573519, new double [] {  });
		assertSequence(comets, sequences, "ALA ASP GLU THR PHE LYS ILE THR", -57.74755836, new double [] {  });
		assertSequence(comets, sequences, "VAL ASP GLU THR PHE LYS ILE THR", -58.87767034, new double [] {  });
		assertSequence(comets, sequences, "ILE ASP GLU THR PHE LYS ILE THR", -58.93722356, new double [] {  });
		assertSequence(comets, sequences, "LEU ASP GLU THR PHE LYS ILE THR", -58.53375350, new double [] {  });
		assertSequence(comets, sequences, "PHE GLU GLU THR PHE LYS ILE THR", -61.11651357, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP ASP THR PHE LYS ILE THR", -60.41594521, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU SER PHE LYS ILE THR", -63.70306141, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU ASN PHE LYS ILE THR", -63.91337074, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU GLN PHE LYS ILE THR", -64.51492231, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU THR TYR LYS ILE THR", -62.35990315, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU THR ALA LYS ILE THR", -58.55409116, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU THR VAL LYS ILE THR", -61.18841895, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU THR ILE LYS ILE THR", -62.45662407, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU THR LEU LYS ILE THR", -58.06841434, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU THR PHE ASP ILE THR", -43.31822569, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU THR PHE GLU ILE THR", -42.28743004, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU THR PHE LYS ALA THR", -56.26178013, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU THR PHE LYS VAL THR", -59.16576214, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU THR PHE LYS LEU THR", -5.30989568, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU THR PHE LYS TYR THR", 1733.69997029, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU THR PHE LYS ILE SER", -62.28347340, new double [] {  });
		assertSequence(comets, sequences, "PHE ASP GLU THR PHE LYS ILE ASN", -60.68742687, new double [] {  });
		assertThat(sequences.size(), is(24));
		checkSequencesOrder(sequences);
	}

	private static void initStates(List<Comets.State> states, boolean boundedMemory) {

		// make the ecalc from all the conf spaces
		List<SimpleConfSpace> confSpaces = states.stream()
			.map(state -> state.confSpace)
			.collect(Collectors.toList());
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			for (Comets.State state : states) {

				// calculate energy matrices
				SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(state.confSpace, ecalc)
					.build()
					.calcReferenceEnergies();

				state.confEcalc = new ConfEnergyCalculator.Builder(state.confSpace, ecalc)
					.setReferenceEnergies(eref)
					.build();

				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(state.confEcalc)
					.build()
					.calcEnergyMatrix();
				state.fragmentEnergies = emat;

				PruningMatrix pmat = new SimpleDEE.Runner()
					.setTypeDependent(true)
					.setGoldsteinDiffThreshold(10.0)
					.run(state.confSpace, emat);

				// make the conf tree factory
				state.confTreeFactory = (rcs) -> new ConfAStarTree.Builder(emat, new RCs(rcs, pmat))
					.setMaxNumNodes(boundedMemory ? 100000L : null)
					.setTraditional()
					.build();
			}
		}
	}

	private static void prepStates(Comets comets, Runnable block) {

		// make the ecalc from all the conf spaces
		List<SimpleConfSpace> confSpaces = comets.states.stream()
			.map(state -> state.confSpace)
			.collect(Collectors.toList());
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces, ffparams)
			.setParallelism(Parallelism.makeCpu(8))
			.build()) {

			// refresh the conf ecalcs
			for (Comets.State state : comets.states) {
				state.confEcalc = new ConfEnergyCalculator(state.confEcalc, ecalc);
			}

			block.run();
		}
	}

	public static void main(String[] args) {
		bruteForce("2RL0 Tiny", make2RL0Tiny(false));
		bruteForce("2RL0 Small", make2RL0Small(false));
		bruteForce("2RL0 PPI", make2RL0PPI(false));
	}

	public static void bruteForce(String name, Comets comets) {

		log("\n\n%s\n", name);

		// explicitly enumerate all the sequences
		List<Sequence> sequences = new ArrayList<>();
		sequences.add(comets.seqSpace.makeWildTypeSequence());
		sequences.addAll(comets.seqSpace.getMutants(1));

		prepStates(comets, () -> {

			for (Sequence sequence : sequences) {

				// calculate all the GMECs
				Map<Comets.State,Double> stateEnergies = new LinkedHashMap<>();
				for (Comets.State state : comets.states) {

					// do a GMEC search
					RCs rcs = sequence
						.filter(state.confSpace.seqSpace)
						.makeRCs(state.confSpace);
					ConfAStarTree astar = state.confTreeFactory.apply(rcs);
					ConfSearch.EnergiedConf gmec = new SimpleGMECFinder.Builder(astar, state.confEcalc)
						.setPrintToConsole(false)
						.build()
						.find();

					stateEnergies.put(state, gmec.getEnergy());
				}

				// compute the objective and constraints, e.g.:
				// assertSequence(comets, sequences, "asp", -15.11711536, new double [] { -10.0 });
				log("assertSequence(comets, sequences, \"%s\", %.8f, new double [] { %s });",
					sequence.toString(Sequence.Renderer.ResType),
					comets.objective.calc(stateEnergies),
					String.join(", ", comets.constraints.stream()
						.map(constraint -> String.format("%.8f", constraint.calc(stateEnergies)))
						.collect(Collectors.toList())
					)
				);
			}
		});
	}

	@Test
	public void tiny2RL0() {
		Comets comets = make2RL0Tiny(false);
		prepStates(comets, () -> check2RL0Tiny(comets));
	}

	@Test
	public void small2RL0() {
		Comets comets = make2RL0Small(false);
		prepStates(comets, () -> check2RL0Small(comets));
	}

	@Test
	public void ppi2RL0() {
		Comets comets = make2RL0PPI(false);
		prepStates(comets, () -> check2RL0PPI(comets));
	}

	@Test
	public void ppi2RL0BoundedMemory() {
		Comets comets = make2RL0PPI(true);
		prepStates(comets, () -> check2RL0PPI(comets));
	}

	public static void assertSequence(Comets comets, List<Comets.SequenceInfo> sequences, String seqStr, double objective, double[] constraints) {

		// reconstruct the sequence
		Sequence sequence = comets.seqSpace.makeSequence(Arrays.asList(seqStr.split(" ")));

		// find the sequence info (or die trying)
		Comets.SequenceInfo info = sequences.stream()
			.filter(i -> i.sequence.equals(sequence))
			.findFirst()
			.orElseThrow(() -> new NoSuchElementException("can't find sequence " + seqStr));

		// check the results
		final double epsilon = 1e-6;
		assertThat(info.sequence, is(sequence));
		assertThat(info.objective, isAbsolutely(objective, epsilon));
		for (int i=0; i<comets.constraints.size(); i++) {
			Comets.LME constraint = comets.constraints.get(i);
			assertThat(info.constraints.get(constraint), isAbsolutely(constraints[i], epsilon));
		}
	}

	public static void checkSequencesOrder(List<Comets.SequenceInfo> sequences) {

		// they should be weakly increasing by objective value
		for (int i=1; i<sequences.size(); i++) {
			assertThat(sequences.get(i - 1).objective, lessThanOrEqualTo(sequences.get(i).objective));
		}
	}
}
