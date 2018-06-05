package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
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
import edu.duke.cs.osprey.tools.MathTools;
import org.junit.BeforeClass;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;


public class TestComets {

	private static ForcefieldParams ffparams;
	private static Comets comets2RL0Tiny;

	@BeforeClass
	public static void beforeClass() {
		ffparams = new ForcefieldParams();
		comets2RL0Tiny = make2RL0Tiny();
	}

	private static Comets make2RL0Tiny() {

		Molecule mol = PDBIO.readResource("/2RL0.min.reduce.pdb");

		// make sure all strands share the same template library
		// especially since all the state conf spaces will add/share wild-type rotamers to/from the library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ffparams.forcefld).build();

		// define the protein strand
		Strand protein = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("G648", "G654")
			.build();
		//protein.flexibility.get("G649").setLibraryRotamers(Strand.WildType, "TYR", "ALA", "VAL", "ILE", "LEU").addWildTypeRotamers().setContinuous();
		protein.flexibility.get("G650").setLibraryRotamers(Strand.WildType, "GLU").addWildTypeRotamers().setContinuous();
		//protein.flexibility.get("G651").setLibraryRotamers(Strand.WildType, "ASP").addWildTypeRotamers().setContinuous();
		//protein.flexibility.get("G654").setLibraryRotamers(Strand.WildType, "SER", "ASN", "GLN").addWildTypeRotamers().setContinuous();

		// define the ligand strand
		Strand ligand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.setResidues("A155", "A194")
			.build();
		ligand.flexibility.get("A156").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		//ligand.flexibility.get("A172").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		//ligand.flexibility.get("A192").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		//ligand.flexibility.get("A193").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

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
			.build();

		// make the ecalc from all the conf spaces
		List<SimpleConfSpace> confSpaces = Arrays.asList(unbound.confSpace, bound.confSpace);
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces, ffparams)
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			for (Comets.State state : comets.states) {

				// calculate energy matrices
				SimpleReferenceEnergies eref = new SimplerEnergyMatrixCalculator.Builder(state.confSpace, ecalc)
					.build()
					.calcReferenceEnergies();

				state.confEcalc = new ConfEnergyCalculator.Builder(state.confSpace, ecalc)
					.setReferenceEnergies(eref)
					.build();

				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(state.confEcalc)
					.setCacheFile(new File(String.format("emat.%s.dat", state.name)))
					.build()
					.calcEnergyMatrix();
				state.fragmentEnergies = emat;

				// do DEE
				PruningMatrix pmat = new SimpleDEE.Runner()
					.setCacheFile(new File(String.format("pmat.%s.dat", state.name)))
					.run(state.confSpace, emat);

				// make the conf tree factory
				state.confTreeFactory = (rcs) -> new ConfAStarTree.Builder(emat, new RCs(rcs, pmat))
					.setTraditional()
					.build();
			}
		}

		return comets;
	}

	public static void main(String[] args) {

		beforeClass();

		for (Comets comets : Arrays.asList(comets2RL0Tiny)) {

			SimpleConfSpace confSpace = comets.states.get(0).confSpace;

			// explicitly enumerate all the sequences
			List<Sequence> sequences = new ArrayList<>();
			sequences.add(confSpace.makeWildTypeSequence());
			for (List<SimpleConfSpace.Position> mutablePositions : MathTools.powerset(confSpace.mutablePositions)) {

				// collect all the mutations away from wild-type
				List<List<String>> resTypes = mutablePositions.stream()
					.map(mpos -> mpos.resTypes)
					.collect(Collectors.toList());

				// enumerate all the combinations of res types
				for (List<String> mutations : MathTools.cartesianProduct(resTypes)) {

					/* TODO
					Sequence seq = confSpace.makeWildTypeSequence();
					for (int i=0; i<mutablePositions.size(); i++) {
						seq.set(mutablePositions.get(i), mutations.get(i));
					}
					sequences.add(seq);
					*/
				}
			}
		}
	}
}
