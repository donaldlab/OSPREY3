package edu.duke.cs.osprey.astar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.pruning.AStarSequencePruner;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.BeforeClass;
import org.junit.Test;

public class TestSequencePruner {

	private static SimpleConfSpace confSpace;
	private static EnergyMatrix emat;

	@BeforeClass
	public static void beforeClass() {

		// get any arbitrary conf space...
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A5").setLibraryRotamers("GLY", "ALA");
		strand.flexibility.get("A7").setLibraryRotamers("GLY", "ALA");
		strand.flexibility.get("A9").setLibraryRotamers("GLY", "ALA");

		confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		// compute an emat
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			 .setParallelism(Parallelism.makeCpu(4))
			 .build()) {

			emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();
		}
	}

	@Test
	public void removeAllSequences() {

		AStarSequencePruner pruner = new AStarSequencePruner(confSpace);

		ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
			.setTraditional()
			.setPruner(pruner)
			.build();

		assertThat(astar.getNumConformations().intValueExact(), is(8));

		// remove the sequence of the first conf 8 times
		// then there should be no confs left

		Sequence sequence = confSpace.makeSequenceFromConf(astar.nextConf());
		assertThat(sequence.toString(Sequence.Renderer.ResType), is("ALA ALA ALA"));
		pruner.add(sequence);

		sequence = confSpace.makeSequenceFromConf(astar.nextConf());
		assertThat(sequence.toString(Sequence.Renderer.ResType), is("GLY ALA ALA"));
		pruner.add(sequence);

		sequence = confSpace.makeSequenceFromConf(astar.nextConf());
		assertThat(sequence.toString(Sequence.Renderer.ResType), is("ALA GLY ALA"));
		pruner.add(sequence);

		sequence = confSpace.makeSequenceFromConf(astar.nextConf());
		assertThat(sequence.toString(Sequence.Renderer.ResType), is("ALA ALA GLY"));
		pruner.add(sequence);

		sequence = confSpace.makeSequenceFromConf(astar.nextConf());
		assertThat(sequence.toString(Sequence.Renderer.ResType), is("GLY GLY ALA"));
		pruner.add(sequence);

		sequence = confSpace.makeSequenceFromConf(astar.nextConf());
		assertThat(sequence.toString(Sequence.Renderer.ResType), is("GLY ALA GLY"));
		pruner.add(sequence);

		sequence = confSpace.makeSequenceFromConf(astar.nextConf());
		assertThat(sequence.toString(Sequence.Renderer.ResType), is("ALA GLY GLY"));
		pruner.add(sequence);

		sequence = confSpace.makeSequenceFromConf(astar.nextConf());
		assertThat(sequence.toString(Sequence.Renderer.ResType), is("GLY GLY GLY"));
		pruner.add(sequence);

		assertThat(astar.nextConf(), is(nullValue()));
	}
}
