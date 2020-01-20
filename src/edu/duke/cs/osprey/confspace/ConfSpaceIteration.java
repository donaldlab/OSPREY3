package edu.duke.cs.osprey.confspace;

import java.util.stream.IntStream;


/**
 * An interface to a conformation space that only exposes
 * enough information to iterate over design positions and conformations.
 *
 * But it doesn't reveal details of how the conformation space is stored,
 * how forcefields are configured, or how atomic coordinates are defined.
 */
public interface ConfSpaceIteration {

	/** counts all the unique pos:conf entries in the conf space */
	int countSingles();

	/** counts all the unique pos:conf - pos:conf pairs in the conf space */
	int countPairs();

	/** returns the number of design positions in the conf space */
	int numPos();

	/** returns the number of conformations at the given design position */
	int numConf(int posi);

	default int[] numConfsByPos() {
		return IntStream.range(0, numPos())
			.map(posi -> numConf(posi))
			.toArray();
	}

	/** returns the name of the position */
	String name(int posi);

	/** returns the identifier (unique at this position) of the conformation */
	String confId(int posi, int confi);

	/** returns the type of the conformation at the position */
	String confType(int posi, int confi);

	/** returns the sequence space for this conformation space */
	SeqSpace seqSpace();

	/** returns the wild type for the position */
	String wildType(int posi);

	/** returns true iff this position has conformations whose types are not the wild type */
	boolean hasMutations(int posi);

	/** creates a new sequence out of assignments from this conf space */
	default Sequence makeSequenceFromAssignments(int[] assignments) {
		Sequence seq = seqSpace().makeUnassignedSequence();
		for (int posi=0; posi<numPos(); posi++) {
			seq.set(name(posi), confType(posi, assignments[posi]));
		}
		return seq;
	}

	/** creates a new sequence out of conformation from this conf space */
	default Sequence makeSequenceFromConf(ConfSearch.ScoredConf conf) {
		return makeSequenceFromAssignments(conf.getAssignments());
	}
}
