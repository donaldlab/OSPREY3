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

	String confResType(int posi, int confi);
}
