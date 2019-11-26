package edu.duke.cs.osprey.confspace.compiled;


/**
 * A single dimension over which a minimizer can explore
 * the energy landscape of a continuous motion.
 *
 * The domain is a finite interval of doubles.
 */
public interface DegreeOfFreedom {

	double min();
	double max();

	double get();
	void set(double val);

	/** which design positions are modified by this DoF? */
	int[] modifiedPosIndices();

	double initialStepSize();
}
