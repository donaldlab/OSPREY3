package edu.duke.cs.osprey.sparse;

import java.util.Iterator;
import edu.duke.cs.osprey.confspace.RCTuple;

/***
 * This class enumerates conformations in order. 
 * @author Jon
 *
 */

public class SparseConformationEnumerator implements ConformationProcessor, Iterator<RCTuple>{
	
	
	public SparseConformationEnumerator(BranchDecomposedProblem searchProblem)
	{
		initializeRecursiveEnumerators(searchProblem);
	}

	private void initializeRecursiveEnumerators (BranchDecomposedProblem searchProblem) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void processConformation (RCTuple conformation) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public boolean hasNext () {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public RCTuple next () {
		// TODO Auto-generated method stub
		return null;
	}

}
