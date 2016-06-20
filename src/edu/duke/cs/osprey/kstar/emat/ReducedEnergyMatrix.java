package edu.duke.cs.osprey.kstar.emat;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

@SuppressWarnings("serial")
public class ReducedEnergyMatrix extends EnergyMatrix {

	protected EnergyMatrix emat;
	SearchProblem sp;

	
	public ReducedEnergyMatrix(SearchProblem sp, EnergyMatrix emat) {
		this.sp = sp;
		this.emat = emat;
	}
	
	
	public EnergyMatrix getEmat() {
		return emat;
	}
	
	
	@Override
    public double getConstTerm() {
        return emat.getConstTerm();
    }
    
    
    @Override
    public Double getOneBody(int res, int index) {
    	
    	Integer pos = sp.posNums.get(res);
    	
        return emat.oneBody.get(pos).get(index);
    }
    
    
    @Override
    public Double getPairwise(int res1, int index1, int res2, int index2) {
    	
    	Integer pos1 = sp.posNums.get(res1);
		Integer pos2 = sp.posNums.get(res2);
		
		return emat.getPairwise(pos1, index1, pos2, index2);
    }
    
    
    @Override
	public HigherTupleFinder<Double> getHigherOrderTerms(int res1, int index1, int res2, int index2) {

		Integer pos1 = sp.posNums.get(res1);
		Integer pos2 = sp.posNums.get(res2);

		return emat.getHigherOrderTerms(pos1, index1, pos2, index2);
	}
	
	
	@Override
	public int numRCsAtPos(int pos) {

		Integer pos1 = sp.posNums.get(pos);

		return emat.numRCsAtPos(pos1);
	}


	@Override
	public int numPos() {
		return sp.posNums.size();
	}
}
