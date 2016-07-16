package edu.duke.cs.osprey.kstar.emat;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.KSSearchProblem;

@SuppressWarnings("serial")
public class ReducedEnergyMatrix extends EnergyMatrix {

	protected KSSearchProblem sp;

	
	public ReducedEnergyMatrix(KSSearchProblem sp, EnergyMatrix emat) {
		super(emat);
		this.sp = sp;
	}
    
    
    @Override
    public Double getOneBody(int res, int index) {
    	
    	Integer pos = sp.posNums.get(res);
    	
        return super.getOneBody(pos, index);
    }
    
    
    @Override
    public Double getPairwise(int res1, int index1, int res2, int index2) {
    	
    	Integer pos1 = sp.posNums.get(res1);
		Integer pos2 = sp.posNums.get(res2);
		
		return super.getPairwise(pos1, index1, pos2, index2);
    }
    
    
    @Override
	public HigherTupleFinder<Double> getHigherOrderTerms(int res1, int index1, int res2, int index2) {

		Integer pos1 = sp.posNums.get(res1);
		Integer pos2 = sp.posNums.get(res2);

		return super.getHigherOrderTerms(pos1, index1, pos2, index2);
	}
	
	
	@Override
	public int getNumConfAtPos(int pos) {

		Integer pos1 = sp.posNums.get(pos);

		return super.getNumConfAtPos(pos1);
	}


	@Override
	public int getNumPos() {
		return sp.posNums.size();
	}
}
