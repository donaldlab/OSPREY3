package edu.duke.cs.osprey.sparse;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.Minimizer.Result;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;

public class PartialConformationEnergyFunction {
	
	private EnergyFunction fullEnergyFunction;
	private ConfSpace conformations;
	
	public PartialConformationEnergyFunction(EnergyFunction termE, ConfSpace conformationSpace)
	{
		fullEnergyFunction = termE;
		conformations = conformationSpace;
	}
	
	
	public double computePartialEnergy(RCTuple priorConformation, RCTuple partialAssignment)
	{
		RCTuple combinedAssignment = new RCTuple();
		combinedAssignment.set(priorConformation);
		combinedAssignment.set(partialAssignment);
		MoleculeModifierAndScorer mof = new MoleculeModifierAndScorer(fullEnergyFunction,conformations,combinedAssignment);

        DoubleMatrix1D bestDOFVals;

        if(mof.getNumDOFs()>0){//there are continuously flexible DOFs to minimize
            CCDMinimizer ccdMin = new CCDMinimizer(mof,true);
            bestDOFVals = ccdMin.minimize().dofValues;
        }
        else//molecule is already in the right, rigid conformation
            bestDOFVals = DoubleFactory1D.dense.make(0);


        return mof.getValue(bestDOFVals);
	}

}
