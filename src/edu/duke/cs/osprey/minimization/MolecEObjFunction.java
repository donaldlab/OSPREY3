/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.EllipseCoordDOF;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.dof.StrandRotation;
import edu.duke.cs.osprey.dof.StrandTranslation;
import edu.duke.cs.osprey.dof.deeper.perts.Backrub;
import edu.duke.cs.osprey.dof.deeper.perts.Shear;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.ematrix.epic.EPICEnergyFunction;
import edu.duke.cs.osprey.structure.Molecule;

import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 *
 * @author mhall44
 */
public class MolecEObjFunction implements ObjectiveFunction {
    //we apply the confDOFs to the molecule and then evaluate the energy function
    //so this objective function maps DOF values to energy function values
    //DOF values bounded by constraints
    
    EnergyFunction efunc;
    Molecule molec;
    ArrayList<DegreeOfFreedom> DOFs;
    DoubleMatrix1D[] constraints;
    
    
    DoubleMatrix1D curDOFVals;
    /*
    DOFDEPENDENCIES;
    (see perturbations);
    
    maybe keep an unmoved molecule;
    */
    
    ArrayList<EnergyFunction> partialEFuncs = null;//if not null, can use when searching along a single DOF
    
    public MolecEObjFunction(EnergyFunction ef, DoubleMatrix1D[] constr, Molecule m, 
            ArrayList<DegreeOfFreedom> DOFList){
        
        efunc = ef;
        constraints = constr;
        molec = m;
        DOFs = DOFList;
        
        curDOFVals = DoubleFactory1D.dense.make(DOFs.size());
        
        initEPIC();
    }
    
    
    public MolecEObjFunction(EnergyFunction ef, ConfSpace cSpace, RCTuple RCTup){
        /*Initialize an objective function to evaluate ef over the portion of cSpace
         * defined by the RCs in RCTup.  Ensure that all confDOFs of residues in RCTup are bounded
         * (if able to vary continuously) or set correctly (if not)
         */
        
        efunc = ef;
        molec = cSpace.m;
        
        
        LinkedHashMap<DegreeOfFreedom,double[]> DOFBounds = new LinkedHashMap<>();//bounds for each conformational DOF
        //LinkedHashMap used to achieve consistency between runs (iterating over a regular HashMap
        //would yield a different order from run to run depending on what DegreeOfFreedom pointers are available)
        
        int numMinDOFs = 0;//number of minimizable confDOFs (bounded but not to a single value)
        
        for(int indexInTup=0; indexInTup<RCTup.RCs.size(); indexInTup++){
            
            int posNum = RCTup.pos.get(indexInTup);
            int RCNum = RCTup.RCs.get(indexInTup);
            RC rc = cSpace.posFlex.get(posNum).RCs.get(RCNum);
            
            //make sure the amino-acid type is set correctly
            ResidueTypeDOF mutDOF = cSpace.mutDOFs.get(posNum);
            if( ! mutDOF.getCurResType().equalsIgnoreCase(rc.AAType) ){
                mutDOF.mutateTo(rc.AAType);
            }
            
            
            for(int dofIndexInRC=0; dofIndexInRC<rc.DOFs.size(); dofIndexInRC++){
                
                //get the DOF bounds
                double maxVal = rc.DOFmax.get(dofIndexInRC);
                double minVal = rc.DOFmin.get(dofIndexInRC);
                
                DegreeOfFreedom curDOF = rc.DOFs.get(dofIndexInRC);
                
                //make sure DOF bounds don't contradict bounds from some other RC
                if(DOFBounds.containsKey(curDOF)){
                    double[] prevBounds = DOFBounds.get(curDOF);
                    if(prevBounds[0]!=minVal || prevBounds[1]!=maxVal){
                        throw new RuntimeException("ERROR: Disagreement in DOF bounds between RCs!");
                    }
                }
                else {//store bounds
                    DOFBounds.put(curDOF, new double[] {minVal,maxVal});
                    numMinDOFs++;
                }
            }
        }
        
        //collect constraints, and apply fixed DOF valuestrue
        DOFs = new ArrayList<>();
        constraints = new DoubleMatrix1D[] 
            {DoubleFactory1D.dense.make(numMinDOFs), DoubleFactory1D.dense.make(numMinDOFs)};
        
        int minDOFCount = 0;
        
        for(DegreeOfFreedom dof : DOFBounds.keySet()){
            double bounds[] = DOFBounds.get(dof);
            
            if(bounds[0]==bounds[1]){//fixed DOF
                dof.apply(bounds[0]);//apply fixed value
            }
            else {//minimizable DOF: record DOF and constraints for this objective function
                constraints[0].set(minDOFCount, bounds[0]);
                constraints[1].set(minDOFCount, bounds[1]);
                DOFs.add(dof);
                minDOFCount++;
            }
        }
        
        curDOFVals = DoubleFactory1D.dense.make(DOFs.size());
        
        
        initEPIC();
    }
    
    private void initEPIC(){//Initialize EPIC stuff if needed
        if(efunc instanceof EPICEnergyFunction){
            ((EPICEnergyFunction)efunc).assignConfReference(curDOFVals,DOFs,molec);
            
            //let's make partial energy functions too for speed...
            partialEFuncs = ((EPICEnergyFunction)efunc).getDOFPartialEFuncs(DOFs,molec);
        }
    }
    
    @Override
    public int getNumDOFs() {
        return DOFs.size();
    }

    
    @Override
    public DoubleMatrix1D[] getConstraints() {
        return constraints;
    }

    
    @Override
    public void setDOFs(DoubleMatrix1D x) {
        
        curDOFVals.assign(x);
        
        if(x.size()!=DOFs.size())
            throw new RuntimeException("ERROR: Trying to set "+DOFs.size()+" DOFs with "+x.size()+" values");
        
        for(int dof=0; dof<x.size(); dof++){
            DOFs.get(dof).apply(x.get(dof));
        }
    }
    

    @Override
    public void setDOF(int dof, double val) {
        
        curDOFVals.set(dof, val);
        
        //if(hasDependencies)//will be needed for DEEPer!  Though a provisional solution would be to set all confDOFs every time
        //    undo();
        
        DOFs.get(dof).apply(val);
        
        //redoDependencies();
    }

    
    @Override
    public double getValue(DoubleMatrix1D x) {
        setDOFs(x);
        return efunc.getEnergy();
    }

    @Override
    public double getValForDOF(int dof, double val) {
        
        setDOF(dof,val);
        
        if(partialEFuncs != null)
            return partialEFuncs.get(dof).getEnergy();
        
        return efunc.getEnergy();
    }

    @Override
    public double getInitStepSize(int dof) {
        DegreeOfFreedom curDOF = DOFs.get(dof);
        if(curDOF instanceof FreeDihedral)
            return 0.25;
        else if (curDOF instanceof EllipseCoordDOF) {
        	EllipseCoordDOF e = (EllipseCoordDOF) curDOF;
        	return (e.getIndex()==0) ? 10 : 0.3; 
        }
        else if(curDOF instanceof StrandRotation)
            return 0.0625;
        else if(curDOF instanceof StrandTranslation)
            return 0.025;
        else if(curDOF instanceof Shear)
            return 0.125;
        else if(curDOF instanceof Backrub)
            return 0.125;
        else if(curDOF instanceof BBFreeDOF)
            return 0.05;
        else
            throw new UnsupportedOperationException("ERROR: DOF type not recognized for step size purposes");
    }

    @Override
    public boolean isDOFAngle(int dof) {
        
        DegreeOfFreedom curDOF = DOFs.get(dof);
        
        if(curDOF instanceof FreeDihedral)
            return true;
        else if(curDOF instanceof StrandRotation)
            return true;
        else if(curDOF instanceof StrandTranslation)
            return false;
        else if(curDOF instanceof Shear)
            return true;
        else if(curDOF instanceof Backrub)
            return true;
        
        throw new UnsupportedOperationException("Degree of freedom type not support here yet: "+DOFs.get(dof).toString());
    }
    
    
    //Will likely want gradient and Hessian too...

    public EnergyFunction getEfunc() {
        return efunc;
    }

    public void setEfunc(EnergyFunction efunc) {
        this.efunc = efunc;
    }

    public Molecule getMolec() {
        return molec;
    }

    public ArrayList<DegreeOfFreedom> getDOFs() {
        return DOFs;
    }
    
    
    
}
