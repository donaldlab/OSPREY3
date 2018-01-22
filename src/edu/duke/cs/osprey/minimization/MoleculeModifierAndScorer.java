/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.minimization;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
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
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

/**
 *
 * @author mhall44
 */
public class MoleculeModifierAndScorer implements ObjectiveFunction {
	//we apply the confDOFs to the molecule and then evaluate the energy function
	//so this objective function maps DOF values to energy function values
	//DOF values bounded by constraints: Motions are intended to be within a voxel

	//Warning: This class has the side effect of modifying the molecule passed to it!
	//It would be very expensive to copy the molecule for each minimization, so we
	//apply the DOFs in place and then evaluate the energy

	private static final long serialVersionUID = 3898313221157632380L;

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

	List<EnergyFunction> partialEFuncs = null;//if not null, can use when searching along a single DOF

	/**
	 * transition adapter, only here temporarily
	 */
	@Deprecated
	public MoleculeModifierAndScorer(MoleculeObjectiveFunction mof) {
		molec = mof.pmol.mol;
		DOFs = new ArrayList<>(mof.pmol.dofs);
		constraints = mof.pmol.dofBounds.getBounds();
		efunc = mof.efunc;
		curDOFVals = mof.curDOFVals;
		partialEFuncs = mof.efuncsByDof;
	}

    public static boolean hasMinimizableDofs(ConfSpace confSpace, RCTuple tuple) {
    
        // for each pos...
        for (int i=0; i<tuple.size(); i++) {
            
            int pos = tuple.pos.get(i);
            int rc = tuple.RCs.get(i);
            RC rcObj = confSpace.posFlex.get(pos).RCs.get(rc);
            
            for (int d=0; d<rcObj.DOFs.size(); d++) {
                double xdmin = rcObj.DOFmin.get(d);
                double xdmax = rcObj.DOFmax.get(d);
                if (xdmax > xdmin) {
                    return true;
                }
            }
        }
        
        return false;
    }
    
    public MoleculeModifierAndScorer(EnergyFunction ef, DoubleMatrix1D[] constr, Molecule m, 
            ArrayList<DegreeOfFreedom> DOFList){
        
        constraints = constr;
        molec = m;
        DOFs = DOFList;
        
        curDOFVals = DoubleFactory1D.dense.make(DOFs.size());
        
        setEfunc(ef);
    }
    
    
    public MoleculeModifierAndScorer(EnergyFunction ef, ConfSpace cSpace, RCTuple RCTup) {
        this(ef, cSpace, RCTup, null);
    }
    
    public MoleculeModifierAndScorer(EnergyFunction ef, ConfSpace cSpace, RCTuple RCTup, 
            ParameterizedMoleculeCopy pmc) {
        /*Initialize an objective function to evaluate ef over the portion of cSpace
         * defined by the RCs in RCTup.  Ensure that all confDOFs of residues in RCTup are bounded
         * (if able to vary continuously) or set correctly (if not)
         */
        
        // TODO: copying the DOFs here and assigning them to the specified mol instance is pretty hacky.
        // in the future, it would be nice to do a bigger refactor so that anything that modifies
        // molecules should expect a molecule as an argument, rather than expecting to find one
        // in a near-global location like the ConfSpace. we'd also have to separate config from state
        // in the ConfSpace object (ie, take out the molecule entirely)
    	
    	// TODO: energy functions often need to have their forcefields immediately rebuilt after
    	// this constructor changes residue templates which wastes a lot of work. Ideally,
    	// we wouldn't build the energy function until after this constructor is done changing
    	// templates, but that refactor will have to wait for another day.
        
        // which molecule are we using?
        if (pmc == null) {
            
            // the one from the conf space
            this.molec = cSpace.m;
            
        } else {
            
            // a separate molecule, so we don't modify the one in the conf space
            this.molec = pmc.getCopiedMolecule();
        }
        
        LinkedHashMap<DegreeOfFreedom,double[]> DOFBounds = new LinkedHashMap<>();//bounds for each conformational DOF
        //LinkedHashMap used to achieve consistency between runs (iterating over a regular HashMap
        //would yield a different order from run to run depending on what DegreeOfFreedom pointers are available)
        
        int numMinDOFs = 0;//number of minimizable confDOFs (bounded but not to a single value)
        
        for(int indexInTup=0; indexInTup<RCTup.RCs.size(); indexInTup++){
            
            int posNum = RCTup.pos.get(indexInTup);
            int RCNum = RCTup.RCs.get(indexInTup);
            RC rc = cSpace.posFlex.get(posNum).RCs.get(RCNum);
            
            // AAO 2016: this code was written for AAs, specifically anything
            // in all_amino_coords.in and not for generic non-AA residues. skipping this
            // step for non AAs (for now).
            //MH now we can do this for non-AAs as long as coords are available for mutation
            //Judge this by whether coords are available for the current template
            Residue res = cSpace.posFlex.get(posNum).res;
            //if(HardCodedResidueInfo.hasAminoAcidBB(res) && !res.fullName.startsWith("FOL")) {
            if(res.template.templateRes.coords!=null){

                ResidueTypeDOF mutDOF = cSpace.mutDOFs.get(posNum);
                
                // if we're not using the conf space molecule, copy the dof
                if (pmc != null) {
                    mutDOF = (ResidueTypeDOF) pmc.getCopiedDOF(mutDOF);
                }
                
                // make sure the residue is using the right template
                ResidueTemplate desiredTemplate;
                if (rc.template != null) {
                    desiredTemplate = rc.template;
                } else {
                    desiredTemplate = mutDOF.getLibraryTemplate(rc.AAType);
                }

                if (mutDOF.isTemplate(desiredTemplate)) {
                    // restore coords from the template to fight roundoff error
                    mutDOF.restoreCoordsFromTemplate();
                } else {
                    // mutate to the new template (which also uses template coords and fights roundoff error)
                    mutDOF.switchToTemplate(desiredTemplate);
                }
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
        
        // if we're not using the conf space molecule, copy the dofs
        if (pmc != null) {
            LinkedHashMap<DegreeOfFreedom,double[]> copiedDofs = new LinkedHashMap<>();
            for (Map.Entry<DegreeOfFreedom,double[]> entry : DOFBounds.entrySet()) {
                DegreeOfFreedom dof = entry.getKey();
                dof = pmc.getCopiedDOF(dof);
                copiedDofs.put(dof, entry.getValue());
            }
            DOFBounds = copiedDofs;
        }
        
        init(numMinDOFs, DOFBounds);
        setEfunc(ef);
    }
    
    public MoleculeModifierAndScorer(EnergyFunction efunc, ConfSpace confSpace) {
    	
        this.molec = confSpace.m;
        
        int numMinDOFs = 0;
        LinkedHashMap<DegreeOfFreedom,double[]> DOFBounds = new LinkedHashMap<>();
        
        // build the DoFs based on the current structure instead of residue conformations
        for (int i=0; i<confSpace.posFlex.size(); i++) {
            PositionConfSpace pos = confSpace.posFlex.get(i);
            for(int j=0; j<pos.res.getNumDihedrals(); j++) {
                DOFBounds.put(new FreeDihedral(pos.res, j), pos.makeDOFBounds(pos.res.getDihedralAngle(j)));
                numMinDOFs++;
            }
        }
        
        init(numMinDOFs, DOFBounds);
        setEfunc(efunc);
    }
    
    private void init(int numMinDOFs, LinkedHashMap<DegreeOfFreedom,double[]> DOFBounds) {
        
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
    
    public double getCurValueOfDOF(int dof){
        //get the value of the specified DOF
        //(rather than of the objective function)
        return curDOFVals.get(dof);
    }

    @Override
    public double getInitStepSize(int dof) {
    	return getInitStepSize(DOFs.get(dof));
    }
    
    // TODO: move this into DegreeOfFreedom
    public static double getInitStepSize(DegreeOfFreedom dof) {
        if(dof instanceof FreeDihedral)
            return 0.25;
        else if (dof instanceof EllipseCoordDOF) {
            EllipseCoordDOF e = (EllipseCoordDOF) dof;
            return (e.getIndex()==0) ? 10 : 0.3; 
        }
        else if(dof instanceof StrandRotation)
            return 0.0625;
        else if(dof instanceof StrandTranslation)
            return 0.025;
        else if(dof instanceof Shear)
            return 0.125;
        else if(dof instanceof Backrub)
            return 0.125;
        else if(dof instanceof BBFreeDOF)
            return 0.05;
        else
            throw new UnsupportedOperationException("ERROR: DOF type not recognized for step size purposes");
    }

    @Override
    public boolean isDOFAngle(int dof) {
    	return isDOFAngle(DOFs.get(dof));
    }
    
    // TODO: move this into DegreeOfFreedom
    public static boolean isDOFAngle(DegreeOfFreedom dof) {
        if(dof instanceof FreeDihedral)
            return true;
        else if(dof instanceof StrandRotation)
            return true;
        else if(dof instanceof StrandTranslation)
            return false;
        else if(dof instanceof Shear)
            return true;
        else if(dof instanceof Backrub)
            return true;
        else if(dof instanceof BBFreeDOF)
            return false;
        
        throw new UnsupportedOperationException("Degree of freedom type not support here yet: "+dof.toString());
    }
    
    
    //Will likely want gradient and Hessian too...

    public EnergyFunction getEfunc() {
        return efunc;
    }
    
    public EnergyFunction getEfunc(int dof) {
        if(partialEFuncs == null) {
            return null;
        }
        return partialEFuncs.get(dof);
    }

    public void setEfunc(EnergyFunction efunc) {
        
        this.efunc = efunc;
        
        // init efunc if needed
        if (efunc instanceof EnergyFunction.NeedsInit) {
            ((EnergyFunction.NeedsInit)efunc).init(molec, DOFs, curDOFVals);
        }
        
        // decompose by dofs if supported
        if (efunc instanceof EnergyFunction.DecomposableByDof) {
            partialEFuncs = ((EnergyFunction.DecomposableByDof)efunc).decomposeByDof(molec, DOFs);
        } else {
        	partialEFuncs = null;
        }
        
		// we might have made chemical changes, explicitly update the efuncs if supported
        if (efunc instanceof EnergyFunction.ExplicitChemicalChanges) {
        	((EnergyFunction.ExplicitChemicalChanges)efunc).handleChemicalChanges();
        }
        if (partialEFuncs != null) {
			for (EnergyFunction dofEfunc : partialEFuncs) {
				if (dofEfunc instanceof EnergyFunction.ExplicitChemicalChanges) {
					((EnergyFunction.ExplicitChemicalChanges)dofEfunc).handleChemicalChanges();
				}
			}
        }
    }
    
    public Molecule getMolec() {
        return molec;
    }

    public ArrayList<DegreeOfFreedom> getDOFs() {
        return DOFs;
    }
    
    
    
    
    public boolean isOutOfRange(DoubleMatrix1D x){
        if(x.size()!=DOFs.size())
            throw new RuntimeException("ERROR: Trying to check range on "+DOFs.size()+" DOFs with "+x.size()+" values");
        
        for(int dof=0; dof<DOFs.size(); dof++){
            if( x.get(dof) < constraints[0].get(dof)-1e-6 )
                return true;
            if( x.get(dof) > constraints[1].get(dof)+1e-6 )
                return true;
        }
        
        return false;
    }
    
    
    
    
    @Override
    public ArrayList<Integer> getInitFixableDOFs(){
        //currently only going to do this for CATS
        ArrayList<Integer> fixableDOFs = new ArrayList<>();
        for(int dof=0; dof<DOFs.size(); dof++){
            if(DOFs.get(dof) instanceof BBFreeDOF)
                fixableDOFs.add(dof);
        }
        return fixableDOFs;
    }
    
    
}
