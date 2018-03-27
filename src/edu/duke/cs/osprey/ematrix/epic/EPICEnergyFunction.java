/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix.epic;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.structure.Molecule;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * An energy function representing a sum of EPIC terms
 * it evaluates them using the DOF values in curDOFVals
 * 
 * @author mhall44
 */
public class EPICEnergyFunction implements EnergyFunction.NeedsInit, EnergyFunction.DecomposableByDof {
    
    private static final long serialVersionUID = -6797584391009212938L;

    DoubleMatrix1D curDOFVals = null;//this needs to be assigned to something
    //(e.g. curDOFVals in a MolecEObjFunction) that will be appropriately adjusted
    
    boolean useSharedMolec = true;//Evaluate any SAPE terms on a shared molecule for all terms
    //(expected to be much faster in terms of doing less conformation setting, so true by default)
    
    boolean includeMinE = false;//by default, just evaluating the continuous part (not the minE's)
    
    ArrayList<EPoly> terms;//the EPIC terms to evaluate
    
    ArrayList<ArrayList<Integer>> termDOFs;//for each term, which degrees of freedom (in curDOFVals) it operates on

    
    public EPICEnergyFunction(ArrayList<EPoly> terms, boolean includeMinE) {
        //create an energy function from some terms, will assign curDOFVals, termDOFs, and (if needed) sharedMolec
        //later
        this.terms = terms;
        this.includeMinE = includeMinE;
    }
    
    
    
    @Override
    public void init(Molecule molec, List<DegreeOfFreedom> DOFs, DoubleMatrix1D DOFVector) {
        //To use this energy function, something (e.g., a MolecEObjFunction)
        //will set the values of curDOFVals and sharedMolec,
        //and then call this energy function on them
        //set these fields up for use
        
        curDOFVals = DOFVector;
        
        //to evaluate EPIC terms based on curDOFVals, need to know what degrees of freedom (from this EPIC matrix's conformational space)
        //are represented by curDOFVals.
        //We now match these DOFs to the DOFs of each individual EPIC term
        termDOFs = new ArrayList<>();
        for(EPoly term : terms){
            ArrayList<Integer> singleTermDOFs = new ArrayList<>();
            
            for(String dofName : term.DOFNames){//DOF referred to by the current polynomial term
                
                int indexOfDOF = -1;//index of dof in DOFs and thus in curDOFVals
                for(int DOFNum=0; DOFNum<DOFs.size(); DOFNum++){
                    //if(dof == DOFs.get(DOFNum)){
                    if(dofName.equalsIgnoreCase(DOFs.get(DOFNum).getName())){
                        indexOfDOF = DOFNum;
                        break;
                    }
                }
                
                if(indexOfDOF==-1){//DOF in the current polynomial term not found in DOFs, so can't interpret term
                    throw new RuntimeException("ERROR: Degree of freedom in EPoly term not found in conf reference DOFs");
                }
                
                singleTermDOFs.add(indexOfDOF);
            }
            
            termDOFs.add(singleTermDOFs);
        }
        
        
        
        //all terms with SAPE need to know about the shared molecule, if we are using one
        //(molec can be null if we aren't or if we aren't using any SAPE)
        if(useSharedMolec){
            for(EPoly term : terms){
                if(term.sapeTerm != null){
                    term.sapeTerm.assignSharedMolecule(molec);
                }
            }
        }
    }
    
    
    public void unassignSharedMolec(){
        //delete the sharedMolecEnergyFunction for each of the SAPE terms used by this EPICEnergyFunction
        //If a new molecule is made for every minimization and this isn't done,
        //there's effectively a memory leak where different SAPE terms start filling up with old molecules
        for(EPoly term : terms){
            if(term.sapeTerm != null){
                term.sapeTerm.sharedMolecEnergyFunction = null;
            }
        }
    }
    
    
    
    @Override
    public double getEnergy() {
        if(curDOFVals==null){
            throw new RuntimeException("ERROR: Trying to evaluate an EPICEnergyFunction "
                    + "before assigning it to a vector of DOF values");
        }
        
        double E = 0;
        for(int termNum=0; termNum<terms.size(); termNum++){
            EPoly term = terms.get(termNum);
            
            DoubleMatrix1D DOFValsForTerm = DoubleFactory1D.dense.make(term.numDOFs);
            for(int DOFCount=0; DOFCount<term.numDOFs; DOFCount++)
                DOFValsForTerm.set( DOFCount, curDOFVals.get(termDOFs.get(termNum).get(DOFCount)) );
            
            double termVal = term.evaluate(DOFValsForTerm, includeMinE, useSharedMolec);
            E += termVal;
        }
        
        return E;
    }
    
    
    
    
    public ArrayList<Double> allTermValues(){
        //values of all epic terms at current curDOFVals
        ArrayList<Double> ans = new ArrayList<>();
        
        if(curDOFVals==null){
            throw new RuntimeException("ERROR: Trying to evaluate an EPICEnergyFunction "
                    + "before assigning it to a vector of DOF values");
        }
        
        for(int termNum=0; termNum<terms.size(); termNum++){
            EPoly term = terms.get(termNum);
            
            DoubleMatrix1D DOFValsForTerm = DoubleFactory1D.dense.make(term.numDOFs);
            for(int DOFCount=0; DOFCount<term.numDOFs; DOFCount++)
                DOFValsForTerm.set( DOFCount, curDOFVals.get(termDOFs.get(termNum).get(DOFCount)) );
            
            double termVal = term.evaluate(DOFValsForTerm, includeMinE, useSharedMolec);
            ans.add(termVal);
        }
        
        return ans;
    }
    
    
    public void printAllTermValues() {
        for(double termVal : allTermValues())
            System.out.println(termVal);
    }
    
    @Override
    public List<EnergyFunction> decomposeByDof(Molecule molec, List<DegreeOfFreedom> DOFs) {
        //make a list of energy functions that only include the terms involving each of the specified DOFs
        //of the specified molecule
        //these are assumed to be the same DOFs and molecule used by this energy function
        //these are suitable for minimizing this energy function with respect to a single DOF
        
        ArrayList<EnergyFunction> ans = new ArrayList<>();
        
        for(DegreeOfFreedom dof : DOFs){
            String dofName = dof.getName();
            ArrayList<EPoly> dofTerms = new ArrayList<>();
            for(EPoly term : terms){
                if(term.DOFNames.contains(dofName)){
                    dofTerms.add(term);
                }
            }
            
            EPICEnergyFunction partial = new EPICEnergyFunction(dofTerms, includeMinE);
            partial.init(molec, DOFs, curDOFVals);
            
            ans.add(partial);
        }
        
        
        return ans;
    }
}
