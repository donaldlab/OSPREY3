/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dof;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.SidechainIdealizer;
import edu.duke.cs.osprey.structure.ConfProblem;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;

/**
 *
 * This degree of freedom specifies the ring pucker of a proline
 * 0 indicates the original pucker, and 1 the "flipped" one
 * 
 * @author mhall44
 */
public class ProlinePucker extends DegreeOfFreedom {
    
	private static final long serialVersionUID = -5326286110623671996L;

	public static boolean UP=true, DOWN=false;
    
    boolean basePucker;//original pucker: param=0 describes this
    boolean curPucker;//pucker that we currently want
    Residue res;
    
    
    ConfProblem puckerProblem = null;//if not null, indicates failure to close ring
    
    
    public ProlinePucker(Residue res){
        this.res = res;
        
        if(res.fullName.toUpperCase().startsWith("PRO"))
            basePucker = computePucker(res);
        else {//pucker will arise by mutation.  Use ideal Pro pucker
            Residue idealPro = EnvironmentVars.resTemplates.getTemplateForMutation("PRO", res, true).templateRes;
            basePucker = computePucker(idealPro);
        }
        
        curPucker = basePucker;
    }
    
    
    public static boolean computePucker(Residue res){
        //Compute the pucker with the current geometry, based on chi2
        //get the coordinates defining chi2
        double CA[] = res.getCoordsByAtomName("CA");
        double CB[] = res.getCoordsByAtomName("CB");
        double CG[] = res.getCoordsByAtomName("CG");
        double CD[] = res.getCoordsByAtomName("CD");

        double chi2 = Protractor.measureDihedral( new double[][] {CA,CB,CG,CD} );

        if( chi2 > 0 )
            return UP;
        else
            return DOWN;
    }
    

    
    @Override
    public void apply(double paramVal) {
        
        if(paramVal==0)//original pucker
            curPucker = basePucker;
        else if(paramVal==1)//flipped pucker
            curPucker = !basePucker;
        else
            throw new RuntimeException("ERROR: Bad parameter value of proline pucker: "+paramVal);
        
        //let's generate the correct geometry for the current pucker.  
        SidechainIdealizer.idealizeSidechain(res);
    }
    
    
    public boolean ringCloseable(){
        //ring can be closed coorectly
        return (puckerProblem==null);
    }
    
    public void setRingCloseable(boolean closeable){
        if( closeable && (puckerProblem!=null) ){//remove problem
            puckerProblem.removeFromRes();
            puckerProblem = null;
        }
        else if( (!closeable) && (puckerProblem==null) ){
            puckerProblem = new ConfProblem(this,res);
        }
        //else stay as-is
    }

    public boolean getCurPucker() {
        return curPucker;
    }

    public boolean getBasePucker() {
        return basePucker;
    }
    
    

    public void setBasePucker(boolean basePucker) {
        this.basePucker = basePucker;
    }

    public void setCurPucker(boolean curPucker) {
        this.curPucker = curPucker;
    }
    
    
    
    @Override
    public Residue getResidue() { return res; }
    
    
    @Override
    public DegreeOfFreedom copy() {
        return new ProlinePucker(res);
    }
    
    @Override
    public void setMolecule(Molecule val) {
        
        // match our residue to the one in the other molecule
        res = val.getResByPDBResNumber(res.getPDBResNumber());
    }
}
