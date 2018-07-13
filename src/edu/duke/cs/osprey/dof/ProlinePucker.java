/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.dof;

import edu.duke.cs.osprey.dof.deeper.SidechainIdealizer;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
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
// TODO: is this really a DOF? or a residue conf?
public class ProlinePucker extends DegreeOfFreedom {
    
    private static final long serialVersionUID = -5326286110623671996L;

    @Override
    public String getName() {
        return "PROPUCKER"+res.getPDBResNumber();
    }

    public static enum Direction {
        
    	// NOTE: this enum value order is important
    	// for residue confs to have the same order as older proline/conf space code
        DOWN {
            @Override
            public Direction flip() {
                return UP;
            }
        },
        UP {
            @Override
            public Direction flip() {
                return DOWN;
            }
        };
        
        public abstract Direction flip();
    }
    
    Direction basePucker;//original pucker: param=0 describes this
    Direction curPucker;//pucker that we currently want
	ResidueTemplateLibrary templateLib;
    Residue res;
    
    
    ConfProblem puckerProblem = null;//if not null, indicates failure to close ring
    
    
    public ProlinePucker(ResidueTemplateLibrary templateLib, Residue res){
        this.templateLib = templateLib;
        this.res = res;
        
        if(res.fullName.toUpperCase().startsWith("PRO"))
            basePucker = computePucker(res);
        else {//pucker will arise by mutation.  Use ideal Pro pucker
            Residue idealPro = this.templateLib.getTemplateForMutation("PRO", res).templateRes;
            basePucker = computePucker(idealPro);
        }
        
        curPucker = basePucker;
    }
    
    
    public static Direction computePucker(Residue res){
        //Compute the pucker with the current geometry, based on chi2
        //get the coordinates defining chi2
        double CA[] = res.getCoordsByAtomName("CA");
        double CB[] = res.getCoordsByAtomName("CB");
        double CG[] = res.getCoordsByAtomName("CG");
        double CD[] = res.getCoordsByAtomName("CD");

        double chi2 = Protractor.measureDihedral( new double[][] {CA,CB,CG,CD} );

        if( chi2 > 0 )
            return Direction.UP;
        else
            return Direction.DOWN;
    }
    

    
    @Override
    public void apply(double paramVal) {
        
        if(paramVal==0)//original pucker
            apply(basePucker);
        else if(paramVal==1)//flipped pucker
            apply(basePucker.flip());
        else
            throw new RuntimeException("ERROR: Bad parameter value of proline pucker: "+paramVal);
    }
    
    public void apply(Direction val) {
        curPucker = val;
        
        //let's generate the correct geometry for the current pucker.  
        SidechainIdealizer.idealizeSidechain(templateLib, res);
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

    public Direction getCurPucker() {
        return curPucker;
    }

    public Direction getBasePucker() {
        return basePucker;
    }
    
    

    public void setBasePucker(Direction basePucker) {
        this.basePucker = basePucker;
    }

    public void setCurPucker(Direction curPucker) {
        this.curPucker = curPucker;
    }
    
    
    
    @Override
    public Residue getResidue() { return res; }
    
    
    @Override
    public DegreeOfFreedom copy() {
        return new ProlinePucker(templateLib, res);
    }
    
    @Override
    public void setMolecule(Molecule val) {
        
        // match our residue to the one in the other molecule
        res = val.getResByPDBResNumber(res.getPDBResNumber());
    }
    
    @Override
    public DOFBlock getBlock(){
        return null;
    }
}
