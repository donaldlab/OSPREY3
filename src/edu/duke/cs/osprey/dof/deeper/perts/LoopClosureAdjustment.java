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

package edu.duke.cs.osprey.dof.deeper.perts;

import edu.duke.cs.osprey.structure.ConfProblem;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.RigidBodyMotion;
import edu.duke.cs.osprey.tools.RotationMatrix;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * This is a discrete perturbation that performs a backbone 
 * conformational change on just 3 residues, preserving bond lengths and angles and omegas
 * The motion is composed of the same 3 rotations as a backrub, but for different angles
 * 
 * @author mhall44
 */
public class LoopClosureAdjustment extends Perturbation {
    
    HashMap<String,RigidBodyMotion[][]> solnCache = new HashMap<>();
    //cache loop closure solutions for different starting-coord values
    
    ArrayList<ConfProblem> problems = new ArrayList<>();
    //if the LCA is currently in an invalid state (solution num >= number of solutions),
    //then the conf problems will be listed here.  If valid state, problems is empty.  

    public LoopClosureAdjustment(ArrayList<Residue> resDirectlyAffected) {
        super(resDirectlyAffected);
    }
    
    
    @Override
    public boolean doPerturbationMotion(double paramVal) {
                
        int solnNum = (int)paramVal;//This is a discrete perturbation taking an integer parameter
        RigidBodyMotion[][] rotations = getRotationsForCurState();
        
        if(solnNum>=rotations.length){
            //paramVal is impossible...there aren't that many solutions
            //declare all directly affected residues to have invalid conf
            //(we will not update their coordinates to perturb at all)
            if(problems.isEmpty()){//currently in a conformationally good state...mark as changed
                for(Residue res : resDirectlyAffected)
                    problems.add( new ConfProblem(this,res) );
            }

            return false;
        }
        
        //If we have a solution, our conf problems are over
        for(ConfProblem prob : problems)
            prob.removeFromRes();
        problems = new ArrayList<>();//remove them all
        
        applyBackrubLikeMotion(rotations[solnNum]);//the motion is the same as a backrub, just for different angles
        return true;
    }
    
    
    private RigidBodyMotion[][] getRotationsForCurState(){
        //The available set of states depends on the coordinates of the starting residue
        //backbones.  Calculate available motions for current coords.  
        //We'll cache these to avoid full tripeptide calculation every time...
        //likely in many cases, the same set of motions will be available
        
        String startCoordHash = hashStartingCoords();//String made from coords of all BB atoms used to calc TC
        
        if(!solnCache.containsKey(startCoordHash)){
            RigidBodyMotion[][] solns = calcSolns();
            solnCache.put(startCoordHash, solns);
        }
        
        return solnCache.get(startCoordHash);
    }
    
    
        
        
    public RigidBodyMotion[][] calcSolns(){
        //Solve the loop closure equations given the current starting BB coordinates
        //Return the options for BB rigid-body motions to perform
        //For each closure solution (first index in answer), we have motions 
        //for first and second peptide planes (second index in answer)
        
        
        //This object stores the bond lengths & angles and omega dihedrals
        //and will be used to compute the alternate conformations
        TripeptideClosure tc = new TripeptideClosure(resDirectlyAffected);


        double r_soln_n[][][] = new double[16][3][3];
        double r_soln_a[][][] = new double[16][3][3];
        double r_soln_c[][][] = new double[16][3][3];

        double firstN[] = resDirectlyAffected.get(0).getCoordsByAtomName("N");//Coordinates of the first residue's N
        double firstCA[] = resDirectlyAffected.get(0).getCoordsByAtomName("CA");//Its CA
        double firstC[] = resDirectlyAffected.get(0).getCoordsByAtomName("C");

        double midN[] = resDirectlyAffected.get(1).getCoordsByAtomName("N");//Starting coordinates of the middle N
        double midCA[] = resDirectlyAffected.get(1).getCoordsByAtomName("CA");
        double midC[] = resDirectlyAffected.get(1).getCoordsByAtomName("C");

        double lastN[] = resDirectlyAffected.get(2).getCoordsByAtomName("N");
        double lastCA[] = resDirectlyAffected.get(2).getCoordsByAtomName("CA");
        double lastC[] = resDirectlyAffected.get(2).getCoordsByAtomName("C");


        int numSoln = tc.solve_3pep_poly( firstN, firstCA, lastCA, lastC, r_soln_n, r_soln_a, r_soln_c);


        RigidBodyMotion solns[][] = new RigidBodyMotion[numSoln][2];

        int unperturbed = -1;//The unperturbed state: will be found by least-squares comparison to the original state
        double lowestSum = Double.POSITIVE_INFINITY;

        
        for(int s=0; s<numSoln; s++){

            //First rotation: Based on midCA and midN
            RotationMatrix rm1 = RotationMatrix.getSuperposingRotMatrix( 
                    VectorAlgebra.subtract( midCA, firstCA ), 
                    VectorAlgebra.subtract( r_soln_a[s][1], firstCA),
                    VectorAlgebra.subtract( midN, firstCA ), 
                    VectorAlgebra.subtract( r_soln_n[s][1], firstCA) );

            //This rotation is about the last CA instead of the first one
            RotationMatrix rm2 = RotationMatrix.getSuperposingRotMatrix( 
                    VectorAlgebra.subtract( midC, lastCA ), 
                    VectorAlgebra.subtract( r_soln_c[s][1], lastCA),
                    VectorAlgebra.subtract( lastN, lastCA ), 
                    VectorAlgebra.subtract( r_soln_n[s][2], lastCA) );

            
            //See if this might be the unperturbed state
            double checkSum = VectorAlgebra.normsq( VectorAlgebra.subtract(firstC, r_soln_c[s][0]) )
                    + VectorAlgebra.normsq( VectorAlgebra.subtract(midN, r_soln_n[s][1]) )
                    + VectorAlgebra.normsq( VectorAlgebra.subtract(midCA, r_soln_a[s][1]) )
                    + VectorAlgebra.normsq( VectorAlgebra.subtract(midC, r_soln_c[s][1]) )
                    + VectorAlgebra.normsq( VectorAlgebra.subtract(lastN, r_soln_n[s][2]) );

            if(checkSum < lowestSum){
                lowestSum = checkSum;
                unperturbed = s;
            }
            
            
            //Now create RigidBodyMotions from the two rotation matrices
            solns[s][0] = new RigidBodyMotion(firstCA,rm1,firstCA);
            solns[s][1] = new RigidBodyMotion(lastCA,rm2,lastCA);
        }

        //Now we will set parameter 0 to indicate the unperturbed state
        //meaning both motions will be the identity
        //There may be another solution in solns[0], which will be moved to the unperturbed soln
        if(numSoln>0)
            solns[unperturbed] = solns[0];
        else//ok there should always be an unperturbed state; this would suggest numerical error
            solns = new RigidBodyMotion[1][2];//we'll put the unperturbed state here
        
        solns[0] = new RigidBodyMotion[2];
        solns[0][0] = new RigidBodyMotion(new double[3], RotationMatrix.identity(), new double[3]);
        solns[0][1] = new RigidBodyMotion(new double[3], RotationMatrix.identity(), new double[3]);
        
        return solns;
    }
    
    
    
    
    
    
    private String hashStartingCoords(){
        //Compute a string from the backbone starting coordinates
        //that we need to compute loops closure solutions
        
        String ans = "";
        for(int resNum=0; resNum<3; resNum++){
            for(String atomName : new String[] {"N","CA","C"}){
                double coords[] = resDirectlyAffected.get(resNum).getCoordsByAtomName(atomName);
                for(int dim=0; dim<3; dim++)
                    ans += "|" + String.valueOf(coords[dim]);
            }
        }
        
        return ans;
    }
    
    
    @Override
    public Perturbation copyForNewMolecule(Molecule mol, PerturbationBlock block){
        LoopClosureAdjustment lca = new LoopClosureAdjustment(Residue.equivalentInMolec(resDirectlyAffected, mol));
        lca.curParamVal = curParamVal;
        lca.indexInBlock = indexInBlock;
        lca.block = block;
        
        for(ConfProblem cp : problems)
            lca.problems.add( new ConfProblem(lca,cp.getBrokenResidue().equivalentInMolec(mol)) );

        solnCache = (HashMap<String,RigidBodyMotion[][]>) ObjectIO.deepCopy(solnCache);
        
        return lca;
    }
    
    
}
