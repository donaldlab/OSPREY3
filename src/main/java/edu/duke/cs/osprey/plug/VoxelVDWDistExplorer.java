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

package edu.duke.cs.osprey.plug;

import edu.duke.cs.osprey.bbfree.BBFreeBlock;
import edu.duke.cs.osprey.bbfree.BBFreeDOF;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.util.ArrayList;
import org.apache.commons.math3.optim.linear.LinearConstraint;

/**
 *
 * Explore distances between atoms w/ significant VDW
 * 
 * @author mhall44
 */
public class VoxelVDWDistExplorer {
    
    //first test: visualize feasible region wrt chi's, and also CATS dofs
    //do for optimal vox, and also vox w/ known clashes (will have to identify both)
    
    
            
    public static void main(String[] args){
        
        //load data
        args = new String[] {"KStar.cfg","System.cfg","DEE.cfg"};
        ConfigFileParser cfp = ConfigFileParser.makeFromFilePaths(args);
        cfp.loadData();

        Molecule m = new Strand.Builder(PDBIO.readFile("1CC8.GMEC.pdb")).build().mol;

        Residue mainRes = m.getResByPDBResNumber("45");//for convenience
        
        DegreeOfFreedom dof1 = new FreeDihedral(mainRes,0);
        
        ArrayList<Residue> bbFreeRes = new ArrayList<>();
        for(int a=43; a<=47; a++)
            bbFreeRes.add(m.getResByPDBResNumber(String.valueOf(a)));
        BBFreeBlock stupid = new BBFreeBlock(bbFreeRes);
        DegreeOfFreedom dof2 = stupid.getDOFs().get(3);

        ArrayList<Atom[]> interactingAtoms = new ArrayList<>();
       
        addAtomPair(m, interactingAtoms, "45", "HG11", "47", "HG22");
        addAtomPair(m, interactingAtoms, "45", "HG12", "22", "HG22");
        addAtomPair(m, interactingAtoms, "45", "HG13", "38", "HG12");
        addAtomPair(m, interactingAtoms, "45", "HG21", "38", "HG12");
        addAtomPair(m, interactingAtoms, "45", "HG22", "22", "HG22");
        addAtomPair(m, interactingAtoms, "45", "HG23", "11", "HB");
        
        //showGridVDWBoundaries(dof1, dof2, 170., -1., 18., 2., interactingAtoms);//wt
        //showGridVDWBoundaries(dof1, dof2, 51., -1., 18., 2., interactingAtoms);//not good
        
        //dof1 = stupid.getDOFs().get(0);
        //showGridVDWBoundaries(dof1, dof2, -0.5, -0.5, 1., 1., interactingAtoms);
        
        //DIH ONLY
        dof2 = new FreeDihedral(m.getResByPDBResNumber("38"), 0);
        showGridVDWBoundaries(dof1, dof2, 170., -69.8, 18., 18., interactingAtoms);
    }
    
    
    static void addAtomPair(Molecule m, ArrayList<Atom[]> atomList, String resName1, String atomName1,
            String resName2, String atomName2){
        
        atomList.add( new Atom[] {
            m.getResByPDBResNumber(resName1).getAtomByName(atomName1),
            m.getResByPDBResNumber(resName2).getAtomByName(atomName2)
        } );
    }
    
    public static double getVDWRadius(Atom at){
        if(at.isHydrogen()){
            switch(at.bonds.get(0).elementNumber){
                case 7:
                case 8:
                    return 1.;
                default:
                    return 1.17;
                //DEBUG!! aromatic sohuld really be 1 as well
            }
        }
        else if(at.isCarbon()){
            if(at.bonds.size()==3)//DEBUG!! assuming carbonyl
                return 1.65;
            return 1.75;
        }
        
        switch(at.elementNumber){
            case 7:
                return 1.55;
            case 8:
                return 1.4;
            case 15:
            case 16:
                return 1.8;
            default:
                throw new RuntimeException("ERROR: Unknown VDW rad for el number: "+at.elementNumber);
        }
    }
    
    
     
    static void showGridVDWBoundaries(DegreeOfFreedom dof1, DegreeOfFreedom dof2, double lb1, double lb2,
            double range1, double range2, ArrayList<Atom[]> interactingAtoms){
        //For specified VDW interactions, grid scan the DOF's and see if they're too close, too far or neither
        //FOR NOW assuming no self clashes
        
        int gridSize = 80;
        double step1 = range1/(gridSize-1);
        double step2 = range2/(gridSize-1);
        
        int numVDWPairs = interactingAtoms.size();
        
        //For each quantity (inter-atom dist that should be correct, or whatevs)
        //Make a grid by outputting 1 for clash, 0 for good, -1 for too far
        int gridVals[][][] = new int[numVDWPairs+1][gridSize][gridSize];
        
        for(int g=0; g<gridSize; g++){
            for(int h=0; h<gridSize; h++){
                dof1.apply(lb1+step1*g);
                dof2.apply(lb2+step2*h);
                
                boolean ok = true;
                for(int p=0; p<numVDWPairs; p++){
                    Atom[] pair = interactingAtoms.get(p);
                    double dist = VectorAlgebra.distance(pair[0].getCoords(), pair[1].getCoords());
                    double targetDist = getVDWRadius(pair[0]) + getVDWRadius(pair[1]);
                    
                    if(dist>targetDist+0.25){//too far
                        //ok = false;//DEBUG!!
                        gridVals[p][g][h] = -1;
                    }
                    else if(dist<targetDist-0.4){//too close
                        ok = false;
                        gridVals[p][g][h] = 1;
                    }
                    else//just right
                        gridVals[p][g][h] = 0;
                }
                
                if(ok)//overall local geometry favorability
                    gridVals[numVDWPairs][g][h] = 0;
                else
                    gridVals[numVDWPairs][g][h] = 1;
            }
        }
        
        System.out.println("VDW BOUNDARIES ON GRID: ");
        
        for(int p=0; p<numVDWPairs+1; p++){
            
            if(p<numVDWPairs)
                System.out.println("VDW PAIR "+p+": ");
            else
                System.out.println("OVERALL FAVORABILITY: ");
            
            for(int g=0; g<gridSize; g++){
                for(int h=0; h<gridSize; h++)
                    System.out.print(gridVals[p][g][h]+" ");
                System.out.println();
            }
        }
        
        
        //Let's check linear constr too
        VoxelVDWListChecker vvlc = new VoxelVDWListChecker();
        vvlc.addDOFInterval(dof1, lb1, range1);
        vvlc.addDOFInterval(dof2, lb2, range2);
        for(Atom[] pair : interactingAtoms)
            vvlc.addAtomPair(pair[0], pair[1]);
        ArrayList<LinearConstraint> polytope = vvlc.calcFeasiblePolytope(null);
        System.out.println("LINEAR-BASED APPROX OF FEASIBLE REGION: ");
        double x[] = new double[2];
        for(int g=0; g<gridSize; g++){
            for(int h=0; h<gridSize; h++){
                x[0] = lb1+step1*g;
                x[1] = lb2+step2*h;
                boolean ok = LPChecks.isPointInPolytope(polytope, x);
                System.out.print(ok ? "0 " : "1 ");
            }
            System.out.println();
        }
    }
    
}
