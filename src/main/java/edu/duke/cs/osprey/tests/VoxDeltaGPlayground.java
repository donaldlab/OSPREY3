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

package edu.duke.cs.osprey.tests;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.gmec.GMECFinder;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.IdealSeparableReference;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.voxq.VoxelsDeltaG;

/**
 *
 * @author mhall44
 */
public class VoxDeltaGPlayground {
    
    
    public static void main(String args[]){
        //Trying to compute delta G's with continuous entropy between voxels.     
        //args like for findGMEC (currently set up for default 1CC8 system)
        
        ConfigFileParser cfp = ConfigFileParser.makeFromFilePaths(args);//args are configuration files
        cfp.loadData();
        
        
        //for these three confs, pairwise differences (by alignment) and separable reference G's check out
        //int conf1[] = new int[] {5,7,12,5,1,7,4};
        //int conf1[] = new int[] {5,7,7,5,0,7,4};
        //int conf1[] = new int[] {5,10,5,7,1,7,4};
        
        //for voxdg
        int conf1[] = new int[] {5,7,12,5,0,7,4};
        String epicMatrixName = "1CC8.EPICMAT.dat";

        //for 1CC8.bbfree
        //int conf1[] = new int[] {5,10,4,6,19,7,4};
        //String epicMatrixName = "1CC8.bbfree.nolute.EPICMAT.dat";
        
        
        EPICMatrix epicMat1 = (EPICMatrix) ObjectIO.readObject(epicMatrixName, true);
        EPICMatrix epicMat2 = (EPICMatrix) ObjectIO.readObject(epicMatrixName, true);
                
        /*MoleculeModifierAndScorer mms1 = new MoleculeModifierAndScorer(sp1.fullConfE,
            sp1.confSpace, new RCTuple(conf1) );*/
        
        System.out.println("Testing IVS...");

        
        MoleculeModifierAndScorer mms1 = new MoleculeModifierAndScorer(
                epicMat1.internalEnergyFunction(new RCTuple(conf1), true),
                epicMat1.getConfSpace(), new RCTuple(conf1) );
        
        //for doing difference by BAR between two confs
        /*MoleculeModifierAndScorer mms2 = new MoleculeModifierAndScorer(
                epicMat2.internalEnergyFunction(new RCTuple(conf2)), 
                epicMat2.getConfSpace(), new RCTuple(conf2) );*/
        
        
        //For G calc by separable reference
        CCDMinimizer ccdMin = new CCDMinimizer(mms1,false);
        DoubleMatrix1D center = ccdMin.minimize().dofValues;
        MoleculeModifierAndScorer mms2 = new IdealSeparableReference(
                epicMat2.internalEnergyFunction(new RCTuple(conf1), true),
                epicMat2.getConfSpace(), new RCTuple(conf1), center );
        
        System.out.println("SEP REF G: "+((IdealSeparableReference)mms2).calcG());
        
        
        /*for(int rep=0; rep<3; rep++){
            System.out.println("TESTING SEPARABLE REFERENCE...");
            double centerE = mms1.getValue(center);
            System.out.println("CENTER VALS "+mms1.getValue(center)+" "+mms2.getValue(center));
            System.out.println("CENTER DOF VALS: ");
            for(int dof=0; dof<center.size(); dof++)
                System.out.println(mms2.getValForDOF(dof, center.get(dof)));
            System.out.println("UBPOINTS ALONG AXES: ");
            DoubleMatrix1D ub = mms1.getConstraints()[1];
            for(int dof=0; dof<center.size(); dof++){
                System.out.println(mms2.getValForDOF(dof, ub.get(dof)));
                DoubleMatrix1D gg = center.copy();
                gg.set(dof, ub.get(dof));
                System.out.println(mms2.getValue(gg) - centerE);
                System.out.println(mms1.getValForDOF(dof,ub.get(dof)) - mms1.getValForDOF(dof,center.get(dof)));
            }
            System.out.println("2-POINT DIFFERENCES: ");
            for(int[] a : new int[][] {new int[]{0,3}, new int[]{1,2}}){
                DoubleMatrix1D q1 = center.copy();
                DoubleMatrix1D q2 = center.copy();
                DoubleMatrix1D q12 = center.copy();
                q12.set(a[0], ub.get(a[0]));
                q12.set(a[1], ub.get(a[1]));
                q1.set(a[0], ub.get(a[0]));
                q2.set(a[1], ub.get(a[1]));
                System.out.println(mms2.getValue(q12)-mms2.getValue(q2)-mms2.getValue(q1)+centerE);
            }
        }*/
        
        
        
        VoxelsDeltaG vdg = new VoxelsDeltaG(mms1,mms2,false);//energy alignment false for sep ref, true for differencing
        double dG = vdg.estDeltaG(0.05);
        System.out.println("delta G: "+dG+" from "+vdg.numSamplesNeeded()+" samples");
        dG = vdg.estDeltaG(0.05);
        System.out.println("delta G: "+dG+" from "+vdg.numSamplesNeeded()+" samples");
        dG = vdg.estDeltaG(0.05);
        System.out.println("delta G: "+dG+" from "+vdg.numSamplesNeeded()+" samples");
        
        System.out.println("New normalization");
        vdg = new VoxelsDeltaG(mms1,mms2,false);
        dG = vdg.estDeltaG(0.05);
        System.out.println("delta G: "+dG+" from "+vdg.numSamplesNeeded()+" samples");
        dG = vdg.estDeltaG(0.05);
        System.out.println("delta G: "+dG+" from "+vdg.numSamplesNeeded()+" samples");
        dG = vdg.estDeltaG(0.05);
        System.out.println("delta G: "+dG+" from "+vdg.numSamplesNeeded()+" samples");
        /*IntraVoxelSampler ivs = new IntraVoxelSampler(mms);
        for(int s=0; s<20; s++){
            System.out.println(ivs.nextSample());
        }*/
        System.exit(0);


        //conf2 = new int[] {5,7,12,5,1,7,4};
        //conf2 = new int[] {5,7,7,5,0,7,4};
    }
    
    
}
