/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.voxq;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.IdealSeparableReference;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 *
 * Calculates the free energy for a voxel.
 * We set zero free energy to indicate zero energy over the entire voxel--this way,
 * the usual pairwise-minimized lower bound (or any other lower bound on the minimized energy)
 * is also a valid lower bound on the free energy
 * This definition is constant across all voxels with the same number of dimensions (e.g., same sequence
 * AS LONG AS width of each is standardized (e.g. 18 degrees)
 * 
 * @author mhall44
 */
public class VoxelGCalculator {
    
    EnergyMatrix emat;
    EPICMatrix epicMat1, epicMat2;
    
    public VoxelGCalculator(SearchProblem sp){
        if(!sp.useEPIC)
            throw new RuntimeException("ERROR NEED EPIC FOR BAR+LUTE");
        
        epicMat1 = sp.epicMat;
        epicMat2 = (EPICMatrix) ObjectIO.deepCopy(epicMat1);
        emat = sp.emat;
    }
    
    
    public double calcG(int[] assignmentList) {
        
        MoleculeModifierAndScorer mms1 = new MoleculeModifierAndScorer(
                epicMat1.internalEnergyFunction(new RCTuple(assignmentList), true),
                epicMat1.getConfSpace(), new RCTuple(assignmentList) );
        
        CCDMinimizer ccdMin = new CCDMinimizer(mms1,false);
        DoubleMatrix1D center = ccdMin.minimize().dofValues;
        MoleculeModifierAndScorer mms2 = new IdealSeparableReference(
                epicMat2.internalEnergyFunction(new RCTuple(assignmentList), true),
                epicMat2.getConfSpace(), new RCTuple(assignmentList), center );
        
        
        VoxelsDeltaG vdg = new VoxelsDeltaG(mms2,mms1,false);
        double E = vdg.estDeltaG(0.05);
        E += ((IdealSeparableReference)mms2).calcG();

        //NOW SUBTRACT OFF ENERGY FOR CONSTANT ZERO VOXEL
        double voxelVolume = computeVoxelVolume(mms1.getConstraints());
        E += IntraVoxelSampler.RT * Math.log(voxelVolume);
        
        return E;
    }
    
    
    double computeVoxelVolume(DoubleMatrix1D voxelBounds[]){
        double vol = 1;
        for(int dim=0; dim<voxelBounds[0].size(); dim++)
            vol *= voxelBounds[1].get(dim) - voxelBounds[0].get(dim);
        return vol;
    }
    
}
