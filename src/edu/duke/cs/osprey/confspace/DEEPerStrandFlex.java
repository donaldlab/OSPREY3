/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/


package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.dof.DOFBlock;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.dof.deeper.perts.Perturbation;
import edu.duke.cs.osprey.dof.deeper.perts.PerturbationBlock;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.structure.Molecule;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

/**
 *
 * StrandFlex representing a perturbation block
 * 
 * @author mhall44
 */
public class DEEPerStrandFlex extends StrandFlex {
    
    DEEPerSettings dset;
    ArrayList<Perturbation> perts;//the perturbations, represented in the strand's molecule
    //will be transferred to voxel-specific molecules for minimization etc
    
    public DEEPerStrandFlex(Strand strand, DEEPerSettings dset){
        //We will start with a DEEPerSettings objects that only has parameters specified,
        //but hasn't loaded the pert file or built perturbations into a molecule
        //(i.e. in the form that is made by the constructor or could easily be made in Python)
        dset.loadPertFile(null);
        perts = dset.makePerturbations(strand.mol);
        this.dset = dset;
    }

    
    @Override
    public List<? extends DegreeOfFreedom> makeDofs(Strand strand, Molecule mol) {
        //copy everything to the new molecule
        
        //keep track of which perturbations and blocks have been copied so far
        LinkedHashMap<DegreeOfFreedom,DegreeOfFreedom> original2CopiedDOF = new LinkedHashMap<>();
        HashMap<DOFBlock,DOFBlock> newBlocks = new HashMap<>();

        
        for(DegreeOfFreedom dof : perts){
            DOFBlock block = dof.getBlock();
            if(!newBlocks.containsKey(block)){//first DOF in its block
                    DOFBlock copiedBlock = block.copyForNewMolecule(mol, original2CopiedDOF);
                    //this will copy all the DOFs in the block at once
                    newBlocks.put(block,copiedBlock);
            }
        }
        
        return new ArrayList(original2CopiedDOF.values());
    }

    
    @Override
    public ObjectiveFunction.DofBounds makeBounds(Strand strand) {
        throw new RuntimeException("ERROR: makeBounds is not well-defined for DEEPerStrandFlex,"
                + " because there may be multiple backbone voxels and DOFs may only affect a few positions");
    }
    
    @Override
    public ArrayList<HashMap<String, double[]>> listBackboneVoxels(SimpleConfSpace.Position pos) {
        ArrayList<ArrayList<double[]>> pertIntervals = dset.getPertIntervals();
        ArrayList<ArrayList<int[]>> pertStates = dset.getPertStates(pos.index);
        
        if(pertStates==null){//no DEEPer flexibility...
            pertStates = new ArrayList<>();
            pertStates.add(null);
        }
        
        ArrayList<HashMap<String, double[]>> ans = new ArrayList<>();
        
        for(ArrayList<int[]> pertState : pertStates){
            HashMap<String,double[]> vox = new HashMap<>();
            
            if(pertState != null) {
                //need to add DEEPer DOFs

                for(int[] singlePertState : pertState){
                    int pertNum = singlePertState[0];//number (in perts) of the perturbation we're adding
                    int pertStateNum = singlePertState[1];//state of this perturbation

                    double[] pertInterval = pertIntervals.get(pertNum).get(pertStateNum);//interval for this perturbation in this state
                    vox.put(perts.get(pertNum).getName(), pertInterval.clone());
                }
            }
            ans.add(vox);
        }
        
        return ans;
    }
    
    
}

