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

import edu.duke.cs.osprey.bbfree.BBFreeBlock;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.minimization.ObjectiveFunction.DofBounds;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

/**
 *
 * StrandFlex representing a range of residues with all dihedral DOFs subject to loop closure
 * //(i.e. a BBFreeBlock)
 * 
 * @author mhall44
 */
public class CATSStrandFlex extends StrandFlex {
    
    BBFreeBlock block;

    public CATSStrandFlex(Strand strand, String firstRes, String lastRes){
        List<Residue> catsRes = strand.mol.getResRangeByPDBResNumber(firstRes, lastRes);
        block = new BBFreeBlock(catsRes);
    }
    
    @Override
    public List<? extends DegreeOfFreedom> makeDofs(Strand strand, Molecule mol) {
        BBFreeBlock molBlock = (BBFreeBlock)block.copyForNewMolecule(mol, new LinkedHashMap<>());
        //can throw away the DOF map
        
        return molBlock.getDOFs();
    }

    
    //For now let's just have CATS do a single backbone voxel
    @Override
    public ObjectiveFunction.DofBounds makeBounds(Strand strand) {
        double freeDOFVoxel[][] = block.getFreeDOFVoxel();
        int numDOFs = freeDOFVoxel[0].length;
        DofBounds ans = new DofBounds(numDOFs);
        for(int d=0; d<numDOFs; d++)
            ans.set(d, freeDOFVoxel[0][d], freeDOFVoxel[1][d]);
        return ans;
    }
    
    @Override
    public ArrayList<HashMap<String, double[]>> listBackboneVoxels(SimpleConfSpace.Position pos) {
            if(doesBlockAffectResidue(pos.resNum))
                return new ArrayList(Arrays.asList(defaultBackboneVoxel(pos.strand)));
            else//No flexibility because pos not affected by these CATS DOFs
                return new ArrayList<>();
    }
    
    
    private boolean doesBlockAffectResidue(String resNum){
        List<Residue> resList = block.getResidues();
        for(Residue res : resList){
            if(res.getPDBResNumber().equalsIgnoreCase(resNum))
                return true;
        }
        return false;
    }
    
    
    
}

