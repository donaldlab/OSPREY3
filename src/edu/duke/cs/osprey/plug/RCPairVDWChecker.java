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

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import static edu.duke.cs.osprey.plug.VoxelVDWDistExplorer.getVDWRadius;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.Residue;
import static edu.duke.cs.osprey.plug.VoxelVDWListChecker.DOFInterval;
import edu.duke.cs.osprey.structure.AtomNeighbors;
import edu.duke.cs.osprey.structure.ProbeAtomNeighbors;
import edu.duke.cs.osprey.tools.VectorAlgebra;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import org.apache.commons.math3.optim.linear.LinearConstraint;

/**
 *
 * This class looks for possible contacts between atoms in a pair of RCs
 * or between an RC and the shell if there's just one RC
 * It makes a polytope representing the geometrically favorable region, or determines there is none
 * For now, assuming a polytope suffices and that pushing any of the good contacts that are possible
 * anywhere in the voxel out of the favorable distance range will create a hole or clash
 * Then higher-order checks on the voxels (up to full confs) can limit themselves to the 
 * intersection of pairwise polytopes, and can check the contacts in each pair to make sure
 * each atom that should have a contact does have one
 * 
 * @author mhall44
 */
public class RCPairVDWChecker {
    
    Residue res1;
    Residue res2 = null;//if null we do res1 intra&shell interactions
    ArrayList<Residue> shellResidues = null;
    ArrayList<DOFInterval> dofIntervals = new ArrayList<>();
    ArrayList<LinearConstraint> voxLinConstr = new ArrayList<>();
    ResidueTypeDOF res1MutDOF=null, res2MutDOF=null;//can leave null if AA type already good
    
    //to be filled in during checkVDW
    ArrayList<Atom[]> interactingAtoms;
    ArrayList<LinearConstraint> feasiblePolytope;
    
    
    public RCPairVDWChecker(Residue res1, RC rc1, Residue res2, RC rc2, ArrayList<Residue> shellResidues){
        this.res1 = res1;
        this.res2 = res2;
        this.shellResidues = shellResidues;
        
        initDOFIntervals(rc1, rc2);
        //assuming AA types already match rc1, rc2!  If not will have null pointer when try to mutate
    }
    
    public RCPairVDWChecker(ConfSpace cSpace, RCTuple tup, ArrayList<Residue> shellResidues){
        int pos1 = tup.pos.get(0);
        res1 = cSpace.posFlex.get(pos1).res;
        RC rc1 = cSpace.posFlex.get(pos1).RCs.get(tup.RCs.get(0));
        res1MutDOF = cSpace.mutDOFs.get(pos1);
        this.shellResidues = shellResidues;
        
        switch(tup.pos.size()){
            case 1:
                initDOFIntervals(rc1,null);
                break;
            case 2:
                int pos2 = tup.pos.get(1);
                res2 = cSpace.posFlex.get(pos2).res;
                RC rc2 = cSpace.posFlex.get(pos2).RCs.get(tup.RCs.get(1));
                res2MutDOF = cSpace.mutDOFs.get(pos2);
                initDOFIntervals(rc1,rc2);
                break;
            default:
                throw new RuntimeException("ERROR: Bad tuple size for RCPairVDWChecker");
        }
    }
        
        
    private void initDOFIntervals(RC rc1, RC rc2){
        
        HashSet<DegreeOfFreedom> res1DOFs = new HashSet<>();
        for(int dof=0; dof<rc1.DOFs.size(); dof++){
            DegreeOfFreedom curDOF = rc1.DOFs.get(dof);
            double lb = rc1.DOFmin.get(dof);
            dofIntervals.add( new DOFInterval(curDOF, lb, rc1.DOFmax.get(dof)-lb) );
            res1DOFs.add(curDOF);
        }
        
        if( ! rc1.AAType.equalsIgnoreCase(res1.template.name) )
            res1MutDOF.mutateTo(rc1.AAType);        
        
        if(rc2!=null){
            for(int dof=0; dof<rc2.DOFs.size(); dof++){
                //DegreeOfFreedom curDOF = rc2.DOFs.get(dof);
                DegreeOfFreedom curDOF = rc2.DOFs.get(dof);
                if( ! res1DOFs.contains(curDOF) ){
                    double lb = rc2.DOFmin.get(dof);
                    dofIntervals.add( new DOFInterval(curDOF, lb, rc2.DOFmax.get(dof)-lb) );
                }
            }

            if( ! rc2.AAType.equalsIgnoreCase(res2.template.name) )
                res2MutDOF.mutateTo(rc2.AAType);
        }
    }

    
    
    public boolean checkVDW(){
        return calcPolytopeConstr()!=null;
    }
    
    public ArrayList<LinearConstraint> calcPolytopeConstr(){
        return calcPolytopeConstr(null);
    }
    
    public ArrayList<LinearConstraint> calcPolytopeConstr(ArrayList<String> atomPairNames){
        //Find possible contacting atoms, see if the contacts work
        //return whether they do or not; next can query about polygon or what atoms have contacts
        
        interactingAtoms = new ArrayList<>();
        //go to center of voxel
        for(DOFInterval di : dofIntervals){
            di.dof.apply(di.lb+di.range/2);
        }
        //figure out how far we can get from the center
        double maxMotion = 3;//DEBUG!!!!  basically want voxelUBMotionFromCenter();
        //it would be very questionable to let this be >3 though
        
        VoxelVDWListChecker distChecker = new VoxelVDWListChecker(dofIntervals,new ArrayList<>(),
                new ArrayList<>(),new ArrayList<>(),voxLinConstr);
        //for checking distances in candidate pairs
        
        for(Atom at1 : res1.atoms){
            double vdwRad1 = getVDWRadius(at1);
            //if at1 bb, calc at1 neighbors to make sure at2 not in them
            
            Iterable<Atom> atoms2;
            if(res2==null)
                atoms2 = new ISAtomsIterable(at1,res1,shellResidues);
            else
                atoms2 = new Res2AtomsIterable(at1,res1,res2);
            
            for(Atom at2 : atoms2){
                double vdwRad2 = getVDWRadius(at2);
                double outerDist = vdwRad1+vdwRad2+0.25;//max dist for VDW contact
                double centerDist = VectorAlgebra.distance(at1.getCoords(), at2.getCoords());
                if(centerDist<=outerDist+maxMotion){//could get close enough...
                    //check if can get close enough, grad method should work ok
                    if(centerDist>outerDist){//skip atom pair if can't
                        if(distChecker.DOFsAtAtomPairDist(new Atom[]{at1,at2},outerDist)==null){
                            continue;
                        }
                    }
                    //ok we can get close enough
                    interactingAtoms.add(new Atom[]{at1,at2});
                }
            }
        }
        
        //OK now check feasibility of the specified interactions
        ArrayList<Residue> knownConfRes = new ArrayList<>();
        knownConfRes.add(res1);
        if(res2==null)
            knownConfRes.addAll(shellResidues);
        else
            knownConfRes.add(res2);
        VoxelVDWListChecker vvlc = new VoxelVDWListChecker(dofIntervals,interactingAtoms,
                knownConfRes,shellResidues,voxLinConstr);
        
        return vvlc.calcFeasiblePolytope(atomPairNames);
        //SHOULD ALSO ID POLYTOPES W/ NO ADDL CONSTR I.E. FARAWAY RES!!
    }
    
    
    public RCTuplePolytope buildPolytope(){
        
        ArrayList<String> atomPairNames = new ArrayList<>();//DEBUG!!!  So can see what's clashing
        
        ArrayList<LinearConstraint> constr = calcPolytopeConstr(atomPairNames);
        ArrayList<DegreeOfFreedom> DOFs = new ArrayList<>();
        for(DOFInterval di : dofIntervals)
            DOFs.add(di.dof);
        
        if(constr==null)
            return null;
        
        RCTuplePolytope ans = new RCTuplePolytope(DOFs,constr);
        ans.atomPairNames = atomPairNames;
        return ans;
    }
    
    
    //method to query whether each atom in a ResContactRequirements gets a contact or not
    //method to return polygon
    
    
    
    
}
