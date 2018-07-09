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

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.TupleEnumerator;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import static edu.duke.cs.osprey.plug.LPChecks.polytopeHasFeasiblePt;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Molecule;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

import edu.duke.cs.osprey.structure.PDBIO;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.Relationship;

/**
 *
 * A matrix of polytopes representing the allowed area
 * According to this matrix, confs for a given RC list will be limited
 * to the intersection of the polytopes for all RC tuples within that list
 * (note: DOF bounds won't be in this list bc too redundant)
 * 
 * @author mhall44
 */

//RC tuple polytope has list of DOFs and then the linear constraints on them

public class PolytopeMatrix extends TupleMatrixGeneric<RCTuplePolytope> {
    
    //DEBUG!!
    public ConfSpace cSpace;//conf space this matrix is for
        
    void writePDBFile(Molecule m, String name){//DEBUG!!!  
        PDBIO.writeFile(m, name);
    }
    
    
    public PolytopeMatrix(ConfSpace confSpace){//empty matrix...
        super(confSpace, Double.POSITIVE_INFINITY, null);
    }
    
    public PolytopeMatrix(SearchProblem sp, boolean doPruning){//actually compute the matrix
        super(sp.confSpace, Double.POSITIVE_INFINITY, null);
        this.cSpace = sp.confSpace;
        int numRCsAtPos[] = cSpace.getNumRCsAtPos();
        
        System.out.println("Building PLUG matrix");
        
        
        
        //DEBUG!!!!
        //new RCPairVDWChecker(cSpace, new RCTuple(3,5), sp.shellResidues).buildPolytope();
        
        for(int pos=0; pos<cSpace.numPos; pos++){
            for(int rc=0; rc<numRCsAtPos[pos]; rc++){
                
                RCTuple single = new RCTuple(pos,rc);
                RCPairVDWChecker rpvc = new RCPairVDWChecker(cSpace, single, sp.shellResidues);
                RCTuplePolytope tope = rpvc.buildPolytope();
                setOneBody(pos, rc, tope);
                if(doPruning && tope==null){//rc impossible
                    sp.pruneMat.setOneBody(pos, rc, true);
                    System.out.println("Geometrically pruned RC "+pos+" "+rc);
                }
                
                for(int pos2=0; pos2<pos; pos2++){
                    for(int rc2=0; rc2<numRCsAtPos[pos2]; rc2++){
                        RCTuple pair = new RCTuple(pos,rc,pos2,rc2);
                        if(sp.pruneMat.getPairwise(pos, rc, pos2, rc2)){//already pruned
                            setPairwise(pos, rc, pos2, rc2, null);
                        }
                        else {
                            RCPairVDWChecker rpvc2 = new RCPairVDWChecker(cSpace, pair, sp.shellResidues);
                            RCTuplePolytope tope2 = rpvc2.buildPolytope();
                            setPairwise(pos, rc, pos2, rc2, tope2);
                            if(doPruning && tope2==null){//rc impossible
                                sp.pruneMat.setPairwise(pos, rc, pos2, rc2, true);
                                System.out.println("Geometrically pruned RC pair "+pos+" "+rc+" , "+pos2+" "+rc2);
                            }
                        }
                    }
                }
            }
        }
    }
    
    
    
    public void doMultiTermPruning(PruningMatrix pruneMat, boolean doTriples){
        //Prune pairs and (if indicated) triples accounting for all the terms within the tuple,
        //not just the interaction term
        //Let's do singles here too so everything incompatible with the PLUG matrix is removed
        System.out.println("Doing full PLUG pruning.");
        int numRCsAtPos[] = cSpace.getNumRCsAtPos();
        
        int numPrunedDirectly = 0;
        
        for(int pos=0; pos<cSpace.numPos; pos++){
            for(int rc=0; rc<numRCsAtPos[pos]; rc++){
                if(isTupleFeasible(new RCTuple(pos,rc))){
                    for(int pos2=0; pos2<pos; pos2++){
                        for(int rc2=0; rc2<numRCsAtPos[pos2]; rc2++){
                            RCTuple pair = new RCTuple(pos,rc,pos2,rc2);
                            if(isTupleFeasible(pair)){
                                //pair is OK, yes, but what about triples?
                                if(doTriples){
                                    for(int pos3=0; pos3<pos2; pos3++){
                                        for(int rc3=0; rc3<numRCsAtPos[pos3]; rc3++){
                                            RCTuple triple = pair.addRC(pos3,rc3);
                                            if(!isTupleFeasible(triple)){
                                                pruneMat.setHigherOrder(triple, true);
                                                numPrunedDirectly++;
                                            }
                                        }
                                    }
                                }
                            }//pair is not OK.  How foolish to have not pruned it just because the polytopes are OK on their own!
                            else{
                                pruneMat.setPairwise(pos, rc, pos2, rc2, true);
                                numPrunedDirectly++;
                            }
                        }
                    }
                }
                else {//single not OK
                    pruneMat.setOneBody(pos, rc, true);
                    numPrunedDirectly++;
                }
            }
        }
        
        System.out.println(numPrunedDirectly+" RCs, pairs, and triples pruned directly");
        
        //What else can be pruned as a result?
        int numPrunedThisRound;
        do {
            numPrunedThisRound = 0;
            TupleEnumerator tupEnum = new TupleEnumerator(pruneMat, null, cSpace.numPos);
            ArrayList<RCTuple> pruneable = tupEnum.enumerateUnprunedTuples(1);
            pruneable.addAll(tupEnum.enumerateUnprunedTuples(2));
            for(RCTuple tup : pruneable){
                for(int pos2=0; pos2<cSpace.numPos; pos2++){
                    if(!tup.pos.contains(pos2)){
                        boolean witnessAvailable = false;
                        for(int rc2=0; rc2<numRCsAtPos[pos2]; rc2++){
                            if(!pruneMat.isPruned(tup.addRC(pos2,rc2))){
                                witnessAvailable = true;
                                break;
                            }
                        }
                        if(!witnessAvailable){
                            numPrunedThisRound++;
                            pruneMat.setTupleValue(tup,true);
                        }
                    }
                }
            }
            System.out.println("Singles/pairs pruned this round because nothing compatible: "+numPrunedThisRound);
        } while(numPrunedThisRound>0);
        
        
        System.out.println("Full PLUG pruning complete");
    }
    
    
    public LinearMultivariateRealFunction[] getJoptPolytope(RCTuple conf, ArrayList<DegreeOfFreedom> DOFs,
            double ineqTol){
        //get the polytope for conf as non-redundant LinearMultivariateRealFunctions
        //relax them all by ineqTol
        final double redundancyTol = 1e-6;//DEBUG!! for checking redundancy of constraints.  From SQPMinimizer
        ArrayList<LinearConstraint> constrList = getFullStericPolytope(conf, DOFs);
        
        ArrayList<LinearMultivariateRealFunction> ajlnv = new ArrayList<>();
        for(LinearConstraint c : constrList){
            LinearMultivariateRealFunction f = LPChecks.toLinearMultivariateRealFunction(c,ineqTol);
            boolean redundant = false;
            for(LinearMultivariateRealFunction g : ajlnv){
                if(Math.abs(f.getR()-g.getR())<redundancyTol){
                    if(Algebra.DEFAULT.norm2(f.getQ().copy().assign(g.getQ(),Functions.minus))<redundancyTol*redundancyTol){
                        redundant = true;
                        break;
                    }
                    //DEBUG!!  This is somewhat redundant with constr cache in getFullStericPolytope
                    //but it will get non-voxel redundant constr (intra constr will appears as (at1,at2) and (at2,at1))
                }
            }
            if(!redundant)
                ajlnv.add(f);
        }
        
        LinearMultivariateRealFunction constr[] = ajlnv.toArray(new LinearMultivariateRealFunction[ajlnv.size()]);
        return constr;
    }
    
    
    public ArrayList<LinearConstraint> getFullStericPolytope(RCTuple conf, ArrayList<DegreeOfFreedom> DOFs){
        //Polytope based on atomic contact constraints
        ArrayList<LinearConstraint> ans = new ArrayList<>();
        
        VoxConstrCache constrCache = new VoxConstrCache();
        
        for(int count=0; count<conf.size(); count++){
            
            RCTuplePolytope singlePolytope = getOneBody(conf.pos.get(count),conf.RCs.get(count));
            if(singlePolytope==null){
                return null;//DEBUG!!  CAN COMMENT OUT TO DE-PRUNE THINGS
            }
            else{
                for(LinearConstraint constr : singlePolytope.expandConstraints(DOFs)){
                    if(!constrCache.checkRedundancy(constr))
                        ans.add(constr);
                }
            }
            
            for(int count2=0; count2<count; count2++){
                RCTuplePolytope pairPolytope = getPairwise(conf.pos.get(count),conf.RCs.get(count),
                        conf.pos.get(count2),conf.RCs.get(count2));
                if(pairPolytope==null){
                    return null;//DEBUG!!  CAN COMMENT OUT TO DE-PRUNE THINGS
                }
                else{
                    for(LinearConstraint constr : pairPolytope.expandConstraints(DOFs)){
                        if(!constrCache.checkRedundancy(constr))
                            ans.add(constr);
                    }
                    //ans.addAll(pairPolytope.expandConstraints(DOFs));
                }
            }
        }
        return ans;
    }
    //Use in isTupleFeasible and in both fitting and getting ready to minimize EPIC
    
    
    public ArrayList<LinearConstraint> getFullPolytope(RCTuple conf){
        //Includes both atomic contact and voxel boundary constraints
        //DEBUG!!!  I think full polytope already has voxel constr??
        LinkedHashMap<DegreeOfFreedom,double[]> DOFBounds = calcDOFBounds(conf);
        ArrayList<DegreeOfFreedom> DOFs = new ArrayList<>(DOFBounds.keySet());

        //Much like VoxelVDWListChecker.voxelPolygon
        ArrayList<LinearConstraint> polytope = getFullStericPolytope(conf, DOFs);
        if(polytope==null)
            return null;
        
        int dof=0;
        for(DegreeOfFreedom curDOF : DOFs){
            double unitVec[] = new double[DOFBounds.size()];
            unitVec[dof] = 1;
            double[] bounds = DOFBounds.get(curDOF);
            polytope.add(new LinearConstraint(unitVec,Relationship.GEQ,bounds[0]));
            polytope.add(new LinearConstraint(unitVec,Relationship.LEQ,bounds[1]));
            dof++;
        }
        
        return polytope;
    }
    
    LinkedHashMap<DegreeOfFreedom,double[]> calcDOFBounds(RCTuple conf){
        LinkedHashMap<DegreeOfFreedom,double[]> DOFBounds = new LinkedHashMap<>();
        
        for(int posCount=0; posCount<conf.pos.size(); posCount++){
            //we may actually need DOF intervals for different RCs to differ here...
            //anyway an RC is just a set of box constr, so we enforce the intersection
            RC curRC = cSpace.posFlex.get(conf.pos.get(posCount)).RCs.get(conf.RCs.get(posCount));
            for(int dofCount=0; dofCount<curRC.DOFs.size(); dofCount++){
                DegreeOfFreedom curDOF = curRC.DOFs.get(dofCount);
                double lb = curRC.DOFmin.get(dofCount);
                double ub = curRC.DOFmax.get(dofCount);

                if(DOFBounds.containsKey(curDOF)){
                    double[] curBounds = DOFBounds.get(curDOF);
                    curBounds[0] = Math.max(lb, curBounds[0]);
                    curBounds[1] = Math.min(ub, curBounds[1]);
                }
                else
                    DOFBounds.put(curDOF, new double[]{lb,ub} );
            }
        }
        
        return DOFBounds;
    }
    
    
    public boolean isTupleFeasible(RCTuple conf){
        ArrayList<LinearConstraint> polytope = getFullPolytope(conf);
        if(polytope==null)//already known to be infeasible
            return false;
        else
            return polytopeHasFeasiblePt(polytope);
    }
    
    
    public boolean isPointFeasible(int[] conf, DoubleMatrix1D pt, ArrayList<DegreeOfFreedom> dofs){
        //Is the point pt in the space of DOFs (assumed to correspond to the RCs in conf)
        //feasible given this PolytopeMatrix?
        HashMap<DegreeOfFreedom,Double> DOFValMapping = new HashMap<>();
        if(pt.size()!=dofs.size())
            throw new RuntimeException("ERROR wrong number of DOFs");
        for(int dof=0; dof<pt.size(); dof++)
            DOFValMapping.put(dofs.get(dof), pt.get(dof));
        
        for(int pos=0; pos<cSpace.numPos; pos++){
            RCTuplePolytope tope = getOneBody(pos,conf[pos]);
            if(!tope.containsPoint(ptForDOFs(DOFValMapping,tope.DOFs)))
                return false;
            for(int pos2=0; pos2<pos; pos2++){
                RCTuplePolytope tope2 = getPairwise(pos,conf[pos],pos2,conf[pos2]);
                if(!tope2.containsPoint(ptForDOFs(DOFValMapping,tope2.DOFs)))
                    return false;
            }
        }
        return true;
    }
    
    private double[] ptForDOFs(HashMap<DegreeOfFreedom,Double> DOFValMapping, ArrayList<DegreeOfFreedom> dofList){
        double[] ans = new double[dofList.size()];
        for(int dof=0; dof<dofList.size(); dof++){
            if(DOFValMapping.containsKey(dofList.get(dof)))
                ans[dof] = DOFValMapping.get(dofList.get(dof));
            else
                System.out.println("DOF UNRECOGNIZED;");
        }
        return ans;
    }
    
    
    public void listClashes(int[] conf, DoubleMatrix1D pt, ArrayList<DegreeOfFreedom> dofs){
        //This version of isPointFeasible prints what atom pairs are clashing in pt
        //(judging by the linearized inequalities to define "clashing" of course)
        HashMap<DegreeOfFreedom,Double> DOFValMapping = new HashMap<>();
        if(pt.size()!=dofs.size())
            throw new RuntimeException("ERROR wrong number of DOFs");
        for(int dof=0; dof<pt.size(); dof++)
            DOFValMapping.put(dofs.get(dof), pt.get(dof));
        
        for(int pos=0; pos<cSpace.numPos; pos++){
            RCTuplePolytope tope = getOneBody(pos,conf[pos]);
            if(tope!=null)
                tope.listClashes(ptForDOFs(DOFValMapping,tope.DOFs));
            for(int pos2=0; pos2<pos; pos2++){
                RCTuplePolytope tope2 = getPairwise(pos,conf[pos],pos2,conf[pos2]);
                if(tope2!=null)
                    tope2.listClashes(ptForDOFs(DOFValMapping,tope2.DOFs));
            }
        }
    }
    
    
    
}
