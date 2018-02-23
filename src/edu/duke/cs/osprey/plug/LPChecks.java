/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.plug;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import java.util.ArrayList;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.NoFeasibleSolutionException;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;

/**
 *
 * @author mhall44
 */
public class LPChecks {
    
    public static boolean canAddConstr(LinearConstraint newConstr, ArrayList<LinearConstraint> oldConstr){
        //oldConstr define a non-empty polytope; can we add newConstr without making it empty?
        SimplexSolver ss = new SimplexSolver();
        LinearObjectiveFunction of = constrAsFunction(newConstr);
        PointValuePair ans = ss.optimize(of, new LinearConstraintSet(oldConstr), GoalType.MINIMIZE);
        return ans.getValue()<-1e-8;//making sure there is numerically significant intersection
    }
    
    
    public static LinearObjectiveFunction constrAsFunction(LinearConstraint newConstr){
        //make an objective function that is below 0 if newConstr is satisfied

        if(newConstr.getRelationship()==Relationship.GEQ)//ax>b so want b-ax<0
            return new LinearObjectiveFunction(newConstr.getCoefficients().mapMultiply(-1), newConstr.getValue());
        else//ax<b so want ax-b<0
            return new LinearObjectiveFunction(newConstr.getCoefficients(), -newConstr.getValue());
    }
    
    
    public static boolean isPointInPolytope(ArrayList<LinearConstraint> polytope, double[] point){
        for(LinearConstraint c : polytope){
            if(constrAsFunction(c).value(point) > 1e-8)//give a little buffer
                return false;
        }
        return true;
    }
    
    
    public static boolean polytopeHasFeasiblePt(ArrayList<LinearConstraint> polytope){
        //check feasibility
        //DEBUG!!  may need more effish way to do this than 1-by-1 opt?
        /*ArrayList<LinearConstraint> reconstructed = new ArrayList<>();
        for(LinearConstraint constr : polytope){
            if(canAddConstr(constr,reconstructed))
                reconstructed.add(constr);
            else
                return false;
        }
        return true;*/
        //above version leads to unboundedness.  stupet
        if(polytope.isEmpty())
            return true;
        try {
            getFeasiblePt(polytope);
        }
        catch(NoFeasibleSolutionException e){
            return false;
        }
        return true;
    }
    
    
    public static double[] getFeasiblePt(ArrayList<LinearConstraint> polytope){
        //assuming polytope has a feasible pt, find one
        if(polytope.isEmpty())
            return new double[0];
        SimplexSolver ss = new SimplexSolver();
        LinearObjectiveFunction of = constrAsFunction(polytope.get(0));//obj func doesn't matter
        PointValuePair ans = ss.optimize(of, new LinearConstraintSet(polytope), GoalType.MINIMIZE);
        return ans.getPoint();
    }
    
    
    public static LinearConstraint toLinearConstraint(LinearMultivariateRealFunction f){
        //convert f<=0 constraint to linear contraint
        return new LinearConstraint(f.getQ().toArray(), Relationship.LEQ, -f.getR());
    }
    
    public static LinearMultivariateRealFunction toLinearMultivariateRealFunction(LinearConstraint curConstr){
        return toLinearMultivariateRealFunction(curConstr,0);
    }
    
    public static LinearMultivariateRealFunction toLinearMultivariateRealFunction(LinearConstraint curConstr, double ineqTol){
        //convert curConstr to f such that f+ineqTol<=0 is equivalent to curConstr
        //(i.e. relaxing by ineqTol)
        if(curConstr.getRelationship()==Relationship.GEQ){
            return new LinearMultivariateRealFunction(curConstr.getCoefficients().mapMultiply(-1).toArray(), curConstr.getValue()-ineqTol);
        }
        else if(curConstr.getRelationship()==Relationship.LEQ){
            return new LinearMultivariateRealFunction(curConstr.getCoefficients().toArray(), -curConstr.getValue()-ineqTol);
        }
        else
            throw new RuntimeException("ERROR: Unsupported relationship type");
    }
    
    
    public static DoubleMatrix1D getInteriorPt(LinearMultivariateRealFunction[] polytope){
        ArrayList<LinearConstraint> constr = new ArrayList<>();
        for(LinearMultivariateRealFunction f : polytope)
            constr.add(toLinearConstraint(f));
        return getInteriorPt(constr);
    }
    
    public static DoubleMatrix1D getInteriorPt(ArrayList<LinearConstraint> polytope){
        //Find a point in the interior of a polytope,
        //First attempt:
        //maximizing and minimizing a random function, and taking the average
        //Then if still any active constraints, maximize to get to corner opposite constraint,
        //will 
        //average these corners to get an interior pt--just need to have no constraint active
        //at all the corners
        
        
        if(polytope.isEmpty())
            return DoubleFactory1D.dense.make(0);
        
        /*ArrayList<DoubleMatrix1D> corners = new ArrayList<>();
        SimplexSolver ss = new SimplexSolver();
        LinearConstraintSet polytopeConstr = new LinearConstraintSet(polytope);
        DoubleMatrix1D rand = DoubleFactory1D.dense.random(polytope.get(0).getCoefficients().getDimension());
        LinearObjectiveFunction func = new LinearObjectiveFunction(rand.toArray(), 0);
        PointValuePair max = ss.optimize(func, polytopeConstr, GoalType.MAXIMIZE);
        PointValuePair min = ss.optimize(func, polytopeConstr, GoalType.MINIMIZE);
        
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(max.getPoint());
        ans.assign(DoubleFactory1D.dense.make(min.getPoint()),Functions.plus);
        ans.assign(Functions.mult(0.5));
        return ans;*/
        
        //DEBUG!!  how about this just maximize all the constraints and average
        //slow but should be in the interior
        ArrayList<DoubleMatrix1D> corners = new ArrayList<>();
        SimplexSolver ss = new SimplexSolver();
        LinearConstraintSet polytopeConstr = new LinearConstraintSet(polytope);
        for(LinearConstraint constr : polytope){
            LinearObjectiveFunction func = constrAsFunction(constr);
            PointValuePair otherCorner;
            try{
                otherCorner = ss.optimize(func, polytopeConstr, GoalType.MINIMIZE);
            }
            catch(NoFeasibleSolutionException e){
                return null;//there are no interior points
            }
            corners.add(DoubleFactory1D.dense.make(otherCorner.getPoint()));
        }
        
        //do the average
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(polytope.get(0).getCoefficients().getDimension());
        for(DoubleMatrix1D corner : corners){
            ans.assign(corner,Functions.plus);
        }
        ans.assign(Functions.mult(1./corners.size()));
        return ans;
    }
    
}
